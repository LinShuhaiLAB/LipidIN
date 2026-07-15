"""Convert vendor MS raw formats (Thermo .raw, SCIEX .wiff/.wiff2, Bruker/Agilent .d) to mzML so the
rest of the LipidIN pipeline (which consumes mzML) works unchanged.

把厂商原始质谱格式（Thermo .raw、SCIEX .wiff/.wiff2、Bruker/Agilent .d）转换为 mzML，让后续
LipidIN 流程（以 mzML 为输入）无需改动即可运行。

The conversion is delegated to ProteoWizard's msconvert, which ships the vendor readers on Windows and
handles both Thermo QE .raw (QE240 etc.) and SCIEX .wiff. Thermo files acquired with polarity switching
carry both polarities in one file, so for those we split by polarity ("--filter polarity positive/negative")
into the pos/ and neg/ subfolders the pipeline expects. Files already organized in pos/ or neg/ subfolders
(or single-polarity acquisitions) are converted in place. `ensure_polarity_mzml` is the lazy entry point
called just before annotation: if the polarity folder has no mzML but vendor files exist, it converts them.

转换工作交给 ProteoWizard 的 msconvert 完成——它在 Windows 上自带厂商读取器，能处理 Thermo QE .raw
（含 QE240 等）与 SCIEX .wiff。采用极性切换采集的 Thermo 文件在同一个文件里包含正负两种极性，因此对这
类文件按极性拆分（"--filter polarity positive/negative"）到流程所需的 pos/、neg/ 子目录。已经放在 pos/、
neg/ 子目录里的文件（或单极性采集）则就地转换。`ensure_polarity_mzml` 是注释前调用的惰性入口：若极性
目录下没有 mzML 但存在原始文件，则触发转换。
"""

import os
import re
import sys
import glob
import shutil
import subprocess
from pathlib import Path


def configure_utf8_stdio():
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


configure_utf8_stdio()


# Vendor raw formats msconvert can read on Windows.
# msconvert 在 Windows 上可读取的厂商原始格式。
VENDOR_EXTS = (".raw", ".wiff", ".wiff2", ".d")
# Sidecar files that accompany a .wiff and must NOT be treated as convertible units.
# 伴随 .wiff 的附属文件，不能当作独立可转换单元。
_WIFF_SIDECARS = (".wiff.scan", ".dad.scan", ".dad.sidx", ".timeseries.data")


def find_msconvert(explicit=None):
    """Locate the msconvert executable. Search order: explicit arg, env var, PATH, ProteoWizard installs.

    定位 msconvert 可执行文件。搜索顺序：显式参数、环境变量、PATH、ProteoWizard 安装目录。
    Returns the path string, or None if not found.
    """
    if explicit and os.path.exists(explicit):
        return explicit
    env = os.environ.get("MSCONVERT_EXE")
    if env and os.path.exists(env):
        return env
    on_path = shutil.which("msconvert")
    if on_path:
        return on_path
    # Common per-user ProteoWizard install locations on Windows.
    # Windows 上常见的按用户安装的 ProteoWizard 位置。
    patterns = []
    local = os.environ.get("LOCALAPPDATA")
    if local:
        patterns.append(os.path.join(local, "Apps", "ProteoWizard*", "msconvert.exe"))
        patterns.append(os.path.join(local, "ProteoWizard*", "msconvert.exe"))
    for base in (r"C:\Program Files", r"C:\Program Files (x86)"):
        patterns.append(os.path.join(base, "ProteoWizard*", "msconvert.exe"))
        patterns.append(os.path.join(base, "ProteoWizard*", "*", "msconvert.exe"))
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[-1]  # prefer the highest-sorted (usually newest) version
    return None


def _base_name(path):
    name = os.path.basename(path)
    low = name.lower()
    if low.endswith(".wiff2"):
        return name[:-6]
    for ext in VENDOR_EXTS:
        if low.endswith(ext):
            return name[: -len(ext)]
    return os.path.splitext(name)[0]


def list_vendor_files(folder):
    """Return the canonical vendor files in a folder (one per acquisition, sidecars removed).

    返回目录中的规范原始文件（每次采集一个，剔除附属文件）。
    """
    if not os.path.isdir(folder):
        return []
    found = {}  # base name -> chosen file (prefer .wiff over .wiff2 over .d over .raw)
    priority = {".raw": 3, ".wiff": 3, ".d": 2, ".wiff2": 1}
    for entry in sorted(os.listdir(folder)):
        full = os.path.join(folder, entry)
        low = entry.lower()
        if any(low.endswith(s) for s in _WIFF_SIDECARS):
            continue
        ext = None
        if low.endswith(".wiff2"):
            ext = ".wiff2"
        else:
            for e in (".raw", ".wiff", ".d"):
                if low.endswith(e):
                    ext = e
                    break
        if ext is None:
            continue
        if ext == ".raw" and not os.path.isfile(full):
            continue
        if ext == ".d" and not os.path.isdir(full):
            continue
        base = _base_name(entry)
        if base not in found or priority.get(ext, 0) > priority.get(os.path.splitext(found[base])[1].lower(), 0):
            found[base] = full
    return [found[b] for b in sorted(found)]


def list_mzml_files(folder):
    files = []
    for pattern in ("*.mzML", "*.mzml", "*.MZML"):
        files.extend(glob.glob(os.path.join(folder, pattern)))
    return sorted(set(files))


# --------------------------------------------------------------------------- #
# Polarity detection / 极性检测
# A vendor file may be positive-only, negative-only, or polarity-switching ("both").
# Single-polarity files must be attributed to the correct polarity, not counted under both.
# 一个原始文件可能是仅正、仅负，或极性切换（"both"）。单极性文件应归到正确极性，而不是两个极性都计。
# --------------------------------------------------------------------------- #
_POLARITY_CACHE = {}  # (path, size) -> "pos"/"neg"/"both"/None


def polarity_from_name(name):
    """Guess polarity from a filename token (neg/pos), or None. Cheap, no file open.

    从文件名中的 neg/pos 记号猜测极性；无则 None。廉价，无需打开文件。
    """
    low = str(name).lower()
    if re.search(r"(?:^|[_\-\. ])neg", low) or "negative" in low:
        return "neg"
    if re.search(r"(?:^|[_\-\. ])pos", low) or "positive" in low:
        return "pos"
    return None


def detect_thermo_polarity(path, sample=150, log=print):
    """Sample a Thermo .raw's scan filter strings to decide pos/neg/both. None on any failure.

    采样 Thermo .raw 的扫描过滤字符串判断 正/负/both；失败返回 None。
    """
    try:
        from fisher_py.raw_file_reader import RawFileReaderAdapter
        from fisher_py.data import Device
    except Exception:
        return None
    raw = None
    try:
        raw = RawFileReaderAdapter.file_factory(path)
        raw.select_instrument(Device.MS, 1)
        header = raw.run_header_ex
        first, last = int(header.first_spectrum), int(header.last_spectrum)
        step = max(1, (last - first) // max(1, sample))
        pos = neg = 0
        for sn in range(first, last + 1, step):
            try:
                fs = str(raw.get_scan_event_string_for_scan_number(sn))
            except Exception:
                continue
            if "FTMS +" in fs or " + " in f" {fs} " or "ITMS +" in fs:
                pos += 1
            elif "FTMS -" in fs or " - " in f" {fs} " or "ITMS -" in fs:
                neg += 1
    except Exception:
        return None
    finally:
        if raw is not None:
            try:
                raw.dispose()
            except Exception:
                pass
    if pos and neg:
        return "both"
    if pos:
        return "pos"
    if neg:
        return "neg"
    return None


def file_polarity(path, log=print):
    """Best-effort polarity for one vendor file: filename token first, then (for .raw) scan sampling.

    尽力判断单个原始文件的极性：先看文件名记号，再（对 .raw）采样扫描。结果缓存。
    """
    by_name = polarity_from_name(os.path.basename(path))
    if by_name:
        return by_name
    ext = os.path.splitext(path)[1].lower()
    if ext == ".raw" and os.path.isfile(path):
        try:
            key = (path, os.path.getsize(path))
        except OSError:
            key = (path, 0)
        if key in _POLARITY_CACHE:
            return _POLARITY_CACHE[key]
        pol = detect_thermo_polarity(path, log=log)
        _POLARITY_CACHE[key] = pol
        return pol
    return None


def folder_polarity_map(files, log=print):
    """Map each file -> polarity. For unnamed .raw, detect from the FIRST .raw and reuse (1 open).

    把每个文件映射到极性。对无极性名的 .raw，仅检测第一个并复用（只开一次文件）。
    """
    result = {}
    detected_raw = "__none__"
    for f in files:
        by_name = polarity_from_name(os.path.basename(f))
        if by_name:
            result[f] = by_name
            continue
        if os.path.splitext(f)[1].lower() == ".raw":
            if detected_raw == "__none__":
                detected_raw = file_polarity(f, log=log)
            result[f] = detected_raw
        else:
            result[f] = None
    return result


def list_polarity_inputs(data_root, polarity, log=print):
    """Resolve the input files for one polarity. Returns (kind, files, src_folder).

    解析某极性的输入文件。返回 (kind, files, src_folder)。
    kind in {"mzml","raw","wiff","none"}. Prefers mzML, then vendor files inside the polarity
    subfolder (subfolder implies polarity), then polarity-matched vendor files at the data root.
    """
    folder = os.path.join(data_root, polarity)
    mzml = list_mzml_files(folder)
    if mzml:
        return ("mzml", mzml, folder)
    in_folder = list_vendor_files(folder)
    if in_folder:
        return (_vendor_kind(in_folder), in_folder, folder)
    root = list_vendor_files(data_root)
    if root:
        pol_map = folder_polarity_map(root, log=log)
        matched = [f for f in root if pol_map.get(f) in (polarity, "both")]
        # If polarity could not be determined for ANY root file, fall back to using them for this
        # polarity (unknown single-polarity data). Determinate data is attributed exactly.
        # 若所有根目录文件都无法判定极性，则回退为本极性使用（未知的单极性数据）；可判定的则精确归属。
        if not matched and all(v is None for v in pol_map.values()):
            matched = list(root)
        if matched:
            return (_vendor_kind(matched), matched, data_root)
    return ("none", [], folder)


def _vendor_kind(files):
    for f in files:
        low = f.lower()
        if low.endswith(".raw"):
            return "raw"
        if low.endswith(".wiff") or low.endswith(".wiff2"):
            return "wiff"
        if low.endswith(".d"):
            return "wiff"  # Bruker/Agilent .d also goes through msconvert
    return "raw"


def has_vendor_data(data_root, polarities=("pos", "neg")):
    """True if the data root (or its pos/neg subfolders) contains any vendor raw files.

    数据根目录（或其 pos/neg 子目录）是否含有任何厂商原始文件。
    """
    if list_vendor_files(data_root):
        return True
    for polarity in polarities:
        if list_vendor_files(os.path.join(data_root, polarity)):
            return True
    return False


def convert_file(msconvert, src, out_dir, polarity=None, centroid=False, timeout=1800, log=print):
    """Run msconvert on one vendor file into out_dir as mzML. Returns the output path or None.

    对单个原始文件运行 msconvert，输出 mzML 到 out_dir。返回输出路径，失败则 None。
    polarity: None keeps all scans; 'pos'/'neg' adds a polarity filter (for switching acquisitions).
    centroid: when True, apply vendor peak-picking so the mzML is compact (not giant profile data).
    centroid：为 True 时用厂商峰提取，使 mzML 精简（避免巨大的 profile 数据）。
    """
    os.makedirs(out_dir, exist_ok=True)
    base = _base_name(src)
    out_path = os.path.join(out_dir, base + ".mzML")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        log(f"[convert] mzML already exists, skip: {out_path}\n[转换] mzML 已存在，跳过：{out_path}")
        return out_path
    cmd = [msconvert, src, "--mzML", "--outdir", out_dir, "--outfile", base + ".mzML"]
    if centroid:
        # Vendor centroiding shrinks SCIEX/Thermo profile spectra by ~10-50x; mzml_reader handles centroids.
        # 厂商峰提取可将 profile 谱缩小约 10-50 倍；mzml_reader 能处理 centroid 谱。
        cmd += ["--filter", "peakPicking vendor msLevel=1-"]
    if polarity in ("pos", "positive"):
        cmd += ["--filter", "polarity positive"]
    elif polarity in ("neg", "negative"):
        cmd += ["--filter", "polarity negative"]
    pol_word = polarity or "all"
    log(f"[convert] {os.path.basename(src)} -> {base}.mzML (polarity={pol_word})\n"
        f"[转换] {os.path.basename(src)} -> {base}.mzML（极性={pol_word}）")
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, encoding="utf-8", errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:
        log(f"[convert] msconvert timed out on {src}\n[转换] msconvert 转换超时：{src}")
        return None
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").strip().splitlines()[-6:])
        log(f"[convert] msconvert failed (code {proc.returncode}) on {os.path.basename(src)}: {tail}\n"
            f"[转换] msconvert 失败（代码 {proc.returncode}）：{os.path.basename(src)}")
        return None
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    log(f"[convert] no output produced for {src}\n[转换] 未生成输出：{src}")
    return None


def ensure_polarity_mzml(data_root, polarity, log=print, msconvert=None):
    """Ensure <data_root>/<polarity> holds mzML; convert vendor files lazily if it does not.

    确保 <data_root>/<polarity> 下有 mzML；若没有则惰性转换厂商原始文件。返回该目录下的 mzML 列表。

    Two layouts are supported:
      - vendor files inside the polarity folder (single-polarity acquisitions, e.g. .wiff in neg/)
        -> converted in place, keeping every scan.
      - vendor files at the data-root level (Thermo polarity-switching .raw)
        -> converted with a polarity filter into the pos/ and neg/ subfolders.
    """
    folder = os.path.join(data_root, polarity)
    existing = list_mzml_files(folder)
    if existing:
        return existing

    in_folder = list_vendor_files(folder)
    at_root = list_vendor_files(data_root)
    if not in_folder and not at_root:
        return []  # nothing to convert; caller handles the "no mzML" case

    exe = find_msconvert(msconvert)
    if not exe:
        raise FileNotFoundError(
            "Vendor raw files detected but msconvert (ProteoWizard) was not found. Install ProteoWizard "
            "or set the MSCONVERT_EXE environment variable.\n"
            "检测到厂商原始文件，但未找到 msconvert（ProteoWizard）。请安装 ProteoWizard 或设置环境变量 "
            "MSCONVERT_EXE 指向 msconvert.exe。")

    produced = []
    if in_folder:
        log(f"[{polarity}] Converting {len(in_folder)} vendor file(s) in place -> mzML "
            f"/ 就地转换 {len(in_folder)} 个原始文件为 mzML")
        for src in in_folder:
            out = convert_file(exe, src, folder, polarity=None, log=log)
            if out:
                produced.append(out)
    else:
        log(f"[{polarity}] Converting {len(at_root)} vendor file(s) with polarity split -> {folder} "
            f"/ 按极性拆分转换 {len(at_root)} 个原始文件到 {folder}")
        for src in at_root:
            out = convert_file(exe, src, folder, polarity=polarity, log=log)
            if out:
                produced.append(out)
    return list_mzml_files(folder) or produced


def prepare_data_root(data_root, polarities=("pos", "neg"), log=print, msconvert=None):
    """Eagerly convert all vendor files for the given polarities. Returns {polarity: [mzml,...]}.

    对给定极性一次性转换全部原始文件。返回 {极性: [mzml,...]}。供 agent 的“导入/转换”命令使用。
    """
    out = {}
    for polarity in polarities:
        try:
            out[polarity] = ensure_polarity_mzml(data_root, polarity, log=log, msconvert=msconvert)
        except FileNotFoundError as e:
            log(str(e))
            out[polarity] = []
    return out


WIFF_MZML_SUBDIR = "_wiff_mzml"  # compact centroided mzML kept next to the pkls for SCIEX inputs


def build_pkl_for_polarity(data_root, polarity, pkl_dir, MS2_filter, max_workers=None,
                           log=print, msconvert=None):
    """Build the standard pkl files for one polarity, choosing the cheapest route per input format.

    为某极性构建标准 pkl，按输入格式选择最省的路径。返回 {kind, files, sample_names, quant_dir}。
      - mzML  -> MzmlToPklAgent (peak-picks profile itself)
      - .raw  -> tm_raw2pkl DIRECT to pkl (+ MS1 centroid cache); NO intermediate mzML
      - .wiff -> msconvert to a compact centroided mzML (kept for quantify), then MzmlToPklAgent
    """
    kind, files, src_folder = list_polarity_inputs(data_root, polarity, log=log)
    if kind == "none" or not files:
        return {"kind": "none", "files": [], "sample_names": [], "quant_dir": src_folder}
    os.makedirs(pkl_dir, exist_ok=True)
    sample_names = [(os.path.splitext(os.path.basename(f))[0] if kind == "mzml" else _base_name(f))
                    for f in files]

    if kind == "mzml":
        from mzml_reader import MzmlToPklAgent
        MzmlToPklAgent(folder_path=src_folder, MS2_filter=MS2_filter, max_workers=max_workers,
                       recursive=0, output_folder=pkl_dir).run()
        return {"kind": kind, "files": files, "sample_names": sample_names, "quant_dir": src_folder}

    if kind == "raw":
        import tm_raw2pkl
        log(f"[{polarity}] Converting {len(files)} Thermo .raw file(s) DIRECTLY to pkl (no mzML) "
            f"/ 直接把 {len(files)} 个 Thermo .raw 转为 pkl（不生成 mzML）")
        for f in files:
            log(f"[{polarity}] raw->pkl: {os.path.basename(f)} / raw 转 pkl：{os.path.basename(f)}")
            tm_raw2pkl.convert_one_file_to_standard_pkl(
                f, os.path.dirname(f), MS2_filter, output_folder=pkl_dir, polarity=polarity)
        return {"kind": kind, "files": files, "sample_names": sample_names, "quant_dir": src_folder}

    # kind == "wiff" (also Bruker/Agilent .d): must go through msconvert; keep only compact centroid mzML.
    exe = find_msconvert(msconvert)
    if not exe:
        raise FileNotFoundError(
            "SCIEX .wiff detected but msconvert (ProteoWizard) was not found. Install ProteoWizard or "
            "set MSCONVERT_EXE.\n检测到 SCIEX .wiff，但未找到 msconvert（ProteoWizard）。请安装或设置 MSCONVERT_EXE。")
    wiff_dir = os.path.join(pkl_dir, WIFF_MZML_SUBDIR)
    os.makedirs(wiff_dir, exist_ok=True)
    produced = []
    log(f"[{polarity}] Converting {len(files)} SCIEX .wiff to compact centroided mzML "
        f"/ 把 {len(files)} 个 SCIEX .wiff 转为精简 centroid mzML")
    for f in files:
        out = convert_file(exe, f, wiff_dir, polarity=None, centroid=True, log=log)
        if out:
            produced.append(out)
    from mzml_reader import MzmlToPklAgent
    MzmlToPklAgent(folder_path=wiff_dir, MS2_filter=MS2_filter, max_workers=max_workers,
                   recursive=0, output_folder=pkl_dir).run()
    return {"kind": kind, "files": produced, "sample_names": sample_names, "quant_dir": wiff_dir}


def resolve_quant_inputs(data_root, polarity, pkl_dir, log=print):
    """Where quantify should read MS1 for this polarity: returns (quant_dir, quant_files).

    定量时该极性应从哪里读取 MS1：返回 (quant_dir, quant_files)。
    Prefers the compact wiff mzML folder, else the polarity's mzML/raw location.
    """
    wiff_dir = os.path.join(pkl_dir, WIFF_MZML_SUBDIR)
    if os.path.isdir(wiff_dir) and list_mzml_files(wiff_dir):
        return wiff_dir, list_mzml_files(wiff_dir)
    kind, files, src = list_polarity_inputs(data_root, polarity, log=log)
    return src, files


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else r"D:\bio_inf\LipidIN_github\common_demo_wiff"
    pols = sys.argv[2].split(",") if len(sys.argv) > 2 else ["pos", "neg"]
    print("msconvert:", find_msconvert())
    print("vendor at root:", [os.path.basename(f) for f in list_vendor_files(root)])
    for p in pols:
        print(f"vendor in {p}:", [os.path.basename(f) for f in list_vendor_files(os.path.join(root, p))])
    result = prepare_data_root(root, pols)
    for p, files in result.items():
        print(p, "->", len(files), "mzML")
