import os
import json
import math
import pickle
import traceback
from bisect import bisect_left
from multiprocessing import cpu_count
import concurrent.futures
import numpy as np
from tqdm import tqdm
C13_MASS_DIFF = 1.0033548378
SUPPORTED_DIRECT_EXTENSIONS = [".raw"]

MS1_CENTROID_CACHE_SUFFIX = "_MS1_centroid.npz"
THERMO_RAW_MAGIC = b"\x01\xa1F\x00i\x00n\x00n\x00i\x00g\x00a\x00n\x00"
STANDARD_PKL_SKILL = {
    "skill_name": "thermo_raw_to_standard_pkl",
    "version": "1.0",
    "scope": "将严格 DDA 的 Thermo .raw 文件转换为两个标准 pickle 文件",
    "input_formats": ["Thermo RAW"],
    "unsupported_acquisition": ["DIA", "SWATH", "AIF", "SRM", "SIM", "PASEF", "缺少母离子 m/z 的 MS2"],
    "ms1_msms_pool_schema": [
        "scan_id：整数",
        "precursor_mz：浮点数，保留 4 位小数",
        "rt：浮点数，单位分钟，保留 2 位小数",
        "fragment_mz_array：浮点数列表，保留 6 位小数",
        "fragment_intensity_array：浮点数列表，归一化后保留 6 位小数",
    ],
    "ms1_only_pool_schema": {
        "ms1_scan_index": "整数",
        "rt": "浮点数，单位分钟，保留 2 位小数",
        "monoisotopic_peaks": [{"mz": "浮点数，保留 6 位小数", "intensity": "浮点数"}],
    },
}
def log_info(message):
    print(f"[INFO] {message}")
def log_warning(message):
    print(f"[WARN] {message}")
def log_error(message):
    print(f"[ERROR] {message}")
def normalize_extension(file_path):
    return os.path.splitext(file_path)[1].lower()
def is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except Exception:
        return False
def assert_thermo_raw_container(file_path):
    if os.path.isdir(file_path):
        raise RuntimeError(f"该 .raw 是文件夹，疑似 Waters 数据，本脚本只支持 Thermo .raw 文件 {file_path}")
    if not os.path.isfile(file_path):
        raise RuntimeError(f"文件不存在 {file_path}")
    with open(file_path, "rb") as f:
        header = f.read(len(THERMO_RAW_MAGIC))
    if header != THERMO_RAW_MAGIC:
        raise RuntimeError(f"文件头不是 Thermo Finnigan RAW 魔数，无法解析 {file_path}")
def import_raw_backend():
    try:
        from fisher_py.raw_file_reader import RawFileReaderAdapter
        from fisher_py.data import Device
        from fisher_py.data.business import Scan
    except Exception as e:
        raise RuntimeError(
            "无法导入 fisher_py，请先执行 pip install fisher_py。"
            f" 导入错误：{type(e).__name__}: {e}"
        )
    return RawFileReaderAdapter, Device, Scan
class ThermoRawReader:
    def __init__(self, file_path):
        self.file_path = file_path
        self._raw = None
        self._scan_cls = None
        self.first_scan = None
        self.last_scan = None
        self._mono_trailer_hint = None
    def __enter__(self):
        assert_thermo_raw_container(self.file_path)
        RawFileReaderAdapter, Device, Scan = import_raw_backend()
        raw = RawFileReaderAdapter.file_factory(self.file_path)
        if raw.is_error:
            raise RuntimeError(f"RawFileReader 打开文件失败 {self.file_path}")
        if raw.in_acquisition:
            raise RuntimeError(f"文件仍在采集中，拒绝解析 {self.file_path}")
        raw.select_instrument(Device.MS, 1)
        header = raw.run_header_ex
        self._raw = raw
        self._scan_cls = Scan
        self.first_scan = int(header.first_spectrum)
        self.last_scan = int(header.last_spectrum)
        return self
    def __exit__(self, exc_type, exc_value, exc_tb):
        if self._raw is not None:
            try:
                self._raw.dispose()
            except Exception:
                pass
            self._raw = None
        return False
    def scan_numbers(self):
        return range(self.first_scan, self.last_scan + 1)
    def get_filter_string(self, scan_number):
        try:
            return str(self._raw.get_scan_event_string_for_scan_number(scan_number))
        except Exception:
            return ""
    def get_ms_order(self, scan_number):
        event = self._raw.get_scan_event_for_scan_number(scan_number)
        return int(event.ms_order.value), event
    def get_polarity(self, scan_number):
        """Return 'pos'/'neg'/None from the scan filter string (+ / -).

        从扫描过滤字符串（+ / -）返回 'pos'/'neg'/None。
        """
        fs = self.get_filter_string(scan_number)
        if "FTMS +" in fs or "ITMS +" in fs or " + " in f" {fs} ":
            return "pos"
        if "FTMS -" in fs or "ITMS -" in fs or " - " in f" {fs} ":
            return "neg"
        return None
    def get_retention_time_minutes(self, scan_number):
        return float(self._raw.retention_time_from_scan_number(scan_number))
    def _get_trailer_map(self, scan_number):
        try:
            trailer = self._raw.get_trailer_extra_information(scan_number)
            return {str(k).strip(): str(v).strip() for k, v in zip(trailer.labels, trailer.values)}
        except Exception:
            return {}
    def _monoisotopic_trailer_indices(self, trailer):
        hint = self._mono_trailer_hint
        if hint is not None:
            labels = trailer.labels
            if all(str(labels[i]).strip() == text for i, text in hint):
                return [i for i, _ in hint]
        
        by_label = {}
        for index, label in enumerate(trailer.labels):
            by_label[str(label).strip()] = index
        found = [(index, text) for text, index in by_label.items() if "monoisotopic" in text.lower()]
        found.sort(key=lambda item: item[0])
        self._mono_trailer_hint = found
        return [i for i, _ in found]

    def get_precursor_mz(self, scan_number, event):
        try:
            trailer = self._raw.get_trailer_extra_information(scan_number)
        except Exception:
            trailer = None
        if trailer is not None:
            values = trailer.values
            for index in self._monoisotopic_trailer_indices(trailer):
                try:
                    mono = float(str(values[index]).strip())
                except Exception:
                    continue
                if math.isfinite(mono) and mono > 0:
                    return mono
        try:
            reaction = event.get_reaction(0)
        except Exception:
            raise ValueError("MS2 扫描没有 reaction 信息，无法读取 precursor m/z")
        precursor = float(reaction.precursor_mass)
        if not math.isfinite(precursor) or precursor <= 0:
            raise ValueError("MS2 扫描的 precursor m/z 无效")
        return precursor
    def get_isolation_width(self, event):
        try:
            return float(event.get_reaction(0).isolation_width)
        except Exception:
            return float("nan")
    def get_peak_arrays(self, scan_number):
        scan = self._scan_cls.from_file(self._raw, scan_number)
        if scan.has_centroid_stream and scan.centroid_scan.length > 0:
            mz_array = np.asarray(scan.centroid_scan.masses, dtype=np.float64)
            intensity_array = np.asarray(scan.centroid_scan.intensities, dtype=np.float64)
            centroided = True
        else:
            segmented = scan.segmented_scan
            mz_array = np.asarray(segmented.positions, dtype=np.float64)
            intensity_array = np.asarray(segmented.intensities, dtype=np.float64)
            centroided = False
        if mz_array.size != intensity_array.size:
            raise ValueError("m/z array 和 intensity array 长度不一致")
        return mz_array, intensity_array, centroided
def infer_acquisition_from_filter(filter_string, ms_order, isolation_width):
    text = f" {filter_string.lower()} "
    if any(term in text for term in [" srm ", " sim ", "q1ms", "q3ms", "cnl", "pr ", " pr@"]):
        return "SRM_OR_SIM"
    if any(term in text for term in ["aif", "bbcid", "sid=", "@cid0.00", "swath"]):
        return "DIA_OR_AIF"
    if ms_order >= 3:
        return "MSN_ABOVE_MS2"
    if ms_order == 2:
        is_dependent = " d " in text
        if is_dependent:
            return "DDA"
        if is_finite_number(isolation_width) and float(isolation_width) >= 4.0:
            return "DIA_OR_AIF"
        return "UNKNOWN_MS2"
    return "MS1"
def detect_acquisition_format_with_reader(reader, max_spectra=500):
    format_counts = {}
    ms1_count = 0
    ms2_count = 0
    precursor_ms2_count = 0
    try:
        for index, scan_number in enumerate(reader.scan_numbers(), start=1):
            if index > max_spectra:
                break
            ms_order, event = reader.get_ms_order(scan_number)
            filter_string = reader.get_filter_string(scan_number)
            isolation_width = reader.get_isolation_width(event) if ms_order == 2 else float("nan")
            if ms_order == 1:
                ms1_count += 1
            if ms_order == 2:
                ms2_count += 1
                try:
                    reader.get_precursor_mz(scan_number, event)
                    precursor_ms2_count += 1
                except Exception:
                    pass
            inferred = infer_acquisition_from_filter(filter_string, ms_order, isolation_width)
            format_counts[inferred] = format_counts.get(inferred, 0) + 1
    except Exception as e:
        return {
            "file_format": "thermo_raw",
            "acquisition_format": "UNKNOWN",
            "format_evidence": {"detection_failed": str(e), "counts": format_counts},
        }
    evidence = {
        "counts": format_counts,
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "precursor_ms2_count": precursor_ms2_count,
    }
    for bad_format in ["SRM_OR_SIM", "DIA_OR_AIF"]:
        if format_counts.get(bad_format, 0) > 0:
            return {"file_format": "thermo_raw", "acquisition_format": bad_format, "format_evidence": evidence}
    if ms2_count > 0 and format_counts.get("DDA", 0) == ms2_count and precursor_ms2_count == ms2_count:
        return {"file_format": "thermo_raw", "acquisition_format": "DDA", "format_evidence": evidence}
    return {"file_format": "thermo_raw", "acquisition_format": "UNKNOWN", "format_evidence": evidence}
def detect_acquisition_format(file_path, max_spectra=500):
    try:
        with ThermoRawReader(file_path) as reader:
            return detect_acquisition_format_with_reader(reader, max_spectra=max_spectra)
    except Exception as e:
        return {
            "file_format": "thermo_raw",
            "acquisition_format": "UNKNOWN",
            "format_evidence": {"detection_failed": str(e), "counts": {}},
        }
def assert_dda_format_info(format_info):
    if format_info.get("acquisition_format") != "DDA":
        raise RuntimeError(
            "当前文件不是严格 DDA 扫描格式，已停止解析。"
            f" acquisition_format={format_info.get('acquisition_format')},"
            f" evidence={format_info.get('format_evidence')}"
        )
    return format_info
def assert_dda_scan_file(file_path):
    return assert_dda_format_info(detect_acquisition_format(file_path))
def find_nearest_peak_index(sorted_mz, target_mz, ppm_tol):
    tolerance = target_mz * ppm_tol * 1e-6
    pos = bisect_left(sorted_mz, target_mz)
    candidates = []
    if pos < len(sorted_mz):
        candidates.append(pos)
    if pos > 0:
        candidates.append(pos - 1)
    best_idx = None
    best_error = None
    for idx in candidates:
        error = abs(sorted_mz[idx] - target_mz)
        if error <= tolerance and (best_error is None or error < best_error):
            best_idx = idx
            best_error = error
    return best_idx
def keep_monoisotopic_peaks(
    mz_array,
    intensity_array,
    isotope_ppm=10,
    max_charge=2,
    max_isotope=3,
    max_isotope_intensity_ratio=1.2,
    min_relative_intensity=0.001,
    max_ms1_peaks=3000,
):
    mz = np.asarray(mz_array, dtype=np.float64).ravel()
    intensity = np.asarray(intensity_array, dtype=np.float64).ravel()
    if mz.size == 0 or intensity.size == 0:
        return []
    valid = np.isfinite(mz) & np.isfinite(intensity) & (intensity > 0)
    if not valid.all():
        if not valid.any():
            return []
        mz = mz[valid]
        intensity = intensity[valid]
    max_intensity = intensity.max()
    if max_intensity <= 0:
        return []
    above_cutoff = intensity > max_intensity * float(min_relative_intensity)
    if not above_cutoff.all():
        if not above_cutoff.any():
            return []
        mz = mz[above_cutoff]
        intensity = intensity[above_cutoff]
    if mz.size > max_ms1_peaks:
        
        brightest = np.argsort(-intensity, kind="stable")[:max_ms1_peaks]
        mz = mz[brightest]
        intensity = intensity[brightest]
    
    if mz.size > 1 and (np.diff(mz) < 0).any():
        by_mz = np.argsort(mz, kind="stable")
        sorted_mz_arr = mz[by_mz]
        sorted_intensity_arr = intensity[by_mz]
    else:
        sorted_mz_arr, sorted_intensity_arr = mz, intensity

    n_peaks = sorted_mz_arr.size
    n_probe = max_charge * max_isotope

    
    isotope_steps = np.arange(1, max_isotope + 1, dtype=np.float64)
    charge_steps = np.arange(1, max_charge + 1, dtype=np.float64)
    spacing = (isotope_steps[None, :] * (C13_MASS_DIFF / charge_steps[:, None])).ravel()

    probe_mz = (sorted_mz_arr[:, None] + spacing[None, :]).ravel()
    pos = np.searchsorted(sorted_mz_arr, probe_mz, side="left")

    right_idx = np.minimum(pos, n_peaks - 1)
    left_idx = pos - 1
    np.maximum(left_idx, 0, out=left_idx)

    right_err = np.abs(sorted_mz_arr[right_idx] - probe_mz)
    left_err = np.abs(sorted_mz_arr[left_idx] - probe_mz)
    tolerance = probe_mz * (isotope_ppm * 1e-6)

    right_ok = (pos < n_peaks) & (right_err <= tolerance)
    left_ok = (pos > 0) & (left_err <= tolerance)

    
    take_left = left_ok & (~right_ok | (left_err < right_err))
    matched = np.where(take_left, left_idx, np.where(right_ok, right_idx, -1))

    found = matched >= 0
    matched_intensity = np.where(found, sorted_intensity_arr[np.maximum(matched, 0)], np.inf)
    ratio_limit = np.repeat(sorted_intensity_arr * float(max_isotope_intensity_ratio), n_probe)
    survives = (found & (matched_intensity <= ratio_limit)).reshape(n_peaks, max_charge, max_isotope)

    
    run_length = np.where(survives.all(axis=2), max_isotope, np.argmin(survives, axis=2))

    matched_flat = matched.tolist()
    run_flat = run_length.ravel().tolist()
    mz_list = sorted_mz_arr.tolist()
    intensity_list = sorted_intensity_arr.tolist()

    is_isotope = bytearray(n_peaks)
    monoisotopic_peaks = []
    append_peak = monoisotopic_peaks.append

    for i in range(n_peaks):
        if is_isotope[i]:
            continue
        probe_base = i * n_probe
        run_base = i * max_charge
        for charge_i in range(max_charge):
            run = run_flat[run_base + charge_i]
            if run:
                start = probe_base + charge_i * max_isotope
                for k in range(start, start + run):
                    is_isotope[matched_flat[k]] = 1
        append_peak({"mz": round(mz_list[i], 6), "intensity": intensity_list[i]})
    return monoisotopic_peaks
def filter_and_normalize_ms2_peaks(mz_array, intensity_array, MS2_filter):
    if len(mz_array) == 0 or len(intensity_array) == 0:
        return None, None
    valid_pairs = []
    for mz, intensity in zip(mz_array, intensity_array):
        try:
            mz_value = float(mz)
            intensity_value = float(intensity)
        except Exception:
            continue
        if math.isfinite(mz_value) and math.isfinite(intensity_value) and intensity_value > 0:
            valid_pairs.append((mz_value, intensity_value))
    if not valid_pairs:
        return None, None
    max_intensity = max(i for _, i in valid_pairs)
    if max_intensity <= 0:
        return None, None
    threshold = float(MS2_filter) * max_intensity
    filtered_peaks = [(mz, intensity) for mz, intensity in valid_pairs if intensity >= threshold]
    if not filtered_peaks:
        return None, None
    filtered_mz_array, filtered_intensity_array = zip(*filtered_peaks)
    total_intensity = float(sum(filtered_intensity_array))
    if total_intensity <= 0:
        return None, None
    return (
        [round(float(mz), 6) for mz in filtered_mz_array],
        [round(float(i) / total_intensity, 6) for i in filtered_intensity_array],
    )
def validate_output_objects(ms1_msms_pool, ms1_only_pool):
    if not isinstance(ms1_msms_pool, list):
        raise RuntimeError("MS1_MSMS 顶层结构不是 list")
    if not isinstance(ms1_only_pool, list):
        raise RuntimeError("MS1_only 顶层结构不是 list")
    for record in ms1_msms_pool:
        if not isinstance(record, list) or len(record) != 5:
            raise RuntimeError("MS1_MSMS 记录结构不符合标准")
        scan_id, precursor_mz, rt, fragment_mz_array, fragment_intensity_array = record
        int(scan_id)
        if not is_finite_number(precursor_mz):
            raise RuntimeError("MS1_MSMS precursor_mz 不是有限数值")
        if not is_finite_number(rt):
            raise RuntimeError("MS1_MSMS rt 不是有限数值")
        if not isinstance(fragment_mz_array, list) or not isinstance(fragment_intensity_array, list):
            raise RuntimeError("MS1_MSMS fragment arrays 必须是 list")
        if len(fragment_mz_array) != len(fragment_intensity_array):
            raise RuntimeError("MS1_MSMS fragment mz 和 intensity 长度不一致")
        if len(fragment_mz_array) == 0:
            raise RuntimeError("MS1_MSMS fragment array 为空")
        try:
            fragment_mz = np.asarray(fragment_mz_array, dtype=np.float64)
        except (TypeError, ValueError):
            raise RuntimeError("MS1_MSMS fragment m/z 中存在非有限数值")
        if not np.isfinite(fragment_mz).all():
            raise RuntimeError("MS1_MSMS fragment m/z 中存在非有限数值")
        try:
            fragment_intensity = np.asarray(fragment_intensity_array, dtype=np.float64)
        except (TypeError, ValueError):
            raise RuntimeError("MS1_MSMS fragment intensity 中存在非有限数值")
        if not np.isfinite(fragment_intensity).all():
            raise RuntimeError("MS1_MSMS fragment intensity 中存在非有限数值")
        intensity_sum = float(fragment_intensity.sum())
        if abs(intensity_sum - 1.0) > 1e-3:
            raise RuntimeError(f"MS1_MSMS fragment intensity 归一化错误，sum={intensity_sum}")
    for record in ms1_only_pool:
        if not isinstance(record, dict):
            raise RuntimeError("MS1_only 记录结构不符合标准")
        for key in ["ms1_scan_index", "rt", "monoisotopic_peaks"]:
            if key not in record:
                raise RuntimeError(f"MS1_only 缺少字段 {key}")
        int(record["ms1_scan_index"])
        if not is_finite_number(record["rt"]):
            raise RuntimeError("MS1_only rt 不是有限数值")
        if not isinstance(record["monoisotopic_peaks"], list):
            raise RuntimeError("MS1_only monoisotopic_peaks 必须是 list")
        
        isfinite = math.isfinite
        for peak in record["monoisotopic_peaks"]:
            if not isinstance(peak, dict):
                raise RuntimeError("MS1_only monoisotopic_peaks 内部元素必须是 dict")
            if "mz" not in peak or "intensity" not in peak:
                raise RuntimeError("MS1_only monoisotopic_peaks 缺少 mz 或 intensity")
            try:
                finite = isfinite(peak["mz"]) and isfinite(peak["intensity"])
            except TypeError:
                finite = False
            if not finite:
                raise RuntimeError("MS1_only monoisotopic_peaks 存在非有限数值")
    return True
def write_pkl_outputs(ms1_msms_pool, ms1_only_pool, ms1_msms_pkl_out, ms1_only_pkl_out):
    os.makedirs(os.path.dirname(ms1_msms_pkl_out), exist_ok=True)
    os.makedirs(os.path.dirname(ms1_only_pkl_out), exist_ok=True)
    with open(ms1_msms_pkl_out, "wb") as fout:
        pickle.dump(ms1_msms_pool, fout, protocol=pickle.HIGHEST_PROTOCOL)
    with open(ms1_only_pkl_out, "wb") as fout:
        pickle.dump(ms1_only_pool, fout, protocol=pickle.HIGHEST_PROTOCOL)
def validate_pkl_outputs(ms1_msms_pkl_out, ms1_only_pkl_out):
    for path, name in ((ms1_msms_pkl_out, "MS1_MSMS"), (ms1_only_pkl_out, "MS1_only")):
        if not os.path.exists(path):
            raise RuntimeError(f"{name} pkl 文件未生成")
        if os.path.getsize(path) == 0:
            raise RuntimeError(f"{name} pkl 文件为空")


def write_ms1_centroid_cache(cache_path, retention_times, mz_chunks, intensity_chunks):

    if not retention_times:
        return None
    counts = np.fromiter((chunk.size for chunk in mz_chunks), dtype=np.int64, count=len(mz_chunks))
    offsets = np.zeros(counts.size + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    mz = np.concatenate(mz_chunks) if mz_chunks else np.empty(0, dtype=np.float64)
    intensity = np.concatenate(intensity_chunks) if intensity_chunks else np.empty(0, dtype=np.float32)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "wb") as fout:
        np.savez(
            fout,
            rt=np.asarray(retention_times, dtype=np.float64),
            offsets=offsets,
            mz=mz.astype(np.float64, copy=False),
            intensity=intensity.astype(np.float32, copy=False),
        )
    return cache_path
def parse_raw_dda_file(
    file_path,
    ms1_msms_pkl_out,
    ms1_only_pkl_out,
    MS2_filter,
    isotope_ppm=10,
    max_charge=2,
    max_isotope=3,
    max_isotope_intensity_ratio=1.2,
    min_ms1_relative_intensity=0.001,
    max_ms1_peaks=3000,
    ms1_centroid_cache_out=None,
    polarity=None,
):
    ms1_msms_pool = []
    ms1_only_pool = []
    total_scan_count = 0
    ms1_count = 0
    ms2_count = 0
    skipped_scan_count = 0
    profile_ms1_count = 0
    ms1_skipped = False
    cache_rt = []
    cache_mz = []
    cache_intensity = []
    
    with ThermoRawReader(file_path) as reader:
        format_info = assert_dda_format_info(detect_acquisition_format_with_reader(reader))
        for scan_number in reader.scan_numbers():
            total_scan_count += 1
            ms_order = None
            try:
                ms_order, event = reader.get_ms_order(scan_number)
                if ms_order not in [1, 2]:
                    skipped_scan_count += 1
                    continue
                # For polarity-switching acquisitions, keep only the requested polarity.
                # 对极性切换采集，只保留请求的极性。
                if polarity is not None and reader.get_polarity(scan_number) not in (polarity, None):
                    skipped_scan_count += 1
                    continue
                rt = reader.get_retention_time_minutes(scan_number)
                mz_array, intensity_array, centroided = reader.get_peak_arrays(scan_number)
                if ms_order == 1:
                    ms1_count += 1
                    if not centroided:
                        profile_ms1_count += 1
                    if ms1_centroid_cache_out is not None:
                        by_mz = np.argsort(mz_array, kind="mergesort")
                        cache_rt.append(float(rt))
                        cache_mz.append(mz_array[by_mz])
                        cache_intensity.append(intensity_array[by_mz].astype(np.float32))
                    monoisotopic_peaks = keep_monoisotopic_peaks(
                        mz_array=mz_array,
                        intensity_array=intensity_array,
                        isotope_ppm=isotope_ppm,
                        max_charge=max_charge,
                        max_isotope=max_isotope,
                        max_isotope_intensity_ratio=max_isotope_intensity_ratio,
                        min_relative_intensity=min_ms1_relative_intensity,
                        max_ms1_peaks=max_ms1_peaks,
                    )
                    ms1_only_pool.append(
                        {
                            "ms1_scan_index": int(ms1_count),
                            "rt": round(float(rt), 2),
                            "monoisotopic_peaks": monoisotopic_peaks,
                        }
                    )
                else:
                    precursor_mz = reader.get_precursor_mz(scan_number, event)
                    ms2_count += 1
                    fragment_mz_array, fragment_intensity_array = filter_and_normalize_ms2_peaks(
                        mz_array, intensity_array, MS2_filter
                    )
                    if fragment_mz_array is None:
                        skipped_scan_count += 1
                        continue
                    ms1_msms_pool.append(
                        [
                            int(ms2_count),
                            round(float(precursor_mz), 4),
                            round(float(rt), 2),
                            fragment_mz_array,
                            fragment_intensity_array,
                        ]
                    )
            except Exception:
                
                
                if ms_order != 2:
                    ms1_skipped = True
                skipped_scan_count += 1
                continue
    validate_output_objects(ms1_msms_pool, ms1_only_pool)
    write_pkl_outputs(ms1_msms_pool, ms1_only_pool, ms1_msms_pkl_out, ms1_only_pkl_out)
    validate_pkl_outputs(ms1_msms_pkl_out, ms1_only_pkl_out)
    ms1_msms_records = len(ms1_msms_pool)
    ms1_only_records = len(ms1_only_pool)
    warnings = []
    if ms1_msms_records == 0:
        warnings.append("MS1_MSMS 输出为空")
    if ms1_only_records == 0:
        warnings.append("MS1_only 输出为空")
    if profile_ms1_count > 0:
        warnings.append(f"{profile_ms1_count} 个 MS1 扫描没有 centroid stream，已回退到 profile 数据")
    cache_written = None
    if ms1_centroid_cache_out is not None and not ms1_skipped:
        cache_written = write_ms1_centroid_cache(ms1_centroid_cache_out, cache_rt, cache_mz, cache_intensity)
    elif ms1_skipped:
        warnings.append("存在读取失败的 MS1 扫描，未写出 MS1 centroid 缓存")
    return {
        "mode": "direct_thermo_raw_dda",
        "file_format": "thermo_raw",
        "acquisition_format": format_info.get("acquisition_format"),
        "total_scan_count": total_scan_count,
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "skipped_scan_count": skipped_scan_count,
        "ms1_msms_records": ms1_msms_records,
        "ms1_only_records": ms1_only_records,
        "ms1_centroid_cache_out": cache_written,
        "warnings": warnings,
    }
def build_output_paths(file_path, root_folder, output_folder=None):
    if output_folder:
        relative_dir = os.path.relpath(os.path.dirname(os.path.abspath(file_path)), os.path.abspath(root_folder))
        out_folder = output_folder if relative_dir == "." else os.path.join(output_folder, relative_dir)
    else:
        out_folder = os.path.dirname(os.path.abspath(file_path))
    os.makedirs(out_folder, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    ms1_msms_output_file = os.path.join(out_folder, f"{base_name}_MS1_MSMS.pkl")
    ms1_only_output_file = os.path.join(out_folder, f"{base_name}_MS1_only.pkl")
    ms1_centroid_cache_file = os.path.join(out_folder, f"{base_name}{MS1_CENTROID_CACHE_SUFFIX}")
    return ms1_msms_output_file, ms1_only_output_file, ms1_centroid_cache_file
def convert_one_file_to_standard_pkl(file_path, root_folder, MS2_filter, output_folder=None, polarity=None):
    if normalize_extension(file_path) not in SUPPORTED_DIRECT_EXTENSIONS:
        raise RuntimeError(f"仅支持 Thermo .raw 文件 {file_path}")
    ms1_msms_output_file, ms1_only_output_file, ms1_centroid_cache_file = build_output_paths(
        file_path, root_folder, output_folder
    )
    result = parse_raw_dda_file(
        file_path=file_path,
        ms1_msms_pkl_out=ms1_msms_output_file,
        ms1_only_pkl_out=ms1_only_output_file,
        MS2_filter=MS2_filter,
        ms1_centroid_cache_out=ms1_centroid_cache_file,
        polarity=polarity,
    )
    result["ms1_msms_pkl_out"] = ms1_msms_output_file
    result["ms1_only_pkl_out"] = ms1_only_output_file
    result["standard_pkl_skill_version"] = STANDARD_PKL_SKILL["version"]
    return result
def collect_candidate_files(folder_path, recursive=True):
    files = []
    if recursive:
        for root, _, file_names in os.walk(folder_path):
            for file_name in file_names:
                file_path = os.path.join(root, file_name)
                if normalize_extension(file_path) in SUPPORTED_DIRECT_EXTENSIONS:
                    files.append(file_path)
    else:
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path) and normalize_extension(file_path) in SUPPORTED_DIRECT_EXTENSIONS:
                files.append(file_path)
    return sorted(files)
def process_single_file(task_args):
    file_path, root_folder, output_folder, MS2_filter = task_args
    file_name = os.path.relpath(file_path, root_folder)
    try:
        result = convert_one_file_to_standard_pkl(
            file_path=file_path,
            root_folder=root_folder,
            MS2_filter=MS2_filter,
            output_folder=output_folder,
        )
        log_info(
            f"File processed {file_name}, MS1_MSMS records={result.get('ms1_msms_records')},"
            f" MS1_only records={result.get('ms1_only_records')}"
            f"\n文件处理完成 {file_name}，MS1_MSMS 记录数={result.get('ms1_msms_records')}，"
            f"MS1_only 记录数={result.get('ms1_only_records')}"
        )
        return {"file": file_path, "status": "success", "result": result}
    except Exception as e:
        error_text = traceback.format_exc(limit=10)
        log_error(f"Failed to process file {file_name}, reason={str(e)}\n处理文件失败 {file_name}，原因={str(e)}")
        return {"file": file_path, "status": "failed", "error": error_text}
def process_raw_files_batch(folder_path, MS2_filter=0.01, max_workers=None, recursive=True, output_folder=None):
    if not os.path.isdir(folder_path):
        log_error(f"Input path is not a valid folder {folder_path}\n输入路径不是有效文件夹 {folder_path}")
        return []
    if output_folder is not None:
        os.makedirs(output_folder, exist_ok=True)
    candidate_files = collect_candidate_files(folder_path, recursive=recursive)
    if not candidate_files:
        log_warning(f"No .raw files found in folder {folder_path}\n文件夹中未找到 .raw 文件 {folder_path}")
        return []
    if max_workers is None:
        max_workers = max(1, int(cpu_count() * 0.85))
    max_workers = max(1, min(int(max_workers), len(candidate_files)))
    log_info(
        f"Found {len(candidate_files)} Thermo .raw files, using {max_workers} processes"
        f"\n发现 {len(candidate_files)} 个 Thermo .raw 文件，使用 {max_workers} 个进程"
    )
    args_list = [(file_path, folder_path, output_folder, MS2_filter) for file_path in candidate_files]
    summary = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_file, task_args) for task_args in args_list]
        with tqdm(total=len(futures), desc="Processing / 处理进度") as pbar:
            for future in concurrent.futures.as_completed(futures):
                try:
                    summary.append(future.result())
                except Exception as e:
                    summary.append({"file": None, "status": "worker_failed", "error": str(e)})
                    log_error(f"Worker process failed, reason={str(e)}\n进程执行失败，原因={str(e)}")
                pbar.update(1)
    summary_path = os.path.join(output_folder if output_folder else folder_path, "raw_conversion_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    success_count = sum(1 for item in summary if item.get("status") == "success")
    failed_count = sum(1 for item in summary if item.get("status") in ("failed", "worker_failed"))
    log_info(
        f"All files processed, success={success_count}, failed={failed_count}, summary={summary_path}"
        f"\n全部文件处理完成，成功={success_count}，失败={failed_count}，汇总={summary_path}"
    )
    return summary
class ThermoRawToPklAgent:
    def __init__(self, folder_path, MS2_filter=0.01, max_workers=1, recursive=True, output_folder=None):
        self.folder_path = folder_path
        self.MS2_filter = MS2_filter
        self.max_workers = max_workers
        self.recursive = bool(recursive)
        self.output_folder = output_folder
    def run(self, folder_path=None):
        return process_raw_files_batch(
            folder_path=folder_path if folder_path else self.folder_path,
            MS2_filter=self.MS2_filter,
            max_workers=self.max_workers,
            recursive=self.recursive,
            output_folder=self.output_folder,
        )
if __name__ == "__main__":
    agent = ThermoRawToPklAgent(
        folder_path=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos",
        MS2_filter=0.001,
        max_workers=8,
        recursive=0,
        output_folder=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos\pkl",
    )
    agent.run()