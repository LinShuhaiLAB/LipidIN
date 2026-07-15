"""Reusable, parameterized runner for the LipidIN qualitative + quantitative pipeline.

可复用、可参数化的 LipidIN 定性+定量流程执行器。

This wraps the exact steps that run.py performs (mzML -> pkl, EQ3 MS2 search, RT inference,
master table, quantification) into a function that can be called for any polarity and any data
root. The interactive agent (lipid_agent.py) uses it; it also runs standalone.

本模块把 run.py 的各步骤（mzML->pkl、EQ3 二级搜库、RT 推断、注释总表、定量）封装成一个可对任意
极性/数据根目录调用的函数。交互式 agent（lipid_agent.py）会调用它，也可单独运行。
"""

import os
import sys
import glob
import time
from pathlib import Path


def configure_utf8_stdio():
    """Force stdout/stderr to UTF-8 so emoji/Chinese never crash on a GBK (cp936) console.

    把 stdout/stderr 强制设为 UTF-8，避免在 GBK(cp936) 控制台下打印 emoji/中文时报错。
    """
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


configure_utf8_stdio()

from mzml_reader import MzmlToPklAgent
import raw_converter
from EQ2 import find_ms1_msms_files
from EQ2_batchsearchv2 import EQ3_batch
from LCI2 import retention_time_inference
from master_table import build_master_table
from quantify import quantify_master_table


CODE_DIR = Path(__file__).resolve().parent


DEFAULT_PARAMS = {
    "MS2_filter": 0.0005,
    "ppm_err1": 10,
    "da_err1": 0.01,
    "ppm_err2": 25,
    "da_err2": 0.05,
    "eq3_workers": 1,
    "max_workers": None,
}


def list_mzml_files(folder):
    files = []
    for pattern in ("*.mzML", "*.mzml", "*.MZML"):
        files.extend(glob.glob(os.path.join(folder, pattern)))
    return sorted(set(files))


def polarity_paths(data_root, polarity, library_dir):
    input_path = os.path.join(data_root, polarity)
    output_path = os.path.join(input_path, "pkl")
    return {
        "polarity": polarity,
        "input_path": input_path,
        "output_path": output_path,
        "sep_library": os.path.join(library_dir, f"{polarity}_common.pkl"),
        "cpp_file": str(CODE_DIR / "eq3_similarity.cpp"),
        "aligned_table": os.path.join(output_path, "rt_ccs_aligned_table.csv"),
        "rt_checked_table": os.path.join(output_path, "rt_ccs_aligned_table_RT_checked.csv"),
        "master_table_csv": os.path.join(output_path, "master_table.csv"),
        "quant_table_csv": os.path.join(output_path, "master_table_quant.csv"),
    }


def validate_polarity_inputs(data_root, polarity, library_dir, log=print):
    """Validate that inputs (mzML OR vendor raw: Thermo .raw / SCIEX .wiff / Bruker .d) and the
    library exist for this polarity. Returns (paths, input_files) where input_files may be raw.

    校验该极性的输入（mzML 或厂商原始文件：Thermo .raw / SCIEX .wiff / Bruker .d）与谱库是否存在。
    返回 (paths, input_files)，input_files 可能是原始文件。
    """
    paths = polarity_paths(data_root, polarity, library_dir)
    if not os.path.exists(paths["sep_library"]):
        raise FileNotFoundError(
            f"The {polarity} spectral library does not exist / {polarity} 谱图库不存在 {paths['sep_library']}"
        )
    # Accept mzML in the polarity folder, or vendor raw files (in the folder or, for single-polarity
    # data, at the data-root level attributed to this polarity).
    # 接受极性目录中的 mzML，或厂商原始文件（在该目录，或对单极性数据在数据根目录并归属到本极性）。
    kind, input_files, _src = raw_converter.list_polarity_inputs(data_root, polarity, log=log)
    if not input_files:
        raise FileNotFoundError(
            f"No mzML/.raw/.wiff files for polarity {polarity} under {data_root}"
            f" / 极性 {polarity} 在 {data_root} 下没有 mzML/.raw/.wiff 文件")
    return paths, input_files


def _merge_params(params):
    merged = dict(DEFAULT_PARAMS)
    if params:
        merged.update({k: v for k, v in params.items() if v is not None})
    return merged


def run_annotation_for_polarity(data_root, polarity, library_dir, params=None, log=print):
    """Qualitative annotation only (steps 1-4): mzML->pkl, EQ3 search, RT inference, master table.

    仅定性注释（第 1-4 步）：mzML->pkl、EQ3 搜库、RT 推断、构建注释总表。不做定量。
    """
    merged = _merge_params(params)
    # EQ3 compiles the C++ helper by relative name, so operate from the code directory.
    # EQ3 用相对名编译 C++ 辅助模块，因此在 code 目录下运行。
    previous_cwd = os.getcwd()
    os.chdir(str(CODE_DIR))
    try:
        paths, input_files = validate_polarity_inputs(data_root, polarity, library_dir, log=log)
        sample_names = [raw_converter._base_name(p) for p in input_files]
        log(f"[{polarity}] Annotation start | input files {len(input_files)} -> {sample_names}")
        t0 = time.perf_counter()

        log(f"[{polarity}] Step 1/4 | Convert input (mzML/.raw/.wiff) to standard pkl / 输入(mzML/.raw/.wiff)转标准 pkl")
        build = raw_converter.build_pkl_for_polarity(
            data_root, polarity, paths["output_path"], MS2_filter=merged["MS2_filter"],
            max_workers=merged["max_workers"], log=log)
        log(f"[{polarity}] Step 1/4 | pkl built from {build['kind']} input / 已从 {build['kind']} 输入构建 pkl")

        log(f"[{polarity}] Step 2/4 | EQ3 MS2 library search / EQ3 二级搜库")
        EQ3_batch(
            cpp_file_input=paths["cpp_file"], ppm_err1_input=merged["ppm_err1"], da_err1_input=merged["da_err1"],
            ppm_err2_input=merged["ppm_err2"], da_err2_input=merged["da_err2"],
            sample_inputs=find_ms1_msms_files(paths["output_path"]),
            library_input=paths["sep_library"], max_workers=merged["eq3_workers"],
        )

        log(f"[{polarity}] Step 3/4 | Retention time inference / 保留时间推断")
        retention_time_inference(INPUT_CSV=paths["aligned_table"], overwrite=False)

        log(f"[{polarity}] Step 4/4 | Build master annotation table / 构建注释总表")
        build_master_table(
            pkl_dir=paths["output_path"], rt_checked_csv=paths["rt_checked_table"],
            library_pkl=paths["sep_library"], ppm_err1=merged["ppm_err1"], da_err1=merged["da_err1"],
        )

        log(f"[{polarity}] Annotation done in {time.perf_counter() - t0:.1f} s / 注释完成")
        paths["sample_names"] = sample_names
        return paths
    finally:
        os.chdir(previous_cwd)


def run_quantification_for_polarity(data_root, polarity, library_dir, params=None, log=print):
    """Quantification only (step 5): EIC extraction + integration. Needs the master table first.

    仅定量（第 5 步）：EIC 提取 + 峰积分。需先完成注释生成 master_table.csv。
    """
    merged = _merge_params(params)
    paths = polarity_paths(data_root, polarity, library_dir)
    if not os.path.exists(paths["master_table_csv"]):
        raise FileNotFoundError(
            f"master_table.csv not found for {polarity}; run annotation first"
            f" / {polarity} 尚无注释总表，请先运行注释 {paths['master_table_csv']}")
    previous_cwd = os.getcwd()
    os.chdir(str(CODE_DIR))
    try:
        log(f"[{polarity}] Quantify (EIC + integration) / 定量（EIC + 积分）")
        # Read MS1 from the right place: the compact wiff mzML folder, or the mzML/.raw location.
        # Thermo .raw is read directly (fisher_py) using the MS1 centroid cache written during annotation.
        # 从正确的位置读取 MS1：精简 wiff mzML 目录，或 mzML/.raw 位置。Thermo .raw 直接读取（fisher_py），
        # 使用注释阶段写出的 MS1 centroid 缓存。
        quant_dir, quant_files = raw_converter.resolve_quant_inputs(
            data_root, polarity, paths["output_path"], log=log)
        quantify_master_table(
            master_table_csv=paths["master_table_csv"], raw_dir=quant_dir,
            raw_files=quant_files or None, cache_dir=paths["output_path"],
            max_workers=merged["max_workers"],
        )
        paths["sample_names"] = [raw_converter._base_name(p) for p in quant_files]
        log(f"[{polarity}] Quantification done / 定量完成")
        return paths
    finally:
        os.chdir(previous_cwd)


def run_pipeline_for_polarity(data_root, polarity, library_dir, params=None, log=print):
    """Full annotation + quantification (steps 1-5) for one polarity. Returns the resolved paths.

    对单个极性运行完整的注释 + 定量（第 1-5 步）。返回解析后的路径。
    """
    run_annotation_for_polarity(data_root, polarity, library_dir, params=params, log=log)
    return run_quantification_for_polarity(data_root, polarity, library_dir, params=params, log=log)


if __name__ == "__main__":
    import sys

    data_root = sys.argv[1] if len(sys.argv) > 1 else r"D:\bio_inf\LipidIN_github\common_demo"
    library_dir = sys.argv[2] if len(sys.argv) > 2 else r"D:\bio_inf\LipidIN_github\library"
    polarities = sys.argv[3].split(",") if len(sys.argv) > 3 else ["pos", "neg"]

    for pol in polarities:
        run_pipeline_for_polarity(data_root, pol, library_dir)
