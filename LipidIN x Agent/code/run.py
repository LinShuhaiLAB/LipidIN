import os
import glob
import time
from pathlib import Path

from mzml_reader import MzmlToPklAgent
from EQ2 import *
from LCI2 import *
from master_table import *
from quantify import *
from qc_analysis import run_qc
# ============================================================
# Only edit these three lines
# 只需修改这三行
# ============================================================
POLARITY = "pos"  # "pos" or "neg" / "pos" 或 "neg"
DATA_ROOT = r"D:\bio_inf\LipidIN_github\common_demo"
LIBRARY_DIR = r"D:\bio_inf\LipidIN_github\library"
EQ3_WORKERS = 1  # depends on your RAM, 16GB -> 1 worker / 看你的计算机运行内存，16GB对应1个线程
MS2_filter = 0.0005  # filter low-abundance MS2 fragments / 过滤低丰度二级碎片
ppm_err1_input = 10  # MS1 ppm error / MS1 ppm误差
da_err1_input = 0.01  # MS1 Dalton error / MS1 道尔顿误差
ppm_err2_input = 25  # MS2 ppm error / MS2 ppm误差
da_err2_input = 0.05  # MS2 Dalton error / MS2 道尔顿误差
max_workers = None  # workers for the other steps / 其他步骤的线程数
# ============================================================
INPUT_EXTENSION = ".mzML"


def list_mzml_files(folder):
    files = []
    for pattern in ("*.mzML", "*.mzml", "*.MZML"):
        files.extend(glob.glob(os.path.join(folder, pattern)))
    return sorted(set(files))


def check_inputs():
    if POLARITY not in ("pos", "neg"):
        raise ValueError(
            f"POLARITY must be 'pos' or 'neg', currently {POLARITY!r}"
            f"\nPOLARITY 只能是 pos 或 neg，当前为 {POLARITY!r}"
        )
    if not os.path.isdir(input_path):
        raise FileNotFoundError(
            f"Input data directory does not exist {input_path}"
            f"\n输入数据目录不存在 {input_path}"
        )
    if not os.path.exists(sep_library):
        raise FileNotFoundError(
            f"The {POLARITY} spectral library does not exist {sep_library}."
            f" Please put {POLARITY}_common.pkl into {LIBRARY_DIR}"
            f"\n{POLARITY} 模式的谱图库不存在 {sep_library}。"
            f" 请把 {POLARITY}_common.pkl 放到 {LIBRARY_DIR}"
        )
    if not os.path.exists(cpp_file):
        raise FileNotFoundError(
            f"The C++ similarity source file does not exist {cpp_file}"
            f"\nC++ 相似度源文件不存在 {cpp_file}"
        )
    mzml_files = list_mzml_files(input_path)
    if len(mzml_files) == 0:
        raise FileNotFoundError(
            f"No {INPUT_EXTENSION} files found under {input_path}"
            f"\n{input_path} 下没有 {INPUT_EXTENSION} 文件"
        )
    print(
        f"[INFO] Polarity {POLARITY}, mzML files {len(mzml_files)}, library {os.path.basename(sep_library)}"
        f"\n[信息] 极性 {POLARITY}，mzML 文件 {len(mzml_files)} 个，库 {os.path.basename(sep_library)}"
    )
    return mzml_files


CODE_DIR = Path(__file__).resolve().parent
input_path = os.path.join(DATA_ROOT, POLARITY)
output_path = os.path.join(input_path, "pkl")
sep_library = os.path.join(LIBRARY_DIR, f"{POLARITY}_common.pkl")
cpp_file = str(CODE_DIR / "eq3_similarity.cpp")
aligned_table = os.path.join(output_path, "rt_ccs_aligned_table.csv")
rt_checked_table = os.path.join(output_path, "rt_ccs_aligned_table_RT_checked.csv")
master_table_csv = os.path.join(output_path, "master_table.csv")
quant_table_csv = os.path.join(output_path, "master_table_quant.csv")
# QC visualization goes into a new subfolder under DATA_ROOT (e.g. common_demo/QC_pos).
# 质控可视化结果放到 DATA_ROOT 下的新子文件夹（例如 common_demo/QC_pos）。
qc_output_dir = os.path.join(DATA_ROOT, f"QC_{POLARITY}")
start = time.perf_counter()
if __name__ == "__main__":
    mzml_files = check_inputs()
    sample_names = [os.path.splitext(os.path.basename(p))[0] for p in mzml_files]

    # Step 1 | Convert mzML files to the standard pkl outputs.
    # 步骤1 | 把 mzML 文件转换为标准 pkl 输出。
    MzmlToPklAgent(folder_path=input_path, MS2_filter=MS2_filter, max_workers=max_workers,
                   recursive=0, output_folder=output_path).run()

    # Step 2 | MS2 library search (qualitative annotation).
    # 步骤2 | MS2 搜库（定性注释）。
    EQ3_batch(cpp_file_input=cpp_file, ppm_err1_input=ppm_err1_input, da_err1_input=da_err1_input,
              ppm_err2_input=ppm_err2_input, da_err2_input=da_err2_input,
              sample_inputs=find_ms1_msms_files(output_path), library_input=sep_library, max_workers=EQ3_WORKERS)

    # Step 3 | Retention time inference / alignment check.
    # 步骤3 | 保留时间推断/校验。
    retention_time_inference(INPUT_CSV=aligned_table, overwrite=False)

    # Step 4 | Build the master annotation table.
    # 步骤4 | 构建注释总表。
    build_master_table(pkl_dir=output_path, rt_checked_csv=rt_checked_table,
                       library_pkl=sep_library, ppm_err1=ppm_err1_input, da_err1=da_err1_input)

    # Step 5 | Quantify (EIC extraction + peak integration) from the mzML files.
    # 步骤5 | 基于 mzML 文件定量（EIC 提取 + 峰积分）。
    quantify_master_table(master_table_csv=master_table_csv, raw_dir=input_path, max_workers=max_workers)

    # Step 6 | Total peak area correction + quality control (CV / PCA / PLS-DA).
    # 步骤6 | 总峰面积校正 + 质控（CV / PCA / PLS-DA）。
    run_qc(quant_csv=quant_table_csv, output_dir=qc_output_dir,
           sample_names=sample_names, polarity=POLARITY)


end = time.perf_counter()
elapsed = end - start
print(f"[INFO] Total pipeline time: {elapsed:.4f} s\n[信息] 全流程耗时: {elapsed:.4f} 秒")
