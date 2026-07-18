import argparse
from pathlib import Path
from EQ2_batchsearchv2 import *
def find_ms1_msms_files(samples_dir):
    samples_dir = Path(samples_dir)
    if not samples_dir.exists():
        raise FileNotFoundError(f"Input samples directory does not exist: {samples_dir}")
    if not samples_dir.is_dir():
        raise NotADirectoryError(f"Input samples path is not a directory: {samples_dir}")
    sample_inputs = sorted(
        str(p) for p in samples_dir.rglob("*MS1_MSMS.pkl*")
        if p.is_file()
    )
    if len(sample_inputs) == 0:
        raise FileNotFoundError(
            f"No files containing MS1_MSMS.pkl were found under: {samples_dir}"
        )
    return sample_inputs
def parse_args():
    parser = argparse.ArgumentParser(
        description="Run EQ3_batch for all MS1_MSMS.pkl files under a sample directory."
    )
    parser.add_argument(
        "--samples_dir",
        required=True,
        help="Directory containing all samples. The script will recursively search for MS1_MSMS.pkl."
    )
    parser.add_argument(
        "--library",
        required=True,
        help="Path to the MS library pkl file."
    )
    parser.add_argument(
        "--cpp_file",
        default="",
        help="Path to search_engine.cpp. Default is empty, consistent with the original code."
    )
    parser.add_argument(
        "--ppm_err1",
        type=float,
        default=5,
        help="ppm_err1_input. Default: 5."
    )
    parser.add_argument(
        "--da_err1",
        type=float,
        default=0.05,
        help="da_err1_input. Default: 0.05."
    )
    parser.add_argument(
        "--ppm_err2",
        type=float,
        default=20,
        help="ppm_err2_input. Default: 20."
    )
    parser.add_argument(
        "--da_err2",
        type=float,
        default=0.02,
        help="da_err2_input. Default: 0.02."
    )
    return parser.parse_args()
if __name__ == "__main__":
    sample_inputs = find_ms1_msms_files(
        r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos\pkl"
    )
    print(f"找到 {len(sample_inputs)} 个 MS1_MSMS.pkl 文件")
    results = EQ3_batch(
        cpp_file_input=r"D:\bio_inf\LipidIN2.0_SZ\code\eq3_similarity.cpp",
        ppm_err1_input=10,
        da_err1_input=0.01,
        ppm_err2_input=20,
        da_err2_input=0.05,
        sample_inputs=sample_inputs,
        library_input=r"D:\bio_inf\02 LipidIN 3.0\pos_common.pkl"
    )