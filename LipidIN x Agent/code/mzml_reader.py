"""Convert mzML files into the same standard pickle outputs the .raw pipeline produces.

将 mzML 文件转换为与 .raw 流程完全一致的标准 pickle 输出。

This module is the mzML counterpart of tm_raw2pkl.py. It is deliberately tolerant of the
very different mzML dialects produced by ProteoWizard/msconvert, Thermo, Waters, Bruker and
Sciex exporters: profile or centroided spectra, retention time in minutes or seconds, missing
optional fields, zlib compressed / 32- or 64-bit binary arrays (decoded transparently by
pyteomics), and precursor m/z stored either as a selected ion or only as an isolation window.

本模块是 tm_raw2pkl.py 的 mzML 版本。它有意兼容各类导出软件（ProteoWizard/msconvert、Thermo、
Waters、Bruker、Sciex）生成的差异极大的 mzML 方言：profile 或 centroid 谱、以分钟或秒为单位的保留
时间、缺失的可选字段、zlib 压缩与 32/64 位二进制数组（由 pyteomics 透明解码）、母离子 m/z 以
selected ion 或仅以隔离窗口形式存储等情况。
"""

import os
import json
import traceback
from multiprocessing import cpu_count
import concurrent.futures

import numpy as np
from tqdm import tqdm

# Reuse the exact peak/monoisotopic/validation/writing helpers from the .raw pipeline so the
# mzML outputs are byte-for-byte compatible with everything downstream.
# 直接复用 .raw 流程里的碎片/单同位素/校验/写盘辅助函数，保证 mzML 输出与下游完全兼容。
from tm_raw2pkl import (
    MS1_CENTROID_CACHE_SUFFIX,
    keep_monoisotopic_peaks,
    filter_and_normalize_ms2_peaks,
    validate_output_objects,
    write_pkl_outputs,
    validate_pkl_outputs,
    write_ms1_centroid_cache,
    build_output_paths,
    log_info,
    log_warning,
    log_error,
)

SUPPORTED_MZML_EXTENSIONS = [".mzml"]


def normalize_extension(file_path):
    return os.path.splitext(file_path)[1].lower()


def retention_time_in_minutes(scan):
    """Read scan start time and normalize it to minutes, honoring the CV unit when present.

    读取扫描开始时间并统一换算为分钟，若存在单位注释则据此换算。
    """
    value = scan.get("scan start time")
    if value is None:
        # Some exporters store it directly on the scan under a different key.
        # 部分导出软件把保留时间放在扫描的其它字段里。
        for key in ("retention time", "scan time"):
            if key in scan:
                value = scan.get(key)
                break
    if value is None:
        return None
    unit = str(getattr(value, "unit_info", "") or "").lower()
    rt = float(value)
    if unit.startswith("second") or unit in ("sec", "s"):
        return rt / 60.0
    # Default assumption is minutes (the most common mzML convention).
    # 默认按分钟处理（这是 mzML 最常见的约定）。
    return rt


def spectrum_is_centroided(spectrum):
    """Decide whether a spectrum is already centroided; default to profile when unstated.

    判断谱图是否已做峰检出；未标注时按 profile 处理（后续自行做峰检出）。
    """
    if "centroid spectrum" in spectrum:
        return True
    if "profile spectrum" in spectrum:
        return False
    return False


def centroid_spectrum(mz_array, intensity_array):
    """Peak-pick a profile spectrum into centroids using local maxima + parabolic apex refinement.

    用局部极大值加抛物线顶点插值，把 profile 谱峰检出为质心峰。
    """
    mz = np.asarray(mz_array, dtype=np.float64).ravel()
    intensity = np.asarray(intensity_array, dtype=np.float64).ravel()
    if mz.size == 0:
        return mz, intensity
    if mz.size < 3:
        keep = intensity > 0
        return mz[keep], intensity[keep]

    if (np.diff(mz) < 0).any():
        order = np.argsort(mz, kind="mergesort")
        mz = mz[order]
        intensity = intensity[order]

    y0 = intensity[:-2]
    y1 = intensity[1:-1]
    y2 = intensity[2:]
    is_peak = (y1 > y0) & (y1 >= y2) & (y1 > 0)
    idx = np.nonzero(is_peak)[0] + 1

    if idx.size == 0:
        # No local maxima found: the data likely was already centroided, keep positive points.
        # 找不到局部极大值：数据很可能本就是质心谱，保留强度大于 0 的点。
        keep = intensity > 0
        return mz[keep], intensity[keep]

    a = intensity[idx - 1]
    b = intensity[idx]
    c = intensity[idx + 1]
    denom = a - 2.0 * b + c

    shift = np.zeros_like(b)
    concave = denom < 0
    shift[concave] = 0.5 * (a[concave] - c[concave]) / denom[concave]
    shift = np.clip(shift, -0.5, 0.5)

    left_spacing = mz[idx] - mz[idx - 1]
    right_spacing = mz[idx + 1] - mz[idx]
    spacing = np.where(shift >= 0, right_spacing, left_spacing)

    peak_mz = mz[idx] + shift * spacing
    peak_intensity = np.where(concave, b - 0.25 * (a - c) * shift, b)
    peak_intensity = np.maximum(peak_intensity, b)
    return peak_mz, peak_intensity


def extract_precursor_mz(spectrum):
    """Pull the precursor m/z from a selected ion, or fall back to the isolation window target.

    从 selected ion 读取母离子 m/z，读不到时回退到隔离窗口中心。
    """
    precursor_list = spectrum.get("precursorList")
    if not precursor_list:
        return None
    precursors = precursor_list.get("precursor") or []
    for precursor in precursors:
        selected_ion_list = precursor.get("selectedIonList") or {}
        for selected_ion in selected_ion_list.get("selectedIon") or []:
            for key in ("selected ion m/z", "m/z", "selected ion mz"):
                if key in selected_ion:
                    try:
                        mz = float(selected_ion[key])
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(mz) and mz > 0:
                        return mz
    # Fallback: use the isolation window target m/z when no selected ion is annotated.
    # 回退：没有 selected ion 时使用隔离窗口的目标 m/z。
    for precursor in precursors:
        window = precursor.get("isolationWindow") or {}
        target = window.get("isolation window target m/z")
        if target is not None:
            try:
                mz = float(target)
            except (TypeError, ValueError):
                continue
            if np.isfinite(mz) and mz > 0:
                return mz
    return None


def infer_ms_level(spectrum):
    level = spectrum.get("ms level")
    if level is not None:
        try:
            return int(level)
        except (TypeError, ValueError):
            pass
    # Infer from spectrum-type cvParams when the ms level cvParam is absent.
    # 当缺少 ms level 参数时，用谱图类型参数推断级别。
    if "MS1 spectrum" in spectrum:
        return 1
    if "MSn spectrum" in spectrum or spectrum.get("precursorList"):
        return 2
    return None


def parse_mzml_file(
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
):
    """Parse one mzML file into MS1_MSMS.pkl, MS1_only.pkl and the MS1 centroid cache.

    将单个 mzML 文件解析为 MS1_MSMS.pkl、MS1_only.pkl 以及 MS1 质心缓存。
    """
    from pyteomics import mzml

    ms1_msms_pool = []
    ms1_only_pool = []
    total_scan_count = 0
    ms1_count = 0
    ms2_count = 0
    ms2_without_precursor = 0
    skipped_scan_count = 0
    unknown_level_count = 0
    missing_rt_count = 0
    cache_rt = []
    cache_mz = []
    cache_intensity = []

    with mzml.MzML(file_path) as reader:
        for spectrum in reader:
            total_scan_count += 1
            try:
                ms_level = infer_ms_level(spectrum)
                if ms_level not in (1, 2):
                    unknown_level_count += 1
                    skipped_scan_count += 1
                    continue

                scan_list = spectrum.get("scanList") or {}
                scans = scan_list.get("scan") or [{}]
                rt = retention_time_in_minutes(scans[0])
                if rt is None:
                    missing_rt_count += 1
                    skipped_scan_count += 1
                    continue

                mz_array = np.asarray(spectrum.get("m/z array", []), dtype=np.float64)
                intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=np.float64)
                if mz_array.size == 0 or intensity_array.size == 0:
                    skipped_scan_count += 1
                    continue

                already_centroided = spectrum_is_centroided(spectrum)

                if ms_level == 1:
                    ms1_count += 1
                    if ms1_centroid_cache_out is not None:
                        by_mz = np.argsort(mz_array, kind="mergesort")
                        cache_rt.append(float(rt))
                        cache_mz.append(mz_array[by_mz])
                        cache_intensity.append(intensity_array[by_mz].astype(np.float32))

                    if already_centroided:
                        peak_mz, peak_intensity = mz_array, intensity_array
                    else:
                        peak_mz, peak_intensity = centroid_spectrum(mz_array, intensity_array)

                    monoisotopic_peaks = keep_monoisotopic_peaks(
                        mz_array=peak_mz,
                        intensity_array=peak_intensity,
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
                    precursor_mz = extract_precursor_mz(spectrum)
                    if precursor_mz is None:
                        ms2_without_precursor += 1
                        skipped_scan_count += 1
                        continue
                    ms2_count += 1

                    if already_centroided:
                        peak_mz, peak_intensity = mz_array, intensity_array
                    else:
                        peak_mz, peak_intensity = centroid_spectrum(mz_array, intensity_array)

                    fragment_mz_array, fragment_intensity_array = filter_and_normalize_ms2_peaks(
                        peak_mz, peak_intensity, MS2_filter
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
                skipped_scan_count += 1
                continue

    validate_output_objects(ms1_msms_pool, ms1_only_pool)
    write_pkl_outputs(ms1_msms_pool, ms1_only_pool, ms1_msms_pkl_out, ms1_only_pkl_out)
    validate_pkl_outputs(ms1_msms_pkl_out, ms1_only_pkl_out)

    warnings = []
    if len(ms1_msms_pool) == 0:
        warnings.append("MS1_MSMS output is empty / MS1_MSMS 输出为空")
    if len(ms1_only_pool) == 0:
        warnings.append("MS1_only output is empty / MS1_only 输出为空")
    if ms2_without_precursor > 0:
        warnings.append(
            f"{ms2_without_precursor} MS2 scans had no precursor m/z and were skipped"
            f" / {ms2_without_precursor} 个 MS2 扫描缺少母离子 m/z，已跳过"
        )
    if unknown_level_count > 0:
        warnings.append(
            f"{unknown_level_count} scans had an unrecognized MS level and were skipped"
            f" / {unknown_level_count} 个扫描 MS 级别无法识别，已跳过"
        )
    if missing_rt_count > 0:
        warnings.append(
            f"{missing_rt_count} scans had no retention time and were skipped"
            f" / {missing_rt_count} 个扫描缺少保留时间，已跳过"
        )

    cache_written = None
    if ms1_centroid_cache_out is not None and cache_rt:
        cache_written = write_ms1_centroid_cache(
            ms1_centroid_cache_out, cache_rt, cache_mz, cache_intensity
        )

    return {
        "mode": "direct_mzml",
        "file_format": "mzml",
        "total_scan_count": total_scan_count,
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "skipped_scan_count": skipped_scan_count,
        "ms1_msms_records": len(ms1_msms_pool),
        "ms1_only_records": len(ms1_only_pool),
        "ms1_centroid_cache_out": cache_written,
        "warnings": warnings,
    }


def convert_one_mzml_to_standard_pkl(file_path, root_folder, MS2_filter, output_folder=None):
    if normalize_extension(file_path) not in SUPPORTED_MZML_EXTENSIONS:
        raise RuntimeError(
            f"Only mzML files are supported {file_path}"
            f"\n仅支持 mzML 文件 {file_path}"
        )
    ms1_msms_output_file, ms1_only_output_file, ms1_centroid_cache_file = build_output_paths(
        file_path, root_folder, output_folder
    )
    result = parse_mzml_file(
        file_path=file_path,
        ms1_msms_pkl_out=ms1_msms_output_file,
        ms1_only_pkl_out=ms1_only_output_file,
        MS2_filter=MS2_filter,
        ms1_centroid_cache_out=ms1_centroid_cache_file,
    )
    result["ms1_msms_pkl_out"] = ms1_msms_output_file
    result["ms1_only_pkl_out"] = ms1_only_output_file
    return result


def collect_candidate_files(folder_path, recursive=True):
    files = []
    if recursive:
        for root, _, file_names in os.walk(folder_path):
            for file_name in file_names:
                file_path = os.path.join(root, file_name)
                if normalize_extension(file_path) in SUPPORTED_MZML_EXTENSIONS:
                    files.append(file_path)
    else:
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path) and normalize_extension(file_path) in SUPPORTED_MZML_EXTENSIONS:
                files.append(file_path)
    return sorted(files)


def process_single_mzml_file(task_args):
    file_path, root_folder, output_folder, MS2_filter = task_args
    file_name = os.path.relpath(file_path, root_folder)
    try:
        result = convert_one_mzml_to_standard_pkl(
            file_path=file_path,
            root_folder=root_folder,
            MS2_filter=MS2_filter,
            output_folder=output_folder,
        )
        log_info(
            f"File converted {file_name}, MS1_MSMS records={result.get('ms1_msms_records')},"
            f" MS1_only records={result.get('ms1_only_records')}"
            f"\n文件转换完成 {file_name}，MS1_MSMS 记录数={result.get('ms1_msms_records')}，"
            f"MS1_only 记录数={result.get('ms1_only_records')}"
        )
        for warning in result.get("warnings", []):
            log_warning(warning)
        return {"file": file_path, "status": "success", "result": result}
    except Exception as e:
        error_text = traceback.format_exc(limit=10)
        log_error(
            f"Failed to convert file {file_name}, reason={str(e)}"
            f"\n处理文件失败 {file_name}，原因={str(e)}"
        )
        return {"file": file_path, "status": "failed", "error": error_text}


def process_mzml_files_batch(folder_path, MS2_filter=0.01, max_workers=None, recursive=True, output_folder=None):
    if not os.path.isdir(folder_path):
        log_error(
            f"Input path is not a valid folder {folder_path}"
            f"\n输入路径不是有效文件夹 {folder_path}"
        )
        return []
    if output_folder is not None:
        os.makedirs(output_folder, exist_ok=True)

    candidate_files = collect_candidate_files(folder_path, recursive=recursive)
    if not candidate_files:
        log_warning(
            f"No mzML files found in folder {folder_path}"
            f"\n文件夹中未找到 mzML 文件 {folder_path}"
        )
        return []

    if max_workers is None:
        max_workers = max(1, int(cpu_count() * 0.85))
    max_workers = max(1, min(int(max_workers), len(candidate_files)))
    log_info(
        f"Found {len(candidate_files)} mzML files, using {max_workers} processes"
        f"\n发现 {len(candidate_files)} 个 mzML 文件，使用 {max_workers} 个进程"
    )

    args_list = [(file_path, folder_path, output_folder, MS2_filter) for file_path in candidate_files]
    summary = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_mzml_file, task_args) for task_args in args_list]
        with tqdm(total=len(futures), desc="Converting mzML / mzML 转换进度") as pbar:
            for future in concurrent.futures.as_completed(futures):
                try:
                    summary.append(future.result())
                except Exception as e:
                    summary.append({"file": None, "status": "worker_failed", "error": str(e)})
                    log_error(
                        f"Worker process failed, reason={str(e)}"
                        f"\n进程执行失败，原因={str(e)}"
                    )
                pbar.update(1)

    summary_path = os.path.join(output_folder if output_folder else folder_path, "mzml_conversion_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    success_count = sum(1 for item in summary if item.get("status") == "success")
    failed_count = sum(1 for item in summary if item.get("status") in ("failed", "worker_failed"))
    log_info(
        f"All files processed, success={success_count}, failed={failed_count}, summary={summary_path}"
        f"\n全部文件处理完成，成功={success_count}，失败={failed_count}，汇总={summary_path}"
    )
    return summary


class MzmlToPklAgent:
    """Drop-in mzML analogue of ThermoRawToPklAgent with the same run() interface.

    与 ThermoRawToPklAgent 接口一致的 mzML 版转换器。
    """

    def __init__(self, folder_path, MS2_filter=0.01, max_workers=1, recursive=True, output_folder=None):
        self.folder_path = folder_path
        self.MS2_filter = MS2_filter
        self.max_workers = max_workers
        self.recursive = bool(recursive)
        self.output_folder = output_folder

    def run(self, folder_path=None):
        return process_mzml_files_batch(
            folder_path=folder_path if folder_path else self.folder_path,
            MS2_filter=self.MS2_filter,
            max_workers=self.max_workers,
            recursive=self.recursive,
            output_folder=self.output_folder,
        )
