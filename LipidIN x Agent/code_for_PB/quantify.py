import os
import glob
import time
import concurrent.futures
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

from tm_raw2pkl import MS1_CENTROID_CACHE_SUFFIX


DROP_COLUMNS = ["fa_composition", "source_raw", "n_merged", "rt_fit_status", "rt_fit_error"]

MZ_TOLERANCE = 0.01
RT_WINDOW = 0.75
APEX_FLOOR_RATIO = 0.05


EVIDENCE_MS2 = "MS2注释"
MISSING_TAG = "缺失"

N_ANNOTATIONS_COLUMN = "n_annotations"
CANDIDATES_COLUMN = "candidates"


def log(message):
    print(f"[信息] {message}")


def trapezoid(y, x):
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def read_ms1_scans(raw_path):
    from fisher_py.raw_file_reader import RawFileReaderAdapter
    from fisher_py.data import Device
    from fisher_py.data.business import Scan

    raw = RawFileReaderAdapter.file_factory(raw_path)

    if raw.is_error:
        raise RuntimeError(f"RawFileReader 打开文件失败 {raw_path}")

    try:
        raw.select_instrument(Device.MS, 1)
        header = raw.run_header_ex
        first, last = int(header.first_spectrum), int(header.last_spectrum)

        retention_times = []
        scans = []

        for scan_number in range(first, last + 1):
            event = raw.get_scan_event_for_scan_number(scan_number)

            if int(event.ms_order.value) != 1:
                continue

            scan = Scan.from_file(raw, scan_number)

            if scan.has_centroid_stream and scan.centroid_scan.length > 0:
                mz = np.asarray(scan.centroid_scan.masses, dtype=np.float64)
                intensity = np.asarray(scan.centroid_scan.intensities, dtype=np.float64)
            else:
                segmented = scan.segmented_scan
                mz = np.asarray(segmented.positions, dtype=np.float64)
                intensity = np.asarray(segmented.intensities, dtype=np.float64)

            order = np.argsort(mz, kind="mergesort")
            retention_times.append(float(raw.retention_time_from_scan_number(scan_number)))
            scans.append((mz[order], intensity[order]))

        return np.asarray(retention_times), scans

    finally:
        try:
            raw.dispose()
        except Exception:
            pass


def ms1_centroid_cache_path(raw_path, cache_dir):
    base_name = os.path.splitext(os.path.basename(raw_path))[0]
    return os.path.join(cache_dir, base_name + MS1_CENTROID_CACHE_SUFFIX)


def read_ms1_centroid_cache(cache_path):
    with np.load(cache_path) as data:
        retention_times = data["rt"]
        offsets = data["offsets"]
        mz = data["mz"]
        intensity = data["intensity"]

    scans = [
        (mz[offsets[i]:offsets[i + 1]], intensity[offsets[i]:offsets[i + 1]])
        for i in range(retention_times.size)
    ]
    return retention_times, scans


def load_ms1_scans(raw_path, cache_dir=None):
    if cache_dir:
        cache_path = ms1_centroid_cache_path(raw_path, cache_dir)
        if os.path.exists(cache_path) and os.path.getmtime(cache_path) >= os.path.getmtime(raw_path):
            try:
                return read_ms1_centroid_cache(cache_path)
            except Exception as e:
                log(f"MS1 缓存不可用，回退直接读 .raw {os.path.basename(cache_path)} {type(e).__name__}: {e}")

    return read_ms1_scans(raw_path)


def build_eic_matrix(target_mz_sorted, scans, mz_tolerance=MZ_TOLERANCE):
    n_target = len(target_mz_sorted)
    eic = np.zeros((n_target, len(scans)), dtype=np.float32)

    lower = target_mz_sorted - mz_tolerance
    upper = target_mz_sorted + mz_tolerance
    target_index = np.arange(n_target)

    for scan_i, (mz, intensity) in enumerate(scans):
        lo = np.searchsorted(mz, lower, side="left")
        hi = np.searchsorted(mz, upper, side="right")
        counts = hi - lo
        total = int(counts.sum())

        if total == 0:
            continue

        
        expanded_target = np.repeat(target_index, counts)
        starts = np.repeat(lo, counts)
        offsets = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)

        column = eic[:, scan_i]
        np.maximum.at(column, expanded_target, intensity[starts + offsets].astype(np.float32))
        eic[:, scan_i] = column

    return eic


def integrate_peak(rt_array, intensity_array, target_rt, rt_window=RT_WINDOW):
    window = (rt_array >= target_rt - rt_window) & (rt_array <= target_rt + rt_window)

    if window.sum() < 3:
        return 0.0

    x = rt_array[window]
    y = intensity_array[window].astype(float)
    smoothed = uniform_filter1d(y, size=3, mode="nearest")

    apex = int(np.argmax(smoothed))

    if smoothed[apex] <= 0:
        return 0.0

    floor = smoothed[apex] * APEX_FLOOR_RATIO

    left = apex
    while left > 0 and smoothed[left - 1] <= smoothed[left] and smoothed[left] > floor:
        left -= 1

    right = apex
    while right < len(smoothed) - 1 and smoothed[right + 1] <= smoothed[right] and smoothed[right] > floor:
        right += 1

    if right <= left:
        return 0.0

    return trapezoid(y[left:right + 1], x[left:right + 1])


def quantify_one_sample(task_args):
    raw_path, target_mz_sorted, target_rt_sorted, cache_dir = task_args
    sample_name = os.path.splitext(os.path.basename(raw_path))[0]

    try:
        start = time.time()
        retention_times, scans = load_ms1_scans(raw_path, cache_dir)
        eic = build_eic_matrix(target_mz_sorted, scans)

        areas = np.zeros(len(target_mz_sorted))
        for i in range(len(target_mz_sorted)):
            areas[i] = integrate_peak(retention_times, eic[i], target_rt_sorted[i])

        log(
            f"{sample_name} 完成，MS1扫描 {len(scans)}，非零峰面积 {(areas > 0).sum()}/{len(areas)}，"
            f"耗时 {time.time() - start:.1f} s"
        )
        return sample_name, areas

    except Exception as e:
        log(f"{sample_name} 定量失败 {type(e).__name__}: {e}")
        return sample_name, np.full(len(target_mz_sorted), np.nan)


def merge_annotations_per_peak(result, mz_col, rt_col, sample_columns, candidate_col="sub_chain_name"):
    """把“严格同一个峰”的多条注释合并成一行。

    同一个峰的判据（从严）：sample_pmz、sample_rt 完全一致，且每个样本的定量峰面积
    也完全一致——用完整精度拼成分组键，杜绝任何四舍五入把 RT/mz 相近但实为不同的峰并到一起。

    每个峰保留打分最高的“代表条”，其 sub_chain_name、total_chain_name 即最可能结果；
    其余候选的注释名称收进 candidates 列，并记录 n_annotations。代表条排序与 master_table
    去重口径一致：MS2证据 > 碎片齐全(无缺失) > IoU > RT拟合分 > cosine。
    """
    result = result.copy()
    sample_columns = list(sample_columns)

    # 完整精度拼键：mz、rt 必须一模一样，且各样本定量结果也必须完全一致，才算同一个峰
    if sample_columns:
        quant_repr = result[sample_columns].astype(str).agg("|".join, axis=1)
    else:
        quant_repr = pd.Series("", index=result.index)

    result["_peak_key"] = (
        result[mz_col].map(repr) + "@" + result[rt_col].map(repr) + "@" + quant_repr
    )

    if "evidence_type" in result.columns:
        result["_r_ev"] = (result["evidence_type"] == EVIDENCE_MS2).astype(int)
    else:
        result["_r_ev"] = 0

    if "fragment_check" in result.columns:
        result["_r_frag"] = (~result["fragment_check"].astype(str).str.contains(MISSING_TAG)).astype(int)
    else:
        result["_r_frag"] = 0

    result["_r_iou"] = pd.to_numeric(result.get("intersection_over_union"), errors="coerce").fillna(-1.0)
    result["_r_rt"] = pd.to_numeric(result.get("rt_fit_score"), errors="coerce").fillna(-1.0)
    result["_r_cos"] = pd.to_numeric(result.get("cosine_similarity"), errors="coerce").fillna(-1.0)

    rank_cols = ["_r_ev", "_r_frag", "_r_iou", "_r_rt", "_r_cos"]

    # 稳定排序后，每个峰内代表条排在最前
    result = result.sort_values(
        ["_peak_key"] + rank_cols,
        ascending=[True] + [False] * len(rank_cols),
        kind="mergesort",
    )

    grouped = result.groupby("_peak_key", sort=False)

    def collect_other_candidates(series):
        # 去重保序；索引 0 是代表条自己，1: 为其余候选
        values = list(dict.fromkeys(series.astype(str).tolist()))
        return "; ".join(values[1:])

    aggregated = pd.DataFrame({
        N_ANNOTATIONS_COLUMN: grouped.size(),
        CANDIDATES_COLUMN: grouped[candidate_col].apply(collect_other_candidates),
    }).reset_index()

    # drop_duplicates(keep="first") 保留完整的代表条整行（groupby.first 会跨行取非空，不可用）
    best = result.drop_duplicates(subset="_peak_key", keep="first").copy()
    best = best.merge(aggregated, on="_peak_key", how="left")
    best = best.drop(columns=rank_cols + ["_peak_key"], errors="ignore")

    # 把两列新信息放到 sub_chain_name 后面，便于查看
    if candidate_col in best.columns:
        cols = [c for c in best.columns if c not in (N_ANNOTATIONS_COLUMN, CANDIDATES_COLUMN)]
        insert_at = cols.index(candidate_col) + 1
        cols = cols[:insert_at] + [N_ANNOTATIONS_COLUMN, CANDIDATES_COLUMN] + cols[insert_at:]
        best = best[cols]

    return best


def quantify_master_table(
    master_table_csv,
    raw_dir,
    output_csv=None,
    max_workers=None,
    cache_dir=None,
):
    if output_csv is None:
        output_csv = os.path.join(os.path.dirname(master_table_csv), "master_table_quant.csv")

    
    if cache_dir is None:
        cache_dir = os.path.dirname(os.path.abspath(master_table_csv))

    log("Step 1 | 读取 master_table 并去除指定列")
    master = pd.read_csv(master_table_csv, encoding="utf-8-sig")
    log(f"原始记录 {len(master)}，原始列 {len(master.columns)}")

    dropped = [c for c in DROP_COLUMNS if c in master.columns]
    master = master.drop(columns=dropped)
    log(f"已去除 {dropped}")

    log("Step 2 | 构建唯一 (mz, rt) 定量目标")
    master["_mz_key"] = master["sample_pmz"].round(4)
    master["_rt_key"] = master["sample_rt"].round(2)

    targets = master[["_mz_key", "_rt_key"]].drop_duplicates().reset_index(drop=True)
    
    targets = targets.sort_values("_mz_key", kind="mergesort").reset_index(drop=True)

    target_mz = targets["_mz_key"].to_numpy(dtype=float)
    target_rt = targets["_rt_key"].to_numpy(dtype=float)
    log(f"唯一定量目标 {len(targets)}（来自 {len(master)} 行注释）")

    raw_files = sorted(glob.glob(os.path.join(raw_dir, "*.raw")))

    if not raw_files:
        raise FileNotFoundError(f"未在 {raw_dir} 找到 .raw 文件")

    if max_workers is None:
        max_workers = max(1, int(cpu_count() * 0.85))
    max_workers = max(1, min(int(max_workers), len(raw_files)))

    cached = sum(1 for path in raw_files if os.path.exists(ms1_centroid_cache_path(path, cache_dir)))
    log(f"Step 3 | 对 {len(raw_files)} 个 raw 文件提取 EIC 并积分，{max_workers} 个进程，命中 MS1 缓存 {cached}/{len(raw_files)}")

    args_list = [(path, target_mz, target_rt, cache_dir) for path in raw_files]
    sample_areas = {}

    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for sample_name, areas in executor.map(quantify_one_sample, args_list):
            sample_areas[sample_name] = areas
    log(f"定量总耗时 {time.time() - start:.1f} s")

    log("Step 4 | 把峰面积拼接到注释列后面")
    for sample_name in sorted(sample_areas):
        targets[sample_name] = sample_areas[sample_name]

    result = master.merge(targets, on=["_mz_key", "_rt_key"], how="left")
    # _mz_key/_rt_key 只是 EIC 提取用的（四舍五入）定量目标，注释合并不能用它，必须用精确 mz/rt
    result = result.drop(columns=["_mz_key", "_rt_key"])

    sample_columns = sorted(sample_areas)

    log("Step 5 | 合并严格同一峰的多条注释（mz、rt 及定量结果完全一致；代表条+候选结果列）")
    n_before = len(result)
    result = merge_annotations_per_peak(
        result,
        mz_col="sample_pmz",
        rt_col="sample_rt",
        sample_columns=sample_columns,
    )
    log(f"合并前注释行 {n_before} -> 合并后峰行 {len(result)}")

    result.to_csv(output_csv, index=False, encoding="utf-8-sig")

    log(f"定量总表已保存 {output_csv}")
    log(f"记录 {len(result)}，注释列 {len(result.columns) - len(sample_columns)}，样本列 {len(sample_columns)}")

    detected = (result[sample_columns] > 0).sum(axis=1)
    log(f"在全部 {len(sample_columns)} 个样本中均检出的行 {int((detected == len(sample_columns)).sum())}")
    log(f"一个样本都未检出的行 {int((detected == 0).sum())}")

    return result


if __name__ == "__main__":
    quantify_master_table(
        master_table_csv=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos\pkl\master_table.csv",
        raw_dir=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos",
        max_workers=4,
    )
