import os
import re
import glob
import pickle
import numpy as np
import pandas as pd


RT_ON_LINE = "在拟合线上"
RT_OFF_LINE = "偏离拟合线"
RT_NOT_APPLICABLE = "无法判断"


EVIDENCE_MS2 = "MS2注释"
EVIDENCE_MS1_UNANNOTATED = "MS1推断(unannotated)"
EVIDENCE_MS1_ONLY = "MS1推断(MS1_only)"


CHECK_ACCEPTED = ("通过", "无需搜索(IoU=1)", "无需搜索(TG IoU>0.75)")


CHECK_MISSING = "缺失"



CHECK_ACCEPTED_OR_MISSING = CHECK_ACCEPTED + (CHECK_MISSING,)

NO_MS2_EVIDENCE = "无MS2碎片证据"

MASTER_COLUMNS = [
    "sample_pmz",
    "sample_rt",
    "total_chain_name",
    "sub_chain_name",
    "subclass",
    "fa_composition",
    "fragment_check",
    "adduct_form",
    "chemical_formula",
    "sample_pmz_da_error",
    "jcr_similarity",
    "cosine_similarity",
    "intersection_over_union",
    "sample_ms2_mz",
    "sample_ms2_int",
    "rt_fit_expected",
    "rt_fit_error",
    "rt_fit_score",
    "rt_fit_status",
    "evidence_type",
    "source_raw",
    "n_merged",
    "all_source_raw",
]


def log(message):
    print(f"[INFO] {message}")


def from_excel_text(value):
    if pd.isna(value):
        return value

    text = str(value)

    if text.startswith('="') and text.endswith('"'):
        return text[2:-1]

    return text


def parse_carbon_and_unsaturation(total_chain_name):
    text = str(total_chain_name).strip()

    if text == "" or text.lower() == "nan":
        return "", np.nan, np.nan

    parts = text.split(maxsplit=1)

    if len(parts) < 2:
        return parts[0], np.nan, np.nan

    match = re.search(r"(\d+)\s*:\s*(\d+)", parts[1])

    if not match:
        return parts[0], np.nan, np.nan

    return parts[0], float(match.group(1)), float(match.group(2))


def calculate_ecn(total_carbon, total_unsaturation, db_weight=2.0):
    return total_carbon - db_weight * total_unsaturation





def fit_rt_lines(rt_checked_csv, min_points=3, trim_sigma=2.0, tol_sigma=2.0, tol_min=0.20):
    df = pd.read_csv(rt_checked_csv, encoding="utf-8-sig")

    for col in ["total_carbon", "total_unsaturation", "aligned_rt"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["lipid_class", "total_carbon", "total_unsaturation", "aligned_rt"]).copy()
    df["ecn"] = calculate_ecn(df["total_carbon"], df["total_unsaturation"])

    lines = {}

    for subclass, group in df.groupby("lipid_class", sort=True):
        x = group["ecn"].to_numpy(dtype=float)
        y = group["aligned_rt"].to_numpy(dtype=float)

        if len(x) < min_points or np.ptp(x) == 0:
            continue

        slope, intercept = np.polyfit(x, y, 1)
        residual = y - (slope * x + intercept)
        sigma = float(np.std(residual))

        
        if sigma > 0:
            keep = np.abs(residual) <= trim_sigma * sigma
            if keep.sum() >= min_points and np.ptp(x[keep]) > 0:
                x, y = x[keep], y[keep]
                slope, intercept = np.polyfit(x, y, 1)
                residual = y - (slope * x + intercept)
                sigma = float(np.std(residual))

        lines[str(subclass)] = {
            "slope": float(slope),
            "intercept": float(intercept),
            "resid_std": sigma,
            "tol": float(max(tol_min, tol_sigma * sigma)),
            "n_anchor": int(len(x)),
        }

    return lines


def judge_rt_on_line(subclass, total_carbon, total_unsaturation, observed_rt, lines):
    line = lines.get(str(subclass))

    if line is None or pd.isna(total_carbon) or pd.isna(total_unsaturation) or pd.isna(observed_rt):
        return np.nan, np.nan, np.nan, RT_NOT_APPLICABLE

    expected = line["slope"] * calculate_ecn(total_carbon, total_unsaturation) + line["intercept"]
    error = abs(float(observed_rt) - expected)
    score = float(np.exp(-error / line["tol"]))
    status = RT_ON_LINE if error <= line["tol"] else RT_OFF_LINE

    return round(expected, 4), round(error, 4), round(score, 4), status


def add_rt_judgement(df, lines):
    expected_list, error_list, score_list, status_list = [], [], [], []

    for subclass, carbon, unsaturation, rt in zip(
        df["subclass"], df["total_carbon"], df["total_unsaturation"], df["sample_rt"]
    ):
        expected, error, score, status = judge_rt_on_line(subclass, carbon, unsaturation, rt, lines)
        expected_list.append(expected)
        error_list.append(error)
        score_list.append(score)
        status_list.append(status)

    df["rt_fit_expected"] = expected_list
    df["rt_fit_error"] = error_list
    df["rt_fit_score"] = score_list
    df["rt_fit_status"] = status_list

    return df





def merge_fragment_checks(headgroup_check, fa_chain_check):
    hg = str(headgroup_check)
    fa = str(fa_chain_check)

    if hg == fa:
        return hg

    return f"头基:{hg}+脂肪链:{fa}"


def source_raw_from_path(path):
    name = os.path.basename(path)

    for suffix in ["_MS1_MSMS_EQ2.csv", "_MS1_MSMS_unannotated.csv", "_MS1_only.pkl"]:
        if name.endswith(suffix):
            return name[: -len(suffix)] + ".raw"

    return name


def load_eq2_rows(pkl_dir, lines):
    rows = []

    for csv_file in sorted(glob.glob(os.path.join(pkl_dir, "*_EQ2.csv"))):
        df = pd.read_csv(csv_file, encoding="utf-8-sig")
        df["fa_composition"] = df["fa_composition"].map(from_excel_text)

        
        
        df = df[
            df["headgroup_check"].isin(CHECK_ACCEPTED_OR_MISSING)
            & df["fa_chain_check"].isin(CHECK_ACCEPTED_OR_MISSING)
        ].copy()

        if df.empty:
            continue

        df["fragment_check"] = [
            merge_fragment_checks(hg, fa) for hg, fa in zip(df["headgroup_check"], df["fa_chain_check"])
        ]

        parsed = df["total_chain_name"].map(parse_carbon_and_unsaturation)
        df["total_carbon"] = [p[1] for p in parsed]
        df["total_unsaturation"] = [p[2] for p in parsed]

        df["evidence_type"] = EVIDENCE_MS2
        df["source_raw"] = source_raw_from_path(csv_file)

        rows.append(df)

    if not rows:
        return pd.DataFrame(columns=MASTER_COLUMNS)

    merged = pd.concat(rows, ignore_index=True)
    merged = add_rt_judgement(merged, lines)

    log(
        f"_EQ2.csv records passing the dual fragment check and RT judgement {len(merged)}"
        f"\n_EQ2.csv 通过碎片双检并完成RT判断的记录 {len(merged)}"
    )

    return merged





# MS1 库索引缓存版本号：只要 _build_library_entries 的解析逻辑变了就要 +1，让旧缓存失效
MS1_INDEX_CACHE_VERSION = 1


def _ms1_index_cache_path(library_pkl):
    return library_pkl + ".ms1index.pkl"


def _build_library_entries(library_pkl):
    """从库 pkl 解析出“与具体数据集无关”的唯一库条目。

    这是 MS1 搜库里最重的一步（载入约 2GB 谱图库 + 解析上千万行标签 + 去重），
    结果只由库文件决定，可以缓存复用。不含 pred_rt/rt_tol（那两列依赖数据集的拟合线）。
    """
    with open(library_pkl, "rb") as f:
        library_df = pickle.load(f)

    label = library_df["lab_mz"].astype(str).str.split("_DTL_", n=1).str[1]

    entries = pd.DataFrame({"pmz": library_df["pmz"].to_numpy(dtype=float), "label": label.to_numpy()})
    entries = entries.drop_duplicates().reset_index(drop=True)

    parts = entries["label"].str.split(r"_CAH_|_ADD_|_FML_", n=3, expand=True, regex=True)
    entries["total_chain_name"] = parts[0]
    entries["sub_chain_name"] = parts[1]
    entries["adduct_form"] = parts[2]
    entries["chemical_formula"] = parts[3]

    # 几十万条目里总链名只有几千种，逐行跑正则纯属重复，按唯一值解析一次再展开
    chain_names, chain_codes = np.unique(entries["total_chain_name"].to_numpy(), return_inverse=True)
    parsed = [parse_carbon_and_unsaturation(name) for name in chain_names]
    entries["subclass"] = np.array([p[0] for p in parsed], dtype=object)[chain_codes]
    entries["total_carbon"] = np.array([p[1] for p in parsed], dtype=float)[chain_codes]
    entries["total_unsaturation"] = np.array([p[2] for p in parsed], dtype=float)[chain_codes]

    # 碳数/不饱和度解析不出来的条目 MS1 无从判别，直接丢（与数据集无关，可随缓存一起固化）
    entries = entries.dropna(subset=["total_carbon", "total_unsaturation"]).reset_index(drop=True)

    return entries


def _load_or_build_library_entries(library_pkl):
    """带磁盘缓存地拿到 _build_library_entries 的结果：缓存比库文件新就直接读，否则重建并写回。"""
    cache_path = _ms1_index_cache_path(library_pkl)
    lib_mtime = os.path.getmtime(library_pkl)

    if os.path.exists(cache_path) and os.path.getmtime(cache_path) >= lib_mtime:
        try:
            with open(cache_path, "rb") as f:
                cached = pickle.load(f)
            if cached.get("version") == MS1_INDEX_CACHE_VERSION and cached.get("lib_mtime") == lib_mtime:
                log(
                    f"Hit MS1 library index cache {os.path.basename(cache_path)}, skipping full-library parse"
                    f"\n命中 MS1 库索引缓存 {os.path.basename(cache_path)}，跳过整库解析"
                )
                return cached["entries"]
            log("MS1 library index cache version or source mismatch, rebuilding\nMS1 库索引缓存版本或来源不匹配，重建")
        except Exception as e:
            log(f"MS1 library index cache unavailable, rebuilding {type(e).__name__}: {e}\nMS1 库索引缓存不可用，重建 {type(e).__name__}: {e}")

    log(
        "Building MS1 library index (slow the first time: must load and parse the whole spectral library)"
        "\n构建 MS1 库索引（首次较慢：需载入并解析整个谱图库）"
    )
    entries = _build_library_entries(library_pkl)

    try:
        with open(cache_path, "wb") as f:
            pickle.dump(
                {"version": MS1_INDEX_CACHE_VERSION, "lib_mtime": lib_mtime, "entries": entries},
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        log(f"MS1 library index cached {cache_path} (reused directly next time)\nMS1 库索引已缓存 {cache_path}（下次直接复用）")
    except Exception as e:
        log(
            f"Failed to write MS1 library index cache (does not affect this run) {type(e).__name__}: {e}"
            f"\nMS1 库索引缓存写入失败（不影响本次结果）{type(e).__name__}: {e}"
        )

    return entries


def build_ms1_library_index(library_pkl, lines):
    # 最重的整库解析走磁盘缓存；下面按数据集拟合线做的筛选与预测 RT 很快，不缓存
    entries = _load_or_build_library_entries(library_pkl)

    # 只保留有拟合线的 subclass：没有拟合线就无法用RT判别，MS1 简并度太高不能凭空收录
    entries = entries[entries["subclass"].isin(lines.keys())].copy()

    ecn = calculate_ecn(entries["total_carbon"].to_numpy(), entries["total_unsaturation"].to_numpy())
    slope = entries["subclass"].map(lambda s: lines[s]["slope"]).to_numpy(dtype=float)
    intercept = entries["subclass"].map(lambda s: lines[s]["intercept"]).to_numpy(dtype=float)

    entries["pred_rt"] = slope * ecn + intercept
    entries["rt_tol"] = entries["subclass"].map(lambda s: lines[s]["tol"]).to_numpy(dtype=float)

    entries = entries.sort_values("pmz", kind="mergesort").reset_index(drop=True)

    log(
        f"Unique library entries {len(entries)} (restricted to the {entries['subclass'].nunique()} subclasses with a fitted line)"
        f"\n库唯一条目 {len(entries)}（已限定在有拟合线的 {entries['subclass'].nunique()} 个 subclass 内）"
    )

    return entries


def search_ms1_peaks(peak_mz, peak_rt, entries, ppm_err1, da_err1):
    entry_pmz = entries["pmz"].to_numpy(dtype=float)
    entry_pred_rt = entries["pred_rt"].to_numpy(dtype=float)
    entry_tol = entries["rt_tol"].to_numpy(dtype=float)

    peak_mz = np.asarray(peak_mz, dtype=float)
    peak_rt = np.asarray(peak_rt, dtype=float)
    empty = (np.array([], dtype=int), np.array([], dtype=int))

    if peak_mz.size == 0:
        return empty

    
    tol = np.maximum(peak_mz * ppm_err1 * 1e-6, da_err1)
    lo = np.searchsorted(entry_pmz, peak_mz - tol, side="left")
    hi = np.searchsorted(entry_pmz, peak_mz + tol, side="right")

    counts = hi - lo
    np.maximum(counts, 0, out=counts)
    total = int(counts.sum())

    if total == 0:
        return empty

    peak_index = np.repeat(np.arange(peak_mz.size), counts)
    within = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
    entry_index = np.repeat(lo, counts) + within

    on_line = np.abs(peak_rt[peak_index] - entry_pred_rt[entry_index]) <= entry_tol[entry_index]

    return peak_index[on_line], entry_index[on_line]


def species_group_codes(entries):
    chain_codes = pd.factorize(entries["total_chain_name"], sort=False)[0].astype(np.int64)
    adduct_codes = pd.factorize(entries["adduct_form"], sort=False)[0].astype(np.int64)
    return chain_codes, adduct_codes


def select_best_rt_fit_per_species(peak_idx, entry_idx, rt_fit_error, chain_codes, adduct_codes):
    if peak_idx.size == 0:
        return np.array([], dtype=np.int64)

    chain = chain_codes[entry_idx]
    adduct = adduct_codes[entry_idx]
    n_chain = int(chain_codes.max()) + 1
    n_adduct = int(adduct_codes.max()) + 1

    group_key = (peak_idx.astype(np.int64) * n_chain + chain) * n_adduct + adduct

    order = np.argsort(rt_fit_error, kind="stable")
    _, first_in_group = np.unique(group_key[order], return_index=True)

    keep = order[first_in_group]
    keep.sort()
    return keep


def collapse_to_species_level(matched_df):
    if matched_df.empty:
        return matched_df

    matched_df = matched_df.sort_values("rt_fit_error", kind="mergesort")

    collapsed = matched_df.groupby(
        ["source_raw", "peak_key", "total_chain_name", "adduct_form"], as_index=False, sort=False
    ).first()

    
    collapsed["sub_chain_name"] = collapsed["total_chain_name"]
    collapsed["fa_composition"] = ""

    return collapsed


def bin_ms1_only_features(ms1_records, entries, mz_tol=0.01, rt_tol=0.2):
    entry_pmz = entries["pmz"].to_numpy(dtype=float)

    
    peak_counts = np.fromiter(
        (len(record["monoisotopic_peaks"]) for record in ms1_records),
        dtype=np.int64,
        count=len(ms1_records),
    )
    total_peaks = int(peak_counts.sum())

    if total_peaks == 0:
        return np.array([]), np.array([]), np.array([])

    record_rt = np.fromiter(
        (float(record["rt"]) for record in ms1_records),
        dtype=np.float64,
        count=len(ms1_records),
    )
    rt_all = np.repeat(record_rt, peak_counts)
    mz_all = np.fromiter(
        (peak["mz"] for record in ms1_records for peak in record["monoisotopic_peaks"]),
        dtype=np.float64,
        count=total_peaks,
    )
    intensity_all = np.fromiter(
        (peak["intensity"] for record in ms1_records for peak in record["monoisotopic_peaks"]),
        dtype=np.float64,
        count=total_peaks,
    )

    
    tol = np.maximum(mz_all * 10e-6, mz_tol)
    left = np.searchsorted(entry_pmz, mz_all - tol, side="left")
    right = np.searchsorted(entry_pmz, mz_all + tol, side="right")
    keep = right > left

    mz_all, rt_all, intensity_all = mz_all[keep], rt_all[keep], intensity_all[keep]

    if mz_all.size == 0:
        return np.array([]), np.array([]), np.array([])

    
    mz_bin = np.round(mz_all / mz_tol).astype(np.int64)
    rt_bin = np.round(rt_all / rt_tol).astype(np.int64)

    
    
    
    order = np.argsort(-intensity_all, kind="stable")

    mz_base = mz_bin.min()
    rt_base = rt_bin.min()
    rt_span = int(rt_bin.max() - rt_base) + 1
    mz_span = int(mz_bin.max() - mz_base) + 1

    if mz_span * rt_span >= 2 ** 62:
        
        feature = pd.DataFrame(
            {"mz": mz_all, "rt": rt_all, "intensity": intensity_all, "mz_bin": mz_bin, "rt_bin": rt_bin}
        )
        feature = feature.sort_values("intensity", ascending=False, kind="mergesort")
        feature = feature.groupby(["mz_bin", "rt_bin"], as_index=False, sort=False).first()
        return feature["mz"].to_numpy(), feature["rt"].to_numpy(), feature["intensity"].to_numpy()

    bin_key = (mz_bin - mz_base) * rt_span + (rt_bin - rt_base)

    
    _, first_in_bin = np.unique(bin_key[order], return_index=True)
    first_in_bin.sort()
    selected = order[first_in_bin]

    return mz_all[selected], rt_all[selected], intensity_all[selected]


def load_unannotated_rows(pkl_dir, entries, lines, ppm_err1, da_err1):
    collected = []
    chain_codes, adduct_codes = species_group_codes(entries)
    entry_pred_rt = entries["pred_rt"].to_numpy(dtype=float)

    for csv_file in sorted(glob.glob(os.path.join(pkl_dir, "*_unannotated.csv"))):
        df = pd.read_csv(csv_file)
        if df.empty:
            continue

        peak_mz = pd.to_numeric(df["sample_pmz"], errors="coerce").to_numpy(dtype=float)
        peak_rt = pd.to_numeric(df["sample_rt"], errors="coerce").to_numpy(dtype=float)

        peak_idx, entry_idx = search_ms1_peaks(peak_mz, peak_rt, entries, ppm_err1, da_err1)

        if peak_idx.size == 0:
            continue

        rt_fit_error = np.abs(peak_rt[peak_idx] - entry_pred_rt[entry_idx])
        keep = select_best_rt_fit_per_species(peak_idx, entry_idx, rt_fit_error, chain_codes, adduct_codes)
        peak_idx, entry_idx = peak_idx[keep], entry_idx[keep]

        matched = entries.iloc[entry_idx].reset_index(drop=True)
        matched["sample_pmz"] = peak_mz[peak_idx]
        matched["sample_rt"] = peak_rt[peak_idx]
        matched["sample_ms2_mz"] = df["sample_ms2_mz"].to_numpy()[peak_idx]
        matched["sample_ms2_int"] = df["sample_ms2_int"].to_numpy()[peak_idx]
        matched["peak_key"] = peak_idx
        matched["source_raw"] = source_raw_from_path(csv_file)
        matched["evidence_type"] = EVIDENCE_MS1_UNANNOTATED

        collected.append(matched)

    if not collected:
        return pd.DataFrame()

    result = pd.concat(collected, ignore_index=True)
    result["rt_fit_error"] = np.abs(result["sample_rt"] - result["pred_rt"])
    result = collapse_to_species_level(result)

    log(
        f"_unannotated.csv records kept after MS1 search and RT-fit filtering {len(result)}"
        f"\n_unannotated.csv 经 MS1 搜库并通过RT拟合筛选的记录 {len(result)}"
    )

    return result


def load_ms1_only_rows(pkl_dir, entries, lines, ppm_err1, da_err1):
    collected = []
    chain_codes, adduct_codes = species_group_codes(entries)
    entry_pred_rt = entries["pred_rt"].to_numpy(dtype=float)

    for pkl_file in sorted(glob.glob(os.path.join(pkl_dir, "*_MS1_only.pkl"))):
        with open(pkl_file, "rb") as f:
            ms1_records = pickle.load(f)

        peak_mz, peak_rt, _ = bin_ms1_only_features(ms1_records, entries)

        if peak_mz.size == 0:
            continue

        peak_idx, entry_idx = search_ms1_peaks(peak_mz, peak_rt, entries, ppm_err1, da_err1)

        if peak_idx.size == 0:
            continue

        rt_fit_error = np.abs(peak_rt[peak_idx] - entry_pred_rt[entry_idx])
        keep = select_best_rt_fit_per_species(peak_idx, entry_idx, rt_fit_error, chain_codes, adduct_codes)
        peak_idx, entry_idx = peak_idx[keep], entry_idx[keep]

        matched = entries.iloc[entry_idx].reset_index(drop=True)
        matched["sample_pmz"] = peak_mz[peak_idx]
        matched["sample_rt"] = peak_rt[peak_idx]
        matched["sample_ms2_mz"] = ""
        matched["sample_ms2_int"] = ""
        matched["peak_key"] = peak_idx
        matched["source_raw"] = source_raw_from_path(pkl_file)
        matched["evidence_type"] = EVIDENCE_MS1_ONLY

        collected.append(matched)

    if not collected:
        return pd.DataFrame()

    result = pd.concat(collected, ignore_index=True)
    result["rt_fit_error"] = np.abs(result["sample_rt"] - result["pred_rt"])
    result = collapse_to_species_level(result)

    log(
        f"MS1_only.pkl records kept after MS1 search and RT-fit filtering {len(result)}"
        f"\nMS1_only.pkl 经 MS1 搜库并通过RT拟合筛选的记录 {len(result)}"
    )

    return result


def finalize_ms1_rows(df):
    if df.empty:
        return pd.DataFrame(columns=MASTER_COLUMNS)

    df = df.copy()
    df["fragment_check"] = NO_MS2_EVIDENCE
    df["sample_pmz_da_error"] = df["sample_pmz"] - df["pmz"]
    df["jcr_similarity"] = np.nan
    df["cosine_similarity"] = np.nan
    df["intersection_over_union"] = np.nan

    df["rt_fit_expected"] = df["pred_rt"].round(4)
    df["rt_fit_error"] = df["rt_fit_error"].round(4)
    df["rt_fit_score"] = np.exp(-df["rt_fit_error"] / df["rt_tol"]).round(4)
    df["rt_fit_status"] = RT_ON_LINE

    return df





def deduplicate_master_table(df, mz_tol=0.01, rt_tol=0.2):
    df = df.copy()

    df["_iou"] = pd.to_numeric(df["intersection_over_union"], errors="coerce").fillna(-1.0)
    df["_rt"] = pd.to_numeric(df["rt_fit_score"], errors="coerce").fillna(-1.0)
    df["_cos"] = pd.to_numeric(df["cosine_similarity"], errors="coerce").fillna(-1.0)
    
    df["_evidence_rank"] = (df["evidence_type"] == EVIDENCE_MS2).astype(int)
    
    df["_frag_complete"] = (~df["fragment_check"].astype(str).str.contains(CHECK_MISSING)).astype(int)

    df = df.sort_values(
        ["sub_chain_name", "_evidence_rank", "_frag_complete", "_iou", "_rt", "_cos"],
        ascending=[True, False, False, False, False, False],
        kind="mergesort",
    )

    grouped = df.groupby("sub_chain_name", sort=False)

    best = grouped.first().reset_index()
    best["n_merged"] = grouped.size().to_numpy()
    best["all_source_raw"] = grouped["source_raw"].apply(lambda s: ";".join(sorted(set(s)))).to_numpy()

    best = best.drop(columns=["_iou", "_rt", "_cos", "_evidence_rank", "_frag_complete"], errors="ignore")

    return best


def drop_ms1_rows_confirmed_by_ms2(df):
    ms2_species = set(df.loc[df["evidence_type"] == EVIDENCE_MS2, "total_chain_name"].astype(str))

    is_ms1 = df["evidence_type"] != EVIDENCE_MS2
    redundant = is_ms1 & df["total_chain_name"].astype(str).isin(ms2_species)

    log(
        f"MS1-inferred rows dropped because MS2 evidence already exists {int(redundant.sum())}"
        f"\n因已有 MS2 证据而丢弃的 MS1 推断行 {int(redundant.sum())}"
    )

    return df[~redundant].copy()


def build_master_table(
    pkl_dir,
    library_pkl,
    rt_checked_csv=None,
    ppm_err1=10,
    da_err1=0.01,
    output_file="master_table.csv",
    include_ms1_only=True,
):
    if rt_checked_csv is None:
        rt_checked_csv = os.path.join(pkl_dir, "rt_ccs_aligned_table_RT_checked.csv")

    log("Step 1 | Fit an RT ~ ECN line per subclass\nStep 1 | 按 subclass 拟合 RT ~ ECN 直线")
    lines = fit_rt_lines(rt_checked_csv)

    summary = pd.DataFrame(lines).T.reset_index().rename(columns={"index": "subclass"})
    print(summary.round(4).to_string(index=False))

    log("Step 2 | Read _EQ2.csv and judge directly\nStep 2 | 读取 _EQ2.csv 并直接判断")
    eq2_rows = load_eq2_rows(pkl_dir, lines)

    log("Step 3 | Build the MS1-level search index\nStep 3 | 构建 MS1 级搜库索引")
    entries = build_ms1_library_index(library_pkl, lines)

    log("Step 4 | _unannotated.csv MS1 search + RT-fit judgement\nStep 4 | _unannotated.csv 做 MS1 搜库 + RT 拟合判断")
    unannotated_rows = finalize_ms1_rows(load_unannotated_rows(pkl_dir, entries, lines, ppm_err1, da_err1))

    if include_ms1_only:
        log("Step 5 | MS1_only.pkl MS1 search + RT-fit judgement\nStep 5 | MS1_only.pkl 做 MS1 搜库 + RT 拟合判断")
        ms1_only_rows = finalize_ms1_rows(load_ms1_only_rows(pkl_dir, entries, lines, ppm_err1, da_err1))
    else:
        ms1_only_rows = pd.DataFrame(columns=MASTER_COLUMNS)

    log("Step 6 | Merge all sources\nStep 6 | 合并所有来源")
    frames = [x for x in [eq2_rows, unannotated_rows, ms1_only_rows] if not x.empty]

    for frame in frames:
        for col in MASTER_COLUMNS:
            if col not in frame.columns:
                frame[col] = np.nan

    combined = pd.concat([frame[[c for c in MASTER_COLUMNS if c != "n_merged" and c != "all_source_raw"]]
                          for frame in frames], ignore_index=True)

    log(f"Combined raw records {len(combined)}\n合并后原始记录 {len(combined)}")

    combined = drop_ms1_rows_confirmed_by_ms2(combined)

    log("Step 7 | Deduplicate by sub_chain_name\nStep 7 | 按 sub_chain_name 去重")
    master = deduplicate_master_table(combined)

    master = master.sort_values(["subclass", "sample_rt", "sample_pmz"], kind="mergesort").reset_index(drop=True)
    master = master[[c for c in MASTER_COLUMNS if c in master.columns]]

    output_path = os.path.join(pkl_dir, output_file)
    master.to_csv(output_path, index=False, encoding="utf-8-sig")

    log(f"Master table saved {output_path}\n总表已保存 {output_path}")
    log(f"Master table records {len(master)}\n总表记录 {len(master)}")
    print("\nEvidence source distribution:\n证据来源分布:")
    print(master["evidence_type"].value_counts().to_string())
    print("\nRT judgement distribution:\nRT判断分布:")
    print(master["rt_fit_status"].value_counts().to_string())

    return master


if __name__ == "__main__":
    build_master_table(
        pkl_dir=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos\pkl",
        library_pkl=r"D:\bio_inf\02 LipidIN 3.0\pos_common.pkl",
        ppm_err1=10,
        da_err1=0.01,
    )
