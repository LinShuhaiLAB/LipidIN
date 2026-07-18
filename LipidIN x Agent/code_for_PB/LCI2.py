import os
import re
from pathlib import Path
import pandas as pd
import numpy as np


def merge_msms_eq2_for_rt_ccs(
    root_dir,
    file_name="MSMS_EQ2.csv",
    mz_col="sample_pmz",
    rt_col="sample_rt",
    name_col="total_chain_name",
    adduct_col="adduct_form",
    ms2mz_col="ms2mz",
    ms2int_col="ms2int",
    mz_tol=0.01,
    rt_tol=0.2,
    output_file="rt_ccs_aligned_table.csv"
):
    root_dir = Path(root_dir)

    similarity_cols = [
        "jcr_similarity",
        "cosine_similarity",
        "intersection_over_union"
    ]

    csv_files = [
        p for p in root_dir.rglob("*.csv")
        if p.name.lower().endswith(file_name.lower())
    ]

    if len(csv_files) == 0:
        raise FileNotFoundError(f"No files ending with {file_name} found under {root_dir}")

    def _split_unique_values(values):
        out = []
        seen = set()

        for value in values:
            if pd.isna(value):
                continue

            for item in str(value).split(";"):
                item = item.strip()
                if item and item not in seen:
                    seen.add(item)
                    out.append(item)

        return out

    def _clean_ms2_text(value):
        if pd.isna(value):
            return ""

        text = str(value).strip()

        if text.lower() == "nan":
            return ""

        return text

    def _add_ms2_pair_to_peak(peak, row):
        ms2mz_text = _clean_ms2_text(row.get(ms2mz_col, np.nan))
        ms2int_text = _clean_ms2_text(row.get(ms2int_col, np.nan))

        if ms2mz_text == "" and ms2int_text == "":
            return

        pair = (ms2mz_text, ms2int_text)

        if pair not in peak["ms2_pair_seen"]:
            peak["ms2_pair_seen"].add(pair)
            peak["ms2_pairs"].append(pair)

    def _align_table(df, id_prefix):
        aligned_rows = []
        peak_index = 0

        df = df.sort_values(
            by=[name_col, rt_col, mz_col],
            ascending=[True, True, True]
        ).reset_index(drop=True)

        for lipid_name, group in df.groupby(name_col, sort=False):
            group = group.reset_index(drop=True)

            peaks = []
            rt_bins = {}

            for _, row in group.iterrows():
                mz = row[mz_col]
                rt = row[rt_col]
                rt_bin = int(np.floor(rt / rt_tol))

                candidate_peak_ids = set()
                for b in range(rt_bin - 1, rt_bin + 2):
                    candidate_peak_ids.update(rt_bins.get(b, set()))

                candidates = []

                for pid in candidate_peak_ids:
                    peak = peaks[pid]

                    mz_err = abs(mz - peak["aligned_mz"])
                    rt_err = abs(rt - peak["aligned_rt"])

                    if mz_err <= mz_tol and rt_err <= rt_tol:
                        candidates.append((pid, rt_err, mz_err))

                if len(candidates) > 0:
                    candidates.sort(key=lambda x: (x[1], x[2]))
                    best_pid = candidates[0][0]
                    peak = peaks[best_pid]

                    old_bin = peak["rt_bin"]

                    row_n = int(row.get("n_records", 1))
                    old_n = peak["n_records"]
                    new_n = old_n + row_n

                    peak["mz_sum"] += mz * row_n
                    peak["rt_sum"] += rt * row_n
                    peak["n_records"] = new_n
                    peak["aligned_mz"] = peak["mz_sum"] / new_n
                    peak["aligned_rt"] = peak["rt_sum"] / new_n

                    peak["mz_min"] = min(peak["mz_min"], row.get("mz_min", mz))
                    peak["mz_max"] = max(peak["mz_max"], row.get("mz_max", mz))
                    peak["rt_min"] = min(peak["rt_min"], row.get("rt_min", rt))
                    peak["rt_max"] = max(peak["rt_max"], row.get("rt_max", rt))

                    for adduct in _split_unique_values([row.get(adduct_col, np.nan)]):
                        peak["adduct_forms"].add(adduct)

                    _add_ms2_pair_to_peak(peak, row)

                    for sim_col in similarity_cols:
                        value = row.get(sim_col, np.nan)
                        if pd.notna(value):
                            if pd.isna(peak[sim_col]):
                                peak[sim_col] = value
                            else:
                                peak[sim_col] = max(peak[sim_col], value)

                    new_bin = int(np.floor(peak["aligned_rt"] / rt_tol))
                    if new_bin != old_bin:
                        rt_bins[old_bin].discard(best_pid)
                        rt_bins.setdefault(new_bin, set()).add(best_pid)
                        peak["rt_bin"] = new_bin

                else:
                    pid = len(peaks)
                    row_n = int(row.get("n_records", 1))

                    peak = {
                        name_col: lipid_name,
                        "aligned_mz": mz,
                        "aligned_rt": rt,
                        "mz_sum": mz * row_n,
                        "rt_sum": rt * row_n,
                        "n_records": row_n,
                        "mz_min": row.get("mz_min", mz),
                        "mz_max": row.get("mz_max", mz),
                        "rt_min": row.get("rt_min", rt),
                        "rt_max": row.get("rt_max", rt),
                        "adduct_forms": set(_split_unique_values([row.get(adduct_col, np.nan)])),
                        "ms2_pairs": [],
                        "ms2_pair_seen": set(),
                        "rt_bin": rt_bin
                    }

                    _add_ms2_pair_to_peak(peak, row)

                    for sim_col in similarity_cols:
                        peak[sim_col] = row.get(sim_col, np.nan)

                    peaks.append(peak)
                    rt_bins.setdefault(rt_bin, set()).add(pid)

            for peak in peaks:
                peak_index += 1

                output_row = {
                    "peak_id": f"{id_prefix}_{peak_index}",
                    name_col: peak[name_col],
                    "aligned_mz": peak["aligned_mz"],
                    "aligned_rt": peak["aligned_rt"],
                    "n_records": peak["n_records"],
                    "mz_min": peak["mz_min"],
                    "mz_max": peak["mz_max"],
                    "rt_min": peak["rt_min"],
                    "rt_max": peak["rt_max"],
                    adduct_col: ";".join(sorted(peak["adduct_forms"])),
                    ms2mz_col: "||".join([pair[0] for pair in peak["ms2_pairs"]]),
                    ms2int_col: "||".join([pair[1] for pair in peak["ms2_pairs"]])
                }

                for sim_col in similarity_cols:
                    output_row[sim_col] = peak[sim_col]

                aligned_rows.append(output_row)

        return pd.DataFrame(aligned_rows)

    per_file_aligned = []

    for csv_file in csv_files:
        df = pd.read_csv(csv_file)

        required_cols = [mz_col, rt_col, name_col, adduct_col]
        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            raise ValueError(f"{csv_file} missing required columns: {missing_cols}")

        if ms2mz_col not in df.columns:
            df[ms2mz_col] = ""

        if ms2int_col not in df.columns:
            df[ms2int_col] = ""

        for col in [mz_col, rt_col] + similarity_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        for sim_col in similarity_cols:
            if sim_col not in df.columns:
                df[sim_col] = np.nan

        df = df.dropna(subset=[mz_col, rt_col, name_col]).copy()

        df["n_records"] = 1
        df["mz_min"] = df[mz_col]
        df["mz_max"] = df[mz_col]
        df["rt_min"] = df[rt_col]
        df["rt_max"] = df[rt_col]

        file_result = _align_table(
            df=df,
            id_prefix=csv_file.stem
        )

        per_file_aligned.append(file_result)

    data = pd.concat(per_file_aligned, ignore_index=True)

    data = data.rename(
        columns={
            "aligned_mz": mz_col,
            "aligned_rt": rt_col
        }
    )

    result = _align_table(
        df=data,
        id_prefix="peak"
    )

    result = result.sort_values(
        by=[name_col, "aligned_rt", "aligned_mz"],
        ascending=[True, True, True]
    ).reset_index(drop=True)

    result["peak_id"] = [f"peak_{i + 1}" for i in range(len(result))]

    output_path = root_dir / output_file
    result.to_csv(output_path, index=False, encoding="utf-8-sig")

    return result

def parse_lipid_info(total_chain_name):
    text = str(total_chain_name).strip()

    if text == "" or text.lower() == "nan":
        return "", np.nan, np.nan, "total_chain_name为空，无法解析"

    parts = text.split(maxsplit=1)

    if len(parts) < 2:
        lipid_class = parts[0]
        return lipid_class, np.nan, np.nan, "未找到总脂肪酸链信息"

    lipid_class = parts[0].strip()
    total_chain = parts[1].strip()

    match = re.search(r"(\d+)\s*:\s*(\d+)", total_chain)

    if match:
        total_carbon = int(match.group(1))
        total_unsaturation = int(match.group(2))
        return lipid_class, total_carbon, total_unsaturation, "按空格和冒号解析成功"

    return lipid_class, np.nan, np.nan, "总脂肪酸链格式复杂，无法解析出碳数与不饱和度"


def add_lipid_features(df):
    df = df.copy()

    if "total_chain_name" not in df.columns:
        raise ValueError("缺少 total_chain_name 列")

    parsed = df["total_chain_name"].apply(parse_lipid_info)

    df["lipid_class"] = parsed.apply(lambda x: x[0])
    df["total_carbon"] = parsed.apply(lambda x: x[1])
    df["total_unsaturation"] = parsed.apply(lambda x: x[2])
    df["lipid_parse_note"] = parsed.apply(lambda x: x[3])

    df["total_carbon"] = pd.to_numeric(df["total_carbon"], errors="coerce")
    df["total_unsaturation"] = pd.to_numeric(df["total_unsaturation"], errors="coerce")

    return df


def get_rt_column(df):
    if "aligned_rt" in df.columns:
        return "aligned_rt"
    if "sample_rt" in df.columns:
        return "sample_rt"
    raise ValueError("缺少 aligned_rt 或 sample_rt 列")


def get_mz_column(df):
    if "aligned_mz" in df.columns:
        return "aligned_mz"
    if "sample_pmz" in df.columns:
        return "sample_pmz"
    return None


def calculate_ecn(total_carbon, total_unsaturation, db_weight=2.0):
    try:
        if pd.isna(total_carbon) or pd.isna(total_unsaturation):
            return np.nan
        return float(total_carbon) - db_weight * float(total_unsaturation)
    except Exception:
        return np.nan


def estimate_rt_from_neighbors(class_df, idx, rt_col):
    row = class_df.loc[idx]

    target_c = row["total_carbon"]
    target_u = row["total_unsaturation"]
    target_rt = row[rt_col]

    if pd.isna(target_c) or pd.isna(target_u) or pd.isna(target_rt):
        return np.nan, False, "缺少RT、总碳数或总不饱和度，无法进行ECN规则推断 | Missing RT or structural features"

    evidence_df = class_df[
        class_df[rt_col].notna() &
        class_df["total_carbon"].notna() &
        class_df["total_unsaturation"].notna()
    ].copy()

    evidence_df = evidence_df.drop(index=idx, errors="ignore")

    if evidence_df.empty:
        return np.nan, False, "同类脂质没有可比记录 | No comparable same-class records"

    if "cosine_similarity" in evidence_df.columns:
        evidence_df["_cos"] = pd.to_numeric(evidence_df["cosine_similarity"], errors="coerce").fillna(-np.inf)
    else:
        evidence_df["_cos"] = 0.0

    evidence_df = evidence_df.sort_values("_cos", ascending=False, kind="mergesort")

    violation = False
    expected_values = []
    notes = []

    same_u = evidence_df[evidence_df["total_unsaturation"] == target_u].copy()

    if not same_u.empty:
        lower_c = same_u[same_u["total_carbon"] < target_c]
        higher_c = same_u[same_u["total_carbon"] > target_c]

        if not lower_c.empty:
            lower_best = lower_c.sort_values(["total_carbon", "_cos"], ascending=[False, False], kind="mergesort").iloc[0]
            if lower_best[rt_col] >= target_rt:
                violation = True
            expected_values.append(float(lower_best[rt_col] + 0.15 * (target_c - lower_best["total_carbon"])))
            notes.append("同不饱和度下低碳数脂质应具有更低RT | Lower carbon number should have lower RT at the same unsaturation")

        if not higher_c.empty:
            higher_best = higher_c.sort_values(["total_carbon", "_cos"], ascending=[True, False], kind="mergesort").iloc[0]
            if higher_best[rt_col] <= target_rt:
                violation = True
            expected_values.append(float(higher_best[rt_col] - 0.15 * (higher_best["total_carbon"] - target_c)))
            notes.append("同不饱和度下高碳数脂质应具有更高RT | Higher carbon number should have higher RT at the same unsaturation")

    same_c = evidence_df[evidence_df["total_carbon"] == target_c].copy()

    if not same_c.empty:
        lower_u = same_c[same_c["total_unsaturation"] < target_u]
        higher_u = same_c[same_c["total_unsaturation"] > target_u]

        if not lower_u.empty:
            lower_best = lower_u.sort_values(["total_unsaturation", "_cos"], ascending=[False, False], kind="mergesort").iloc[0]
            if lower_best[rt_col] <= target_rt:
                violation = True
            expected_values.append(float(lower_best[rt_col] - 0.25 * (target_u - lower_best["total_unsaturation"])))
            notes.append("同碳数下低不饱和度脂质应具有更高RT | Lower unsaturation should have higher RT at the same carbon number")

        if not higher_u.empty:
            higher_best = higher_u.sort_values(["total_unsaturation", "_cos"], ascending=[True, False], kind="mergesort").iloc[0]
            if higher_best[rt_col] >= target_rt:
                violation = True
            expected_values.append(float(higher_best[rt_col] + 0.25 * (higher_best["total_unsaturation"] - target_u)))
            notes.append("同碳数下高不饱和度脂质应具有更低RT | Higher unsaturation should have lower RT at the same carbon number")

    if not expected_values:
        class_df_tmp = evidence_df.copy()
        class_df_tmp["ecn"] = class_df_tmp.apply(
            lambda x: calculate_ecn(x["total_carbon"], x["total_unsaturation"]),
            axis=1
        )
        target_ecn = calculate_ecn(target_c, target_u)

        class_df_tmp = class_df_tmp[class_df_tmp["ecn"].notna()].copy()

        if class_df_tmp.shape[0] >= 2 and pd.notna(target_ecn):
            class_df_tmp["ecn_distance"] = (class_df_tmp["ecn"] - target_ecn).abs()
            near = class_df_tmp.sort_values(["ecn_distance", "_cos"], ascending=[True, False], kind="mergesort").head(3)
            expected_values = near[rt_col].astype(float).tolist()
            notes.append("缺少严格同碳数或同不饱和度证据，使用邻近ECN记录作为弱证据 | Nearby ECN records were used as weak evidence")

    if not expected_values:
        return np.nan, False, "同类脂质缺少可用于RT梯度推断的记录 | Insufficient same-class RT-gradient evidence"

    expected_rt = float(np.median(expected_values))
    note = "；".join(notes)

    return expected_rt, violation, note


def calculate_rule_based_rt_score(df):
    df = df.copy()
    rt_col = get_rt_column(df)

    df["rt_expected_by_rule"] = np.nan
    df["rt_rule_error"] = np.nan
    df["rt_rule_score"] = np.nan
    df["rt_rule_note"] = ""

    if "cosine_similarity" not in df.columns:
        df["cosine_similarity"] = np.nan

    df["_cos_for_rule"] = pd.to_numeric(df["cosine_similarity"], errors="coerce").fillna(-np.inf)

    for lipid_class, class_df in df.groupby("lipid_class", sort=False):
        class_df = class_df.sort_values("_cos_for_rule", ascending=False, kind="mergesort")

        if class_df.shape[0] < 2:
            df.loc[class_df.index, "rt_rule_score"] = 0.5
            df.loc[class_df.index, "rt_rule_note"] = "同类脂质记录不足，规则判断为中性 | Insufficient same-class records"
            continue

        for idx in class_df.index:
            expected_rt, violation, note = estimate_rt_from_neighbors(class_df, idx, rt_col)
            current_rt = df.loc[idx, rt_col]

            if pd.isna(expected_rt) or pd.isna(current_rt):
                df.loc[idx, "rt_rule_score"] = 0.5
                df.loc[idx, "rt_rule_note"] = note
                continue

            rt_error = abs(float(current_rt) - float(expected_rt))

            if violation:
                score = min(0.099, np.exp(-rt_error / 0.2) * 0.099)
                rule_note = f"不满足ECN方向性规则，得分限制在0.1以下 | ECN-direction rule was violated. {note}"
            else:
                score = float(np.exp(-rt_error / 0.5))
                rule_note = f"满足ECN方向性规则，按预测RT距离给分 | ECN-direction rule was satisfied. {note}"

            df.loc[idx, "rt_expected_by_rule"] = expected_rt
            df.loc[idx, "rt_rule_error"] = rt_error
            df.loc[idx, "rt_rule_score"] = round(float(np.clip(score, 0.0, 1.0)), 4)
            df.loc[idx, "rt_rule_note"] = rule_note

    df = df.drop(columns=["_cos_for_rule"], errors="ignore")

    return df


def combine_rt_scores(df):
    df = df.copy()

    rule_score = pd.to_numeric(df["rt_rule_score"], errors="coerce").fillna(0.5)

    df["rt_reliability_score"] = rule_score.clip(0.0, 1.0).round(4)

    return df


def reorder_output_columns(df):
    rt_col = get_rt_column(df)
    mz_col = get_mz_column(df)

    first_columns = [
        "peak_id",
        mz_col,
        rt_col,
        "ms2mz",
        "ms2int",
        "total_chain_name",
        "lipid_class",
        "total_carbon",
        "total_unsaturation",
        "adduct_form",
        "n_records",
        "mz_min",
        "mz_max",
        "rt_min",
        "rt_max",
        "jcr_similarity",
        "cosine_similarity",
        "intersection_over_union",
        "rt_expected_by_rule",
        "rt_rule_error",
        "rt_rule_score",
        "rt_reliability_score",
        "lipid_parse_note",
        "rt_rule_note"
    ]

    first_columns = [col for col in first_columns if col is not None and col in df.columns]
    remaining_columns = [col for col in df.columns if col not in first_columns]

    return df[first_columns + remaining_columns]


def retention_time_inference(
    INPUT_CSV,
    source_mzml_dir=None,
    overwrite=True
):
    print("开始RT可信度评估流程")
    print(f"输入文件 {INPUT_CSV}")

    input_dir = os.path.dirname(INPUT_CSV)
    input_name = os.path.basename(INPUT_CSV)

    if not os.path.exists(INPUT_CSV):
        print("Step 0 | 未找到输入RT判断表，开始从MSMS_EQ2.csv结尾文件重新生成")

        if source_mzml_dir is None or not os.path.exists(source_mzml_dir):
            raise FileNotFoundError(
                f"未找到输入文件 {INPUT_CSV}，且 source_mzml_dir 不存在，无法自动生成RT判断表"
            )

        merged_table = merge_msms_eq2_for_rt_ccs(
            root_dir=source_mzml_dir,
            output_file=input_name
        )

        generated_path = os.path.join(source_mzml_dir, input_name)

        if os.path.abspath(generated_path) != os.path.abspath(INPUT_CSV):
            os.makedirs(input_dir, exist_ok=True)
            merged_table.to_csv(INPUT_CSV, index=False, encoding="utf-8-sig")

        print(f"Step 0 | 已生成新的RT判断输入表 {INPUT_CSV}")

    print("Step 1 | 正在读取输入表格")
    df = pd.read_csv(INPUT_CSV)
    print(f"Step 1 | 输入表格读取完成，记录数 {len(df)}")

    rt_col = get_rt_column(df)
    mz_col = get_mz_column(df)

    df[rt_col] = pd.to_numeric(df[rt_col], errors="coerce")

    if mz_col is not None:
        df[mz_col] = pd.to_numeric(df[mz_col], errors="coerce")

    if "cosine_similarity" in df.columns:
        df["cosine_similarity"] = pd.to_numeric(df["cosine_similarity"], errors="coerce")
    else:
        df["cosine_similarity"] = np.nan

    print("Step 2 | 正在解析 total_chain_name")
    df = add_lipid_features(df)

    print("Step 3 | 正在进行ECN规则RT推断")
    df = calculate_rule_based_rt_score(df)

    print("Step 4 | 正在整合RT判断得分")
    df = combine_rt_scores(df)

    df = reorder_output_columns(df)

    if overwrite:
        output_csv = INPUT_CSV
    else:
        input_base, input_ext = os.path.splitext(INPUT_CSV)
        output_csv = f"{input_base}_RT_checked{input_ext}"

    print("Step 5 | 正在保存结果表格")
    df.to_csv(output_csv, index=False, encoding="utf-8-sig")

    print(f"完成，RT判断结果已保存到 {output_csv}")

    return df
if __name__ == "__main__":
    retention_time_inference(
        INPUT_CSV=r"D:\bio_inf\LipidIN2.0_SZ\test_data\Lipid_15-30-45\pos\pkl\rt_ccs_aligned_table.csv",
        overwrite=False
    )