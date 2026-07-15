"""Total-peak-area correction and quality control (CV / PCA / PLS-DA) for the quantification table.

对定量结果做总峰面积校正，并进行质控（CV / PCA / PLS-DA）。

This runs after quantify_master_table. It takes master_table_quant.csv, normalizes every sample to a
common total peak area (correcting for injection-volume / response drift), then evaluates data quality:
coefficient of variation (CV) on the pooled QC samples, a PCA overview of all samples, and a PLS-DA
supervised separation of the biological groups. All figures and tables are written to output_dir.

本模块在 quantify_master_table 之后运行。它读取 master_table_quant.csv，把每个样本归一化到统一的总
峰面积（校正进样量/响应漂移），再评估数据质量：基于混合 QC 样本的变异系数（CV）、所有样本的 PCA
总览，以及各生物学分组的 PLS-DA 有监督区分。所有图表与表格都会写入 output_dir。
"""

import os
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


def _configure_cjk_font():
    """Pick a CJK-capable font so bilingual figure labels render instead of showing boxes.

    选择支持中日韩字符的字体，让图中的中英文标签正常显示（而非方块）。
    """
    for name in ("Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC",
                 "Source Han Sans SC", "Arial Unicode MS"):
        try:
            font_manager.findfont(name, fallback_to_default=False)
        except Exception:
            continue
        plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        return name
    return None


_CJK_FONT = _configure_cjk_font()

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression


# Annotation columns that quantify.py / master_table.py may emit. Anything not in here and not a
# known sample name is treated cautiously; run_qc prefers an explicit sample_names list from run.py.
# quantify.py / master_table.py 可能产生的注释列。QC 优先使用 run.py 显式传入的样本名列表。
KNOWN_ANNOTATION_COLUMNS = {
    "sample_pmz", "sample_rt", "total_chain_name", "sub_chain_name", "subclass",
    "fa_composition", "fragment_check", "adduct_form", "chemical_formula",
    "sample_pmz_da_error", "jcr_similarity", "cosine_similarity", "intersection_over_union",
    "sample_ms2_mz", "sample_ms2_int", "rt_fit_expected", "rt_fit_error", "rt_fit_score",
    "rt_fit_status", "evidence_type", "source_raw", "n_merged", "all_source_raw",
    "n_annotations", "candidates",
}

QC_GROUP_TOKENS = ("qc", "pool", "pooled")
CV_PASS_THRESHOLD = 30.0  # percent / 百分数


def log(message):
    print(f"[INFO] {message}")


def group_of(sample_name):
    """Derive a group label from a sample name by stripping trailing digits/separators.

    通过去掉样本名末尾的数字/分隔符推断分组标签。
    """
    name = str(sample_name).strip()
    match = re.match(r"^(.*?)[\s_\-]*\d+$", name)
    group = match.group(1) if match and match.group(1) else name
    group = group.strip("_-. ")
    return group if group else name


def is_qc_group(group_label):
    lowered = str(group_label).lower()
    return any(token in lowered for token in QC_GROUP_TOKENS)


def detect_sample_columns(df, sample_names=None):
    """Return the list of sample (quantification) columns, preferring an explicit sample_names list.

    返回样本（定量）列的列表，优先采用显式传入的 sample_names。
    """
    if sample_names:
        present = [name for name in sample_names if name in df.columns]
        missing = [name for name in sample_names if name not in df.columns]
        if missing:
            log(
                f"These expected sample columns are absent from the table and are ignored {missing}"
                f"\n以下预期样本列不在表中，已忽略 {missing}"
            )
        return present
    # Fallback: any numeric column that is not a known annotation column.
    # 回退：所有非已知注释列的数值列。
    candidates = []
    for col in df.columns:
        if col in KNOWN_ANNOTATION_COLUMNS:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            candidates.append(col)
    return candidates


def build_feature_matrix(df, sample_columns):
    """Build a features x samples numeric matrix (NaN -> 0), indexed by a readable feature id.

    构建 特征 x 样本 的数值矩阵（NaN 记为 0），行索引为可读的特征标识。
    """
    if "sub_chain_name" in df.columns:
        labels = df["sub_chain_name"].astype(str)
    elif "total_chain_name" in df.columns:
        labels = df["total_chain_name"].astype(str)
    else:
        labels = pd.Series([f"feature_{i}" for i in range(len(df))], index=df.index)

    # Make the feature id unique so downstream tables never collide.
    # 让特征 id 唯一，避免下游表格冲突。
    feature_ids = [f"{label}__{i}" for i, label in enumerate(labels)]

    matrix = df[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix.index = feature_ids
    return matrix


def total_area_correction(matrix):
    """Scale every sample to the median total peak area. Returns (corrected, factors, totals).

    将每个样本缩放到中位总峰面积。返回 (校正后矩阵, 缩放因子, 原始总峰面积)。
    """
    totals = matrix.sum(axis=0)
    positive = totals[totals > 0]
    target = float(positive.median()) if len(positive) else 0.0

    factors = pd.Series(1.0, index=totals.index)
    if target > 0:
        for sample in totals.index:
            factors[sample] = target / totals[sample] if totals[sample] > 0 else 0.0

    corrected = matrix.mul(factors, axis=1)
    return corrected, factors, totals


def compute_cv(corrected, sample_columns, groups, output_dir, polarity):
    """Compute per-feature CV on the pooled QC samples and plot the CV distribution.

    在混合 QC 样本上计算每个特征的 CV，并绘制 CV 分布图。
    """
    qc_samples = [s for s in sample_columns if is_qc_group(groups[s])]
    if len(qc_samples) < 2:
        log(
            f"Fewer than 2 QC samples ({len(qc_samples)}), skipping CV evaluation"
            f"\nQC 样本少于 2 个（{len(qc_samples)}），跳过 CV 评估"
        )
        return {"qc_samples": qc_samples, "n_features_evaluated": 0}

    qc_matrix = corrected[qc_samples]
    mean = qc_matrix.mean(axis=1)
    std = qc_matrix.std(axis=1, ddof=1)

    evaluated = mean > 0
    cv_percent = pd.Series(np.nan, index=corrected.index)
    cv_percent[evaluated] = 100.0 * std[evaluated] / mean[evaluated]

    cv_table = pd.DataFrame({
        "feature_id": corrected.index,
        "qc_mean_area": mean.values,
        "qc_std_area": std.values,
        "cv_percent": cv_percent.values,
    })
    cv_csv = os.path.join(output_dir, f"qc_cv_{polarity}.csv")
    cv_table.to_csv(cv_csv, index=False, encoding="utf-8-sig")

    valid_cv = cv_percent[evaluated].dropna()
    n_eval = int(len(valid_cv))
    pass_30 = float((valid_cv < CV_PASS_THRESHOLD).mean() * 100.0) if n_eval else float("nan")
    pass_20 = float((valid_cv < 20.0).mean() * 100.0) if n_eval else float("nan")
    median_cv = float(valid_cv.median()) if n_eval else float("nan")

    # Plot the CV distribution with the 30% acceptance line.
    # 绘制 CV 分布图，并标注 30% 接受线。
    fig, ax = plt.subplots(figsize=(7, 5))
    if n_eval:
        clipped = valid_cv.clip(upper=100)
        ax.hist(clipped, bins=40, color="#4C72B0", edgecolor="white", alpha=0.9)
        ax.axvline(CV_PASS_THRESHOLD, color="#C44E52", linestyle="--", linewidth=1.5,
                   label=f"CV = {CV_PASS_THRESHOLD:.0f}%")
        ax.legend(frameon=False)
    ax.set_xlabel("QC CV (%)  |  QC 变异系数 (%)")
    ax.set_ylabel("Number of features  |  特征数")
    ax.set_title(f"QC CV distribution ({polarity})  |  QC 变异系数分布 ({polarity})")
    fig.tight_layout()
    cv_png = os.path.join(output_dir, f"qc_cv_distribution_{polarity}.png")
    fig.savefig(cv_png, dpi=200)
    plt.close(fig)

    log(
        f"CV: {len(qc_samples)} QC samples, {n_eval} features evaluated, median CV {median_cv:.1f}%,"
        f" CV<30% {pass_30:.1f}%, CV<20% {pass_20:.1f}%"
        f"\nCV：{len(qc_samples)} 个 QC 样本，评估特征 {n_eval} 个，CV 中位数 {median_cv:.1f}%，"
        f"CV<30% 占 {pass_30:.1f}%，CV<20% 占 {pass_20:.1f}%"
    )
    log(f"CV table saved {cv_csv}\nCV 表已保存 {cv_csv}")
    log(f"CV distribution figure saved {cv_png}\nCV 分布图已保存 {cv_png}")

    return {
        "qc_samples": qc_samples,
        "n_features_evaluated": n_eval,
        "median_cv_percent": median_cv,
        "fraction_cv_below_30": pass_30,
        "fraction_cv_below_20": pass_20,
        "cv_csv": cv_csv,
        "cv_png": cv_png,
    }


def _scale_for_multivariate(corrected, sample_columns):
    """log10 transform + drop all-zero features + z-score, returning a samples x features array.

    log10 变换 + 去除全零特征 + z-score 标准化，返回 样本 x 特征 数组。
    """
    matrix = corrected[sample_columns]
    keep = matrix.sum(axis=1) > 0
    matrix = matrix[keep]
    if matrix.shape[0] == 0:
        return None, None
    samples_by_features = np.log10(matrix.to_numpy().T + 1.0)
    scaled = StandardScaler().fit_transform(samples_by_features)
    return scaled, keep


def _group_colors(unique_groups):
    palette = plt.get_cmap("tab10")
    return {group: palette(i % 10) for i, group in enumerate(sorted(unique_groups))}


def pca_analysis(corrected, sample_columns, groups, output_dir, polarity):
    """Run a 2-component PCA over all samples and scatter the scores by group.

    对所有样本做 2 主成分 PCA，并按分组绘制得分散点图。
    """
    if len(sample_columns) < 3:
        log(
            f"Fewer than 3 samples ({len(sample_columns)}), skipping PCA"
            f"\n样本少于 3 个（{len(sample_columns)}），跳过 PCA"
        )
        return {}

    scaled, _ = _scale_for_multivariate(corrected, sample_columns)
    if scaled is None or scaled.shape[1] < 2:
        log("Not enough non-zero features for PCA\n非零特征不足，跳过 PCA")
        return {}

    n_components = min(2, scaled.shape[0], scaled.shape[1])
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(scaled)
    variance = pca.explained_variance_ratio_ * 100.0

    scores_df = pd.DataFrame({
        "sample": sample_columns,
        "group": [groups[s] for s in sample_columns],
        "PC1": scores[:, 0],
        "PC2": scores[:, 1] if n_components > 1 else np.zeros(len(sample_columns)),
    })
    scores_csv = os.path.join(output_dir, f"qc_pca_scores_{polarity}.csv")
    scores_df.to_csv(scores_csv, index=False, encoding="utf-8-sig")

    colors = _group_colors(scores_df["group"].unique())
    fig, ax = plt.subplots(figsize=(7, 6))
    for group, part in scores_df.groupby("group"):
        marker = "D" if is_qc_group(group) else "o"
        ax.scatter(part["PC1"], part["PC2"], s=90, color=colors[group], marker=marker,
                   edgecolor="black", linewidth=0.5, label=group, alpha=0.9)
    for _, row in scores_df.iterrows():
        ax.annotate(row["sample"], (row["PC1"], row["PC2"]), fontsize=8,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="grey", linewidth=0.6, linestyle=":")
    ax.axvline(0, color="grey", linewidth=0.6, linestyle=":")
    ax.set_xlabel(f"PC1 ({variance[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance[1]:.1f}%)" if n_components > 1 else "PC2")
    ax.set_title(f"PCA of all samples ({polarity})  |  全样本 PCA ({polarity})")
    ax.legend(frameon=False, title="Group | 分组")
    fig.tight_layout()
    pca_png = os.path.join(output_dir, f"qc_pca_{polarity}.png")
    fig.savefig(pca_png, dpi=200)
    plt.close(fig)

    log(
        f"PCA: PC1 {variance[0]:.1f}%"
        + (f", PC2 {variance[1]:.1f}%" if n_components > 1 else "")
        + f"\nPCA：PC1 {variance[0]:.1f}%"
        + (f"，PC2 {variance[1]:.1f}%" if n_components > 1 else "")
    )
    log(f"PCA figure saved {pca_png}\nPCA 图已保存 {pca_png}")
    log(f"PCA scores saved {scores_csv}\nPCA 得分已保存 {scores_csv}")

    return {"pca_png": pca_png, "pca_scores_csv": scores_csv,
            "explained_variance_percent": variance.tolist()}


def plsda_analysis(corrected, sample_columns, groups, output_dir, polarity):
    """Run PLS-DA on the non-QC biological groups and report R2Y / LOO-Q2 / accuracy.

    对去除 QC 后的生物学分组做 PLS-DA，并报告 R2Y / 留一 Q2 / 准确率。
    """
    bio_samples = [s for s in sample_columns if not is_qc_group(groups[s])]
    bio_groups = [groups[s] for s in bio_samples]
    classes = sorted(set(bio_groups))

    if len(classes) < 2:
        log(
            f"Fewer than 2 non-QC groups ({len(classes)}), skipping PLS-DA"
            f"\n去除 QC 后的分组少于 2 个（{len(classes)}），跳过 PLS-DA"
        )
        return {}
    if len(bio_samples) < 4:
        log(
            f"Too few non-QC samples for PLS-DA ({len(bio_samples)}), skipping"
            f"\n去除 QC 后样本太少（{len(bio_samples)}），跳过 PLS-DA"
        )
        return {}

    scaled_all, keep = _scale_for_multivariate(corrected, sample_columns)
    if scaled_all is None or scaled_all.shape[1] < 2:
        log("Not enough non-zero features for PLS-DA\n非零特征不足，跳过 PLS-DA")
        return {}

    sample_index = {s: i for i, s in enumerate(sample_columns)}
    x = scaled_all[[sample_index[s] for s in bio_samples], :]

    class_index = {c: i for i, c in enumerate(classes)}
    y_onehot = np.zeros((len(bio_samples), len(classes)))
    for row, group in enumerate(bio_groups):
        y_onehot[row, class_index[group]] = 1.0

    n_components = min(2, x.shape[0] - 1, x.shape[1])
    n_components = max(1, n_components)

    pls = PLSRegression(n_components=n_components, scale=False)
    pls.fit(x, y_onehot)
    x_scores = pls.x_scores_

    y_fit = pls.predict(x)
    ss_res = float(np.sum((y_onehot - y_fit) ** 2))
    ss_tot = float(np.sum((y_onehot - y_onehot.mean(axis=0)) ** 2))
    r2y = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # Leave-one-out cross-validation for Q2 and classification accuracy.
    # 留一交叉验证，计算 Q2 与分类准确率。
    press = 0.0
    correct = 0
    n = len(bio_samples)
    for holdout in range(n):
        train = [i for i in range(n) if i != holdout]
        cv_components = min(n_components, len(train) - 1, x.shape[1])
        cv_components = max(1, cv_components)
        model = PLSRegression(n_components=cv_components, scale=False)
        model.fit(x[train], y_onehot[train])
        pred = model.predict(x[holdout:holdout + 1])[0]
        press += float(np.sum((y_onehot[holdout] - pred) ** 2))
        if int(np.argmax(pred)) == int(np.argmax(y_onehot[holdout])):
            correct += 1
    q2 = 1.0 - press / ss_tot if ss_tot > 0 else float("nan")
    accuracy = 100.0 * correct / n

    scores_df = pd.DataFrame({
        "sample": bio_samples,
        "group": bio_groups,
        "LV1": x_scores[:, 0],
        "LV2": x_scores[:, 1] if n_components > 1 else np.zeros(len(bio_samples)),
    })
    scores_csv = os.path.join(output_dir, f"qc_plsda_scores_{polarity}.csv")
    scores_df.to_csv(scores_csv, index=False, encoding="utf-8-sig")

    colors = _group_colors(classes)
    fig, ax = plt.subplots(figsize=(7, 6))
    for group, part in scores_df.groupby("group"):
        ax.scatter(part["LV1"], part["LV2"], s=90, color=colors[group],
                   edgecolor="black", linewidth=0.5, label=group, alpha=0.9)
    for _, row in scores_df.iterrows():
        ax.annotate(row["sample"], (row["LV1"], row["LV2"]), fontsize=8,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="grey", linewidth=0.6, linestyle=":")
    ax.axvline(0, color="grey", linewidth=0.6, linestyle=":")
    ax.set_xlabel("LV1")
    ax.set_ylabel("LV2" if n_components > 1 else "LV2 (n/a)")
    ax.set_title(
        f"PLS-DA ({polarity})  R2Y={r2y:.2f} Q2={q2:.2f} ACC={accuracy:.0f}%"
        f"  |  PLS-DA ({polarity})"
    )
    ax.legend(frameon=False, title="Group | 分组")
    fig.tight_layout()
    plsda_png = os.path.join(output_dir, f"qc_plsda_{polarity}.png")
    fig.savefig(plsda_png, dpi=200)
    plt.close(fig)

    log(
        f"PLS-DA: {len(classes)} groups, {len(bio_samples)} samples, R2Y {r2y:.3f},"
        f" LOO-Q2 {q2:.3f}, LOO accuracy {accuracy:.1f}%"
        f"\nPLS-DA：{len(classes)} 个分组，{len(bio_samples)} 个样本，R2Y {r2y:.3f}，"
        f"留一 Q2 {q2:.3f}，留一准确率 {accuracy:.1f}%"
    )
    log(f"PLS-DA figure saved {plsda_png}\nPLS-DA 图已保存 {plsda_png}")
    log(f"PLS-DA scores saved {scores_csv}\nPLS-DA 得分已保存 {scores_csv}")

    return {"plsda_png": plsda_png, "plsda_scores_csv": scores_csv,
            "r2y": r2y, "q2": q2, "accuracy_percent": accuracy, "classes": classes}


def run_qc(quant_csv, output_dir, sample_names=None, polarity="", annotation_df=None):
    """Full QC entry point: total-area correction, then CV / PCA / PLS-DA, writing to output_dir.

    QC 总入口：先做总峰面积校正，再进行 CV / PCA / PLS-DA，结果写入 output_dir。
    """
    os.makedirs(output_dir, exist_ok=True)
    polarity = polarity or "data"

    log("Step QC-1 | Load quantification table and locate sample columns"
        "\n质控步骤1 | 读取定量表并定位样本列")
    df = pd.read_csv(quant_csv, encoding="utf-8-sig")
    sample_columns = detect_sample_columns(df, sample_names)
    if len(sample_columns) < 2:
        raise ValueError(
            f"Fewer than 2 sample columns detected ({sample_columns}); cannot run QC"
            f"\n检测到的样本列少于 2 个（{sample_columns}），无法进行质控"
        )
    log(
        f"Samples ({len(sample_columns)}) {sample_columns}"
        f"\n样本（{len(sample_columns)} 个）{sample_columns}"
    )

    groups = {s: group_of(s) for s in sample_columns}
    log(
        f"Group assignment {groups}"
        f"\n分组结果 {groups}"
    )

    matrix = build_feature_matrix(df, sample_columns)
    log(
        f"Feature matrix {matrix.shape[0]} features x {matrix.shape[1]} samples"
        f"\n特征矩阵 {matrix.shape[0]} 个特征 x {matrix.shape[1]} 个样本"
    )

    log("Step QC-2 | Total peak area correction (normalize each sample to the median total)"
        "\n质控步骤2 | 总峰面积校正（把每个样本归一到中位总峰面积）")
    corrected, factors, totals = total_area_correction(matrix)

    corrected_out = df.copy()
    for sample in sample_columns:
        corrected_out[sample] = corrected[sample].to_numpy()
    corrected_csv = os.path.join(output_dir, f"quant_total_area_corrected_{polarity}.csv")
    corrected_out.to_csv(corrected_csv, index=False, encoding="utf-8-sig")

    factor_table = pd.DataFrame({
        "sample": sample_columns,
        "group": [groups[s] for s in sample_columns],
        "raw_total_area": [float(totals[s]) for s in sample_columns],
        "correction_factor": [float(factors[s]) for s in sample_columns],
    })
    factor_csv = os.path.join(output_dir, f"total_area_correction_factors_{polarity}.csv")
    factor_table.to_csv(factor_csv, index=False, encoding="utf-8-sig")
    log(f"Corrected quantification table saved {corrected_csv}\n校正后定量表已保存 {corrected_csv}")
    log(f"Correction factors saved {factor_csv}\n校正因子已保存 {factor_csv}")

    log("Step QC-3 | Coefficient of variation (CV) on the QC samples"
        "\n质控步骤3 | 基于 QC 样本的变异系数（CV）")
    cv_result = compute_cv(corrected, sample_columns, groups, output_dir, polarity)

    log("Step QC-4 | PCA overview of all samples\n质控步骤4 | 全样本 PCA 总览")
    pca_result = pca_analysis(corrected, sample_columns, groups, output_dir, polarity)

    log("Step QC-5 | PLS-DA of the biological groups\n质控步骤5 | 生物学分组的 PLS-DA")
    plsda_result = plsda_analysis(corrected, sample_columns, groups, output_dir, polarity)

    # Write a compact bilingual summary.
    # 写出简洁的中英文汇总。
    summary_path = os.path.join(output_dir, f"qc_summary_{polarity}.txt")
    lines = []
    lines.append(f"QC summary ({polarity})  |  质控汇总 ({polarity})")
    lines.append(f"Samples 样本: {len(sample_columns)} -> {sample_columns}")
    lines.append(f"Groups 分组: {groups}")
    if cv_result.get("n_features_evaluated"):
        lines.append(
            f"CV 变异系数: QC={len(cv_result['qc_samples'])}, "
            f"features 特征={cv_result['n_features_evaluated']}, "
            f"median 中位={cv_result['median_cv_percent']:.1f}%, "
            f"CV<30%={cv_result['fraction_cv_below_30']:.1f}%, "
            f"CV<20%={cv_result['fraction_cv_below_20']:.1f}%"
        )
    else:
        lines.append("CV 变异系数: skipped 已跳过 (需要至少 2 个 QC 样本 / needs >=2 QC samples)")
    if pca_result:
        var = pca_result["explained_variance_percent"]
        lines.append("PCA: " + ", ".join(f"PC{i+1}={v:.1f}%" for i, v in enumerate(var)))
    else:
        lines.append("PCA: skipped 已跳过")
    if plsda_result:
        lines.append(
            f"PLS-DA: groups 分组={plsda_result['classes']}, "
            f"R2Y={plsda_result['r2y']:.3f}, Q2={plsda_result['q2']:.3f}, "
            f"accuracy 准确率={plsda_result['accuracy_percent']:.1f}%"
        )
    else:
        lines.append("PLS-DA: skipped 已跳过")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    log(f"QC summary saved {summary_path}\n质控汇总已保存 {summary_path}")
    log(f"All QC results are under {output_dir}\n全部质控结果位于 {output_dir}")

    return {
        "output_dir": output_dir,
        "corrected_csv": corrected_csv,
        "factor_csv": factor_csv,
        "cv": cv_result,
        "pca": pca_result,
        "plsda": plsda_result,
        "summary": summary_path,
    }
