"""Lipid differential analysis and visualization.

脂质差异分析与可视化。

Given a corrected quantification table (annotation columns + one column per sample), a group map and
a case/control pair, this produces: a differential table (log2FC, p, BH-FDR), a volcano plot, a
grouped heatmap of top features, fatty-acid chain-change panels, a lipid-class change bubble plot,
and a random-forest feature-selection + classifier with cross-validated performance. All figures go
to output_dir; a compact structured summary is returned for the conclusion engine.

输入校正后的定量表（注释列 + 每个样本一列）、分组映射与病例/对照两组，产出：差异表（log2FC、p、
BH-FDR）、火山图、Top 特征分组热图、脂肪酸链变化图、脂质分类变化气泡图，以及随机森林特征选择+分类器
（含交叉验证性能）。所有图写入 output_dir，并返回紧凑的结构化摘要供结论引擎使用。
"""

import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from scipy.cluster.hierarchy import linkage, leaves_list

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict, LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score

from master_table import parse_carbon_and_unsaturation, calculate_ecn
from qc_analysis import _configure_cjk_font
import lipid_pathways

_configure_cjk_font()

LOG2FC_THRESHOLD = 1.0
FDR_THRESHOLD = 0.05
P_THRESHOLD = 0.05
# With small replicate numbers (e.g. 3 vs 3) BH-FDR over thousands of lipids is essentially
# unreachable, so the default significance rule follows common small-n lipidomics practice:
# raw p < 0.05 AND |log2FC| >= 1. FDR is still computed and reported. Set use_fdr=True to switch.
# 小重复数（如 3 vs 3）下，对上千脂质做 BH-FDR 几乎无法显著，故默认沿用小样本脂质组学惯例：
# 原始 p < 0.05 且 |log2FC| >= 1；FDR 仍会计算并报告。将 use_fdr 设为 True 可切换为 FDR 判据。
USE_FDR_FOR_SIGNIFICANCE = False


def log(message):
    print(f"[INFO] {message}")


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p down
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(ranked, 0, 1)
    return fdr


def _feature_label(row):
    for col in ("sub_chain_name", "total_chain_name"):
        value = row.get(col)
        if isinstance(value, str) and value.strip() and value.strip().lower() != "nan":
            return value.strip()
    return str(row.name)


def differential_analysis(df, sample_names, groups, case, control,
                          log2fc_threshold=LOG2FC_THRESHOLD, fdr_threshold=FDR_THRESHOLD,
                          p_threshold=P_THRESHOLD, use_fdr=USE_FDR_FOR_SIGNIFICANCE):
    """Compute log2FC, t-test p and BH-FDR for case vs control; flag significant features.

    计算病例 vs 对照的 log2FC、t 检验 p 值与 BH-FDR，并标记显著特征。
    """
    case_cols = [s for s in sample_names if groups.get(s) == case]
    control_cols = [s for s in sample_names if groups.get(s) == control]
    if len(case_cols) < 2 or len(control_cols) < 2:
        raise ValueError(
            f"Need >=2 samples per group; case={case_cols}, control={control_cols}"
            f" / 每组至少 2 个样本"
        )

    values = df[sample_names].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    case_mat = values[case_cols].to_numpy(dtype=float)
    control_mat = values[control_cols].to_numpy(dtype=float)

    positive = values.to_numpy()
    positive = positive[positive > 0]
    eps = float(np.percentile(positive, 1)) if positive.size else 1.0

    mean_case = case_mat.mean(axis=1)
    mean_control = control_mat.mean(axis=1)
    log2fc = np.log2((mean_case + eps) / (mean_control + eps))

    # t-test on log10 space (variance-stabilized), guard zero-variance rows.
    # 在 log10 空间做 t 检验（方差稳定），并对零方差行做保护。
    case_log = np.log10(case_mat + 1.0)
    control_log = np.log10(control_mat + 1.0)
    with np.errstate(all="ignore"):
        _, pvalues = ttest_ind(case_log, control_log, axis=1, equal_var=False)
    pvalues = np.where(np.isfinite(pvalues), pvalues, 1.0)
    fdr = benjamini_hochberg(pvalues)

    result = df.copy()
    result["mean_case"] = mean_case
    result["mean_control"] = mean_control
    result["log2fc"] = log2fc
    result["pvalue"] = pvalues
    result["fdr"] = fdr
    result["neg_log10_fdr"] = -np.log10(np.clip(fdr, 1e-300, 1.0))
    result["neg_log10_p"] = -np.log10(np.clip(pvalues, 1e-300, 1.0))
    # rank_score orders features for the volcano labels, heatmap and top tables.
    # rank_score 用于火山图标注、热图与 Top 表的排序。
    result["rank_score"] = result["neg_log10_fdr"] if use_fdr else result["neg_log10_p"]
    if use_fdr:
        result["significant"] = (fdr < fdr_threshold) & (np.abs(log2fc) >= log2fc_threshold)
    else:
        result["significant"] = (pvalues < p_threshold) & (np.abs(log2fc) >= log2fc_threshold)
    result["direction"] = np.where(log2fc > 0, "up", "down")
    result.loc[~result["significant"], "direction"] = "ns"

    parsed = result.get("total_chain_name", pd.Series([""] * len(result))).map(parse_carbon_and_unsaturation)
    result["parsed_subclass"] = [p[0] for p in parsed]
    result["total_carbon"] = [p[1] for p in parsed]
    result["total_unsaturation"] = [p[2] for p in parsed]
    if "subclass" not in result.columns:
        result["subclass"] = result["parsed_subclass"]

    result["case_group"] = case
    result["control_group"] = control
    return result


def volcano_plot(diff, output_dir, tag, log2fc_threshold=LOG2FC_THRESHOLD, fdr_threshold=FDR_THRESHOLD):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    colors = {"up": "#C44E52", "down": "#4C72B0", "ns": "#BBBBBB"}
    y_col = "rank_score" if "rank_score" in diff.columns else "neg_log10_fdr"
    for direction, part in diff.groupby("direction"):
        ax.scatter(part["log2fc"], part[y_col], s=18, alpha=0.7,
                   color=colors.get(direction, "#BBBBBB"), label=direction, edgecolor="none")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1)
    ax.axvline(log2fc_threshold, color="grey", linestyle="--", linewidth=1)
    ax.axvline(-log2fc_threshold, color="grey", linestyle="--", linewidth=1)

    top = diff[diff["significant"]].reindex(
        diff[diff["significant"]][y_col].sort_values(ascending=False).index
    ).head(12)
    for _, row in top.iterrows():
        ax.annotate(_feature_label(row), (row["log2fc"], row[y_col]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel(f"log2 fold change ({diff['case_group'].iloc[0]} / {diff['control_group'].iloc[0]})")
    ax.set_ylabel("-log10 p-value")
    ax.set_title(f"Volcano plot ({tag})")
    ax.legend(frameon=False)
    fig.tight_layout()
    path = os.path.join(output_dir, f"volcano_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def heatmap_plot(diff, sample_names, groups, case, control, output_dir, tag, top_n=40):
    sig = diff[diff["significant"]].copy()
    rank_col = "rank_score" if "rank_score" in diff.columns else "neg_log10_fdr"
    ranked = sig.reindex(sig[rank_col].sort_values(ascending=False).index).head(top_n)
    if len(ranked) < 2:
        ranked = diff.reindex(diff[rank_col].sort_values(ascending=False).index).head(top_n)
    if len(ranked) < 2:
        return None

    cols = [s for s in sample_names if groups.get(s) in (case, control)]
    cols = sorted(cols, key=lambda s: (groups.get(s) != case, s))
    matrix = ranked[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    z = np.log10(matrix + 1.0)
    z = (z - z.mean(axis=1, keepdims=True)) / (z.std(axis=1, keepdims=True) + 1e-9)

    try:
        row_order = leaves_list(linkage(z, method="average"))
    except Exception:
        row_order = np.argsort(-ranked["log2fc"].to_numpy())
    z = z[row_order]
    labels = [ _feature_label(row) for _, row in ranked.iloc[row_order].iterrows() ]

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(cols) + 3), max(5, 0.28 * len(labels) + 1.5)))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([f"{s}\n{groups.get(s)}" for s in cols], fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(f"Top differential lipids ({tag})")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z-score")
    fig.tight_layout()
    path = os.path.join(output_dir, f"heatmap_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def fatty_acid_change_plot(diff, output_dir, tag):
    data = diff.dropna(subset=["total_carbon", "total_unsaturation"]).copy()
    if data.empty:
        return None, pd.DataFrame()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    color = data["direction"].map({"up": "#C44E52", "down": "#4C72B0", "ns": "#CCCCCC"})
    axes[0].scatter(data["total_carbon"], data["log2fc"], c=color, s=22, alpha=0.7, edgecolor="none")
    axes[0].axhline(0, color="grey", linewidth=0.8)
    axes[0].set_xlabel("Total acyl carbon number")
    axes[0].set_ylabel("log2 fold change")
    axes[0].set_title(f"Chain length vs change ({tag})")

    axes[1].scatter(data["total_unsaturation"], data["log2fc"], c=color, s=22, alpha=0.7, edgecolor="none")
    axes[1].axhline(0, color="grey", linewidth=0.8)
    axes[1].set_xlabel("Total double bonds")
    axes[1].set_ylabel("log2 fold change")
    axes[1].set_title(f"Unsaturation vs change ({tag})")
    fig.tight_layout()
    path = os.path.join(output_dir, f"fatty_acid_change_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)

    fa_summary = (
        data.groupby("total_unsaturation")
        .agg(n=("log2fc", "size"), mean_log2fc=("log2fc", "mean"))
        .reset_index()
    )
    return path, fa_summary


def class_bubble_plot(diff, output_dir, tag):
    data = diff.copy()
    data["subclass"] = data["subclass"].astype(str)
    grouped = data.groupby("subclass").agg(
        n=("log2fc", "size"),
        mean_log2fc=("log2fc", "mean"),
        n_sig=("significant", "sum"),
        median_fdr=("fdr", "median"),
    ).reset_index()
    grouped = grouped[grouped["n"] >= 1].sort_values("mean_log2fc")
    if grouped.empty:
        return None, grouped

    fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(grouped) + 1.5)))
    sizes = 40 + 40 * grouped["n"].to_numpy()
    color = -np.log10(np.clip(grouped["median_fdr"].to_numpy(), 1e-300, 1.0))
    scatter = ax.scatter(grouped["mean_log2fc"], range(len(grouped)), s=sizes,
                         c=color, cmap="viridis", alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(grouped)))
    ax.set_yticklabels(grouped["subclass"])
    ax.set_xlabel("Mean log2 fold change")
    ax.set_title(f"Lipid class change bubble ({tag})  (size=#features)")
    fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02, label="-log10 median FDR")
    fig.tight_layout()
    path = os.path.join(output_dir, f"class_bubble_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path, grouped


def random_forest_analysis(diff, sample_names, groups, case, control, output_dir, tag, top_k=20):
    """RF feature selection (importance) + RF classifier with cross-validated accuracy/AUC.

    随机森林特征选择（重要性）+ 随机森林分类器（交叉验证准确率/AUC）。
    """
    cols = [s for s in sample_names if groups.get(s) in (case, control)]
    y = np.array([1 if groups.get(s) == case else 0 for s in cols])
    if len(np.unique(y)) < 2 or len(cols) < 4:
        log(f"RF skipped: too few samples/classes ({cols}) / 随机森林跳过：样本或类别太少")
        return {}

    x = np.log10(diff[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float).T + 1.0)
    labels = [_feature_label(row) for _, row in diff.iterrows()]

    n_splits = min(5, np.bincount(y).min())
    forest = RandomForestClassifier(n_estimators=400, random_state=0, n_jobs=-1)
    forest.fit(x, y)
    importances = forest.feature_importances_
    order = np.argsort(-importances)
    top_idx = order[:top_k]
    selected = [(labels[i], float(importances[i])) for i in top_idx]

    # Cross-validated performance (LOO if a fold would be too small, else stratified k-fold).
    # 交叉验证性能（折太小用留一，否则用分层 k 折）。
    try:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0) if n_splits >= 2 else LeaveOneOut()
        proba = cross_val_predict(
            RandomForestClassifier(n_estimators=400, random_state=0, n_jobs=-1),
            x, y, cv=cv, method="predict_proba"
        )[:, 1]
        pred = (proba >= 0.5).astype(int)
        cv_accuracy = float(accuracy_score(y, pred))
        cv_auc = float(roc_auc_score(y, proba)) if len(np.unique(y)) == 2 else float("nan")
    except Exception as e:
        log(f"RF cross-validation failed {type(e).__name__}: {e} / 随机森林交叉验证失败")
        proba, cv_accuracy, cv_auc = None, float("nan"), float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sel_sorted = sorted(selected, key=lambda t: t[1])
    axes[0].barh([s[0] for s in sel_sorted], [s[1] for s in sel_sorted], color="#55A868")
    axes[0].set_xlabel("RF importance")
    axes[0].set_title(f"Top {len(selected)} features ({tag})")
    axes[0].tick_params(axis="y", labelsize=7)

    if proba is not None and len(np.unique(y)) == 2:
        fpr, tpr, _ = roc_curve(y, proba)
        axes[1].plot(fpr, tpr, color="#C44E52", linewidth=2, label=f"AUC={cv_auc:.2f}")
        axes[1].plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1)
        axes[1].set_xlabel("False positive rate")
        axes[1].set_ylabel("True positive rate")
        axes[1].set_title(f"CV ROC ({tag}) ACC={cv_accuracy:.2f}")
        axes[1].legend(frameon=False)
    else:
        axes[1].axis("off")
    fig.tight_layout()
    path = os.path.join(output_dir, f"random_forest_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)

    return {
        "selected_features": selected,
        "cv_accuracy": cv_accuracy,
        "cv_auc": cv_auc,
        "n_samples": len(cols),
        "plot": path,
    }


def run_analysis(df, sample_names, groups, case, control, output_dir, tag="all", language="en"):
    """Run the full differential analysis + all visualizations + LIPID MAPS pathway enrichment.

    运行完整差异分析 + 全部可视化 + LIPID MAPS 通路富集。返回结构化结果字典。
    """
    os.makedirs(output_dir, exist_ok=True)
    log(f"Differential analysis {case} vs {control} on {len(df)} features / 差异分析 {case} vs {control}，特征 {len(df)}")

    diff = differential_analysis(df, sample_names, groups, case, control)
    diff_csv = os.path.join(output_dir, f"differential_{tag}.csv")
    save_cols = [c for c in (
        "polarity", "sub_chain_name", "total_chain_name", "subclass", "adduct_form", "evidence_type",
        "fragment_check", "cosine_similarity", "intersection_over_union", "rt_fit_score", "rt_fit_status",
        "sample_pmz", "sample_rt", "total_carbon", "total_unsaturation",
        "mean_case", "mean_control", "log2fc", "pvalue", "fdr", "significant", "direction",
    ) if c in diff.columns]
    diff[save_cols].to_csv(diff_csv, index=False, encoding="utf-8-sig")

    n_sig = int(diff["significant"].sum())
    n_up = int(((diff["direction"] == "up")).sum())
    n_down = int(((diff["direction"] == "down")).sum())
    log(f"Significant features {n_sig} (up {n_up}, down {n_down}) / 显著特征 {n_sig}（上调 {n_up}，下调 {n_down}）")

    volcano = volcano_plot(diff, output_dir, tag)
    heatmap = heatmap_plot(diff, sample_names, groups, case, control, output_dir, tag)
    fa_plot, fa_summary = fatty_acid_change_plot(diff, output_dir, tag)
    bubble, class_summary = class_bubble_plot(diff, output_dir, tag)
    rf = random_forest_analysis(diff, sample_names, groups, case, control, output_dir, tag)
    pathways = lipid_pathways.analyze_pathways(diff, output_dir, tag, language=language, log=log)

    def _top(direction, k=10):
        subset = diff[(diff["significant"]) & (diff["direction"] == direction)]
        rank_col = "rank_score" if "rank_score" in subset.columns else "neg_log10_fdr"
        subset = subset.reindex(subset[rank_col].sort_values(ascending=False).index).head(k)
        return [
            {
                "name": _feature_label(row),
                "subclass": str(row.get("subclass", "")),
                "log2fc": round(float(row["log2fc"]), 3),
                "fdr": float(row["fdr"]),
                "evidence_type": str(row.get("evidence_type", "")),
                "rt_fit_score": (None if pd.isna(row.get("rt_fit_score")) else float(row.get("rt_fit_score"))),
                "cosine_similarity": (None if pd.isna(row.get("cosine_similarity")) else float(row.get("cosine_similarity"))),
                "polarity": str(row.get("polarity", "")),
            }
            for _, row in subset.iterrows()
        ]

    return {
        "tag": tag,
        "case": case,
        "control": control,
        "n_features": int(len(diff)),
        "n_significant": n_sig,
        "n_up": n_up,
        "n_down": n_down,
        "top_up": _top("up"),
        "top_down": _top("down"),
        "class_summary": class_summary.to_dict(orient="records") if class_summary is not None else [],
        "fa_summary": fa_summary.to_dict(orient="records") if fa_summary is not None else [],
        "random_forest": rf,
        "pathways": pathways,
        "differential_csv": diff_csv,
        "differential_table": diff,
        "plots": {
            "volcano": volcano,
            "heatmap": heatmap,
            "fatty_acid_change": fa_plot,
            "class_bubble": bubble,
            "random_forest": rf.get("plot"),
            "pathway_enrichment": pathways.get("plot"),
        },
    }
