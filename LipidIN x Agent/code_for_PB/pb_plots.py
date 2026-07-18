# -*- coding: utf-8 -*-
"""
第六阶段：PB 结果可视化。

图表契约（先立后画）：
  核心主张 —— PB 提供了常规脂质组学不可见的“物质内部双键位置构成”维度；本数据中该维度
      在 A/B 间呈现描述性的位置偏好重排，但受 ~50% 的 MS2 缺失与 n=5/6 限制，未达统计
      显著。图必须同时呈现这两面，不能只画“发现”而藏起功效不足。
  证据链 ——
      Fig1 (PB 独有维度)  a: 位置偏好斜率图(hero，交叉=偏好翻转)  b: 解耦散点(证明位置轴
                          与总量轴正交)  c: 代表物质逐样本点图(个体级证据)  d: Δ份额热图(全局)
      Fig2 (质控)         a: 每样品覆盖度  b: 进样重复性(同一样品两针)  c: 覆盖度分布
  原型 —— quantitative grid + hero panel。
  诚实性入图 —— “两组都有实测”的可信翻转与“0%/100%”的缺失伪翻转必须在图上可分辨；
      物质级检验 FDR 无显著，故 b 图不画显著性阈线，只画效应量。
  导出 —— SVG/PDF 可编辑文字 + 600dpi TIFF。图内标签用英文(投稿规范, 且避免 CJK 缺字)。
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

PALETTE = {
    "blue_main": "#0F4D92", "blue_secondary": "#3775BA",
    "red_strong": "#B64342", "red_2": "#E9A6A1",
    "neutral_light": "#CFCECE", "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D", "neutral_black": "#272727",
    "teal": "#42949E", "violet": "#9A4D8E", "gold": "#FFD700",
}
GROUP_COLOR = {"A": PALETTE["blue_main"], "B": PALETTE["red_strong"]}
# 双键位置 n-x 是**有序**变量(离甲基端的距离)，故用连续色标而非分类色，让 n-3->n-21
# 的顺序在视觉上自然可读。
POS_CMAP = LinearSegmentedColormap.from_list(
    "omega", ["#0F4D92", "#42949E", "#8BCF8B", "#FFD700", "#B64342"])


def panel_letters(fig, mapping):
    """布局定稿后再按 axes 实际 bbox 放面板字母。

    直接用 ax.transAxes 放字母会与 y 轴标签相撞(标签越长撞得越狠)，而 constrained_layout
    会在 draw 时才决定 axes 最终位置，故必须先 draw、再取 bbox、最后用 fig 坐标落字。
    """
    fig.canvas.draw()
    for ax, letter in mapping.items():
        bb = ax.get_position()
        fig.text(bb.x0 - 0.052, bb.y1 + 0.018, letter, fontsize=10,
                 fontweight="bold", va="bottom", ha="left")


def save_pub(fig, path, dpi=600):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(f"{path}.svg", bbox_inches="tight")
    fig.savefig(f"{path}.pdf", bbox_inches="tight")
    fig.savefig(f"{path}.png", dpi=dpi, bbox_inches="tight")
    print(f"[PB-fig] 已导出 {path}.svg / .pdf / .png")


def _pos_num(s):
    return int(str(s).split("-")[1])


def load_results(d):
    prof = pd.read_csv(os.path.join(d, "PB_position_profile.csv"))
    summ = pd.read_csv(os.path.join(d, "PB_position_preference.csv"))
    sp = pd.read_csv(os.path.join(d, "PB_species_diff.csv"))
    qc = pd.read_csv(os.path.join(d, "PB_qc.csv"))
    rep = pd.read_csv(os.path.join(d, "PB_qc_replicate.csv"))
    quant = pd.read_csv(os.path.join(d, "PB_dbposition_quant.csv"))

    # 可信翻转：两组都对 ≥2 个位置有非零实测。0%/100% 那种翻转多半是该组 DDA 没采到，
    # 不是真的不存在——必须能在图上区分，否则读者会把缺失读成生物学差异。
    def both_detected(g):
        a, b = g["share_A"].to_numpy(), g["share_B"].to_numpy()
        return bool(((a > 0.02) & (b > 0.02)).sum() >= 2)
    rel = (prof.groupby(["species_key", "f"]).apply(both_detected, include_groups=False)
           .rename("both_detected").reset_index())
    summ = summ.merge(rel, on=["species_key", "f"], how="left")
    return prof, summ, sp, qc, rep, quant


# ============================================================
# Figure 1 | PB 独有的位置维度
# ============================================================
def figure1(prof, summ, sp, quant, out):
    # 三行面板各自带标题/轴标签/脚注，间距压太紧会让 (a) 的 A/B 刻度撞上 (b)(c) 的标题、
    # (b) 的 x 标签撞上 (d) 的标题。constrained_layout 的 hspace/wspace 是相对 axes 的比例，
    # h_pad/w_pad 是英寸——两者都要给足。
    fig = plt.figure(figsize=(7.2, 7.6), layout="constrained")
    fig.get_layout_engine().set(hspace=0.22, wspace=0.10, h_pad=0.14, w_pad=0.06)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.05, 1.25])

    # ---- (a) HERO：位置偏好斜率图。线交叉 = 优势位置在 A/B 间翻转 ----
    top = (summ[summ["preference_flip"] & summ["both_detected"]]
           .sort_values("max_abs_delta", ascending=False).head(4))
    axes_a = [fig.add_subplot(gs[0, i]) for i in range(4)]
    for ax, (_, r) in zip(axes_a, top.iterrows()):
        sub = prof[(prof.species_key == r.species_key) & (prof.f == r.f)].copy()
        sub["pn"] = sub["double_bond_position"].map(_pos_num)
        sub = sub.sort_values("pn")
        norm = Normalize(2, 22)
        for _, p in sub.iterrows():
            c = POS_CMAP(norm(p["pn"]))
            ax.plot([0, 1], [p["share_A"], p["share_B"]], "-o", color=c,
                    lw=1.6, ms=3.5, mec="white", mew=0.5, zorder=3, clip_on=False)
            # 直接标注优于图例：位置固定，图例只会增加眼动
            ax.annotate(p["double_bond_position"], (0, p["share_A"]),
                        xytext=(-4, 0), textcoords="offset points",
                        ha="right", va="center", fontsize=6, color=c)
            ax.annotate(p["double_bond_position"], (1, p["share_B"]),
                        xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=6, color=c)
        ax.set_xlim(-0.42, 1.42)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["A", "B"])
        ax.tick_params(axis="x", length=0, pad=2)
        for lbl, col in zip(ax.get_xticklabels(), [GROUP_COLOR["A"], GROUP_COLOR["B"]]):
            lbl.set_color(col); lbl.set_fontweight("bold")
        ax.set_yticks([0, 0.5, 1.0])
        ax.spines["bottom"].set_visible(False)
        ax.set_title(f"{r.total_chain_name}\nf = {int(r.f)}", fontsize=6.5, pad=3)
        if ax is axes_a[0]:
            ax.set_ylabel("Share within f layer")
        else:
            ax.set_yticklabels([])
    axes_a[1].text(1.06, 1.42,
                   "Crossing lines = the dominant double-bond position flips between groups",
                   transform=axes_a[1].transAxes, ha="center", fontsize=6.5,
                   color=PALETTE["neutral_dark"], style="italic")

    # ---- (b) 解耦散点：常规轴(总量) vs PB 轴(位置份额) ----
    ax = fig.add_subplot(gs[1, :2])
    m = summ.merge(sp[["species_key", "log2FC"]], on="species_key", how="left").dropna(
        subset=["log2FC", "max_abs_delta"])
    m = m[m["log2FC"].abs() < 5]
    base = m[~m["preference_flip"]]
    flip = m[m["preference_flip"] & m["both_detected"]]
    flip_lo = m[m["preference_flip"] & ~m["both_detected"]]
    ax.scatter(base["log2FC"], base["max_abs_delta"], s=9, c=PALETTE["neutral_light"],
               edgecolors="none", label="No flip", zorder=2)
    ax.scatter(flip_lo["log2FC"], flip_lo["max_abs_delta"], s=16, facecolors="none",
               edgecolors=PALETTE["neutral_mid"], lw=0.7,
               label="Flip, one group undetected", zorder=3)
    ax.scatter(flip["log2FC"], flip["max_abs_delta"], s=22, c=PALETTE["violet"],
               edgecolors="white", lw=0.4, label="Flip, both groups measured", zorder=4)
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.6, ls="--", zorder=1)
    # 轴标签压成单行：两行标签会顶到下一行面板的标题上（constrained_layout 也救不回来）
    ax.set_xlabel("Species-level log$_2$FC (A/B)  —  conventional axis")
    ax.set_ylabel("Max |Δ share| within f layer\n—  PB-only axis")
    ax.legend(loc="upper left", fontsize=5.8, handletextpad=0.4, borderpad=0.2)
    ax_b = ax
    # 点云集中在 x≈0、y<0.6，右上角是空的 —— 注释放这里不压数据
    ax.text(0.98, 0.97, "No species reached FDR < 0.05 on the\nconventional axis (n=5 vs 6)",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
            color=PALETTE["neutral_mid"], style="italic")

    # ---- (c) 代表物质逐样本点图 ----
    ax = fig.add_subplot(gs[1, 2:])
    _panel_example(ax, prof, summ, quant, top)
    ax_c = ax

    # ---- (d) Δshare 热图：物质 × 位置 ----
    ax = fig.add_subplot(gs[2, :])
    _panel_heatmap(ax, prof, summ, fig)
    panel_letters(fig, {axes_a[0]: "a", ax_b: "b", ax_c: "c", ax: "d"})

    save_pub(fig, out)
    plt.close(fig)


def _panel_example(ax, prof, summ, quant, top):
    """最强的可信翻转组：逐生物样本的份额，点=样本，横线=组均值。"""
    r = top.iloc[0]
    sub = prof[(prof.species_key == r.species_key) & (prof.f == r.f)].copy()
    sub["pn"] = sub["double_bond_position"].map(_pos_num)
    sub = sub.sort_values("pn")
    q = quant[(quant.species_key == r.species_key) & (quant.f == r.f)]
    acols = [c for c in quant.columns if c.endswith("__area")]
    samples = [c[:-6] for c in acols]
    import re
    design = {s: re.match(r"^PB\d+_([AB])\d+$", s).group(1)
              for s in samples if re.match(r"^PB\d+_([AB])\d+$", s)}
    rng = np.random.default_rng(0)
    for i, (_, p) in enumerate(sub.iterrows()):
        row = q[q.double_bond_position == p["double_bond_position"]]
        for grp, dx in [("A", -0.17), ("B", 0.17)]:
            vals = []
            for s, g in design.items():
                if g != grp or row.empty:
                    continue
                tot = q[[f"{s}__area"]].to_numpy(dtype=float).ravel()
                if np.all(np.isnan(tot)):
                    continue
                tot = np.nan_to_num(tot, nan=0.0).sum()
                v = row[f"{s}__area"].to_numpy(dtype=float)[0]
                if tot > 0:
                    vals.append((0.0 if np.isnan(v) else v) / tot)
            if not vals:
                continue
            x = i + dx + rng.uniform(-0.045, 0.045, len(vals))
            ax.scatter(x, vals, s=11, c=GROUP_COLOR[grp], alpha=0.75,
                       edgecolors="none", zorder=3)
            ax.hlines(np.mean(vals), i + dx - 0.11, i + dx + 0.11,
                      color=GROUP_COLOR[grp], lw=1.6, zorder=4)
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["double_bond_position"])
    ax.set_ylim(-0.05, 1.08)
    ax.set_ylabel("Share within f layer")
    ax.set_title(f"{r.total_chain_name}  (f = {int(r.f)})", fontsize=6.5, pad=3)
    # 图例放轴外会撞标题；两组的点分居两侧、位置固定，直接色标注文字即可，省掉图例
    ax.text(0.02, 0.97, "Group A", transform=ax.transAxes, ha="left", va="top",
            fontsize=6, fontweight="bold", color=GROUP_COLOR["A"])
    ax.text(0.02, 0.89, "Group B", transform=ax.transAxes, ha="left", va="top",
            fontsize=6, fontweight="bold", color=GROUP_COLOR["B"])
    ax.text(0.98, 0.97, "Each dot = one sample\n(2 injections merged);\nbar = group mean",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
            color=PALETTE["neutral_mid"], style="italic")


def _panel_heatmap(ax, prof, summ, fig):
    keep = summ[summ["both_detected"]].sort_values("max_abs_delta", ascending=False).head(26)
    sub = prof.merge(keep[["species_key", "f"]], on=["species_key", "f"])
    sub["row"] = sub["total_chain_name"] + "  f" + sub["f"].astype(int).astype(str)
    sub["pn"] = sub["double_bond_position"].map(_pos_num)
    M = sub.pivot_table(index="row", columns="pn", values="delta_share_A_minus_B",
                        aggfunc="first")
    order = (sub.groupby("row")["delta_share_A_minus_B"].apply(lambda x: np.nanmax(np.abs(x)))
             .sort_values(ascending=False).index)
    M = M.loc[order]
    vmax = float(np.nanmax(np.abs(M.to_numpy())))
    im = ax.imshow(M.to_numpy().T, aspect="auto", cmap="RdBu_r",
                   norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax))
    ax.set_xticks(range(len(M.index)))
    ax.set_xticklabels(M.index, rotation=55, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(M.columns)))
    ax.set_yticklabels([f"n-{c}" for c in M.columns], fontsize=6)
    ax.set_ylabel("Double-bond position")
    # 不给 (d) 加标题：它横跨整幅宽度，标题会与上一行 (b) 的 x 轴标签撞在同一高度上，
    # 而面板字母 + y 轴标签 + colorbar 标签已经把这张图说清楚了。
    ax.set_facecolor(PALETTE["neutral_light"])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.022, aspect=12)
    cb.set_label("Δ share (A − B) within f layer\ngrey = not detected", fontsize=6)
    cb.ax.tick_params(labelsize=5.5, width=0.6)
    cb.outline.set_linewidth(0.6)


# ============================================================
# Figure 2 | 质控
# ============================================================
def figure2(qc, rep, quant, out):
    fig = plt.figure(figsize=(7.2, 2.7), layout="constrained")
    fig.get_layout_engine().set(wspace=0.08, w_pad=0.05)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.05, 1.0])

    # ---- (a) 每样本覆盖度 ----
    # 只画 coverage，不画 missing rate：missing_rate = 1 − n_quantified/总单元数，两者是
    # 同一个量的两种写法，画两遍既冗余又逼出一个让人分不清哪个点读哪个轴的双 x 轴。
    ax = fig.add_subplot(gs[0, 0])
    q = qc.copy()
    q["g"] = q["group"].fillna("QC")
    q = q.sort_values(["g", "sample"])
    y = np.arange(len(q))
    ax.barh(y, q["n_dbpos_quantified"],
            color=[GROUP_COLOR.get(g, PALETTE["neutral_mid"]) for g in q["g"]],
            height=0.74, lw=0)
    med = q["n_dbpos_quantified"].median()
    ax.axvline(med, color=PALETTE["neutral_black"], lw=0.8, ls="--", zorder=3)
    # 中位线标注钉在数据 x、轴顶 y（blended 变换），避免掉到轴外去撞刻度标签
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(med, 1.005, f"median {med:.0f}", transform=tr, fontsize=5.8,
            ha="center", va="bottom", color=PALETTE["neutral_black"])
    ax.set_yticks(y)
    ax.set_yticklabels(q["sample"], fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xlabel("Double-bond positions quantified (of 1,197)")
    ax.set_title("Coverage per sample", fontsize=6.5, pad=10)
    # 柱子铺满整个绘图区，轴内没有放图例的空位 -> 放到 x 轴标签下方
    ax.legend(handles=[mpl.patches.Patch(fc=GROUP_COLOR[g], label=f"Group {g}") for g in "AB"]
              + [mpl.patches.Patch(fc=PALETTE["neutral_mid"], label="Not in A/B design")],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=5.5,
              handlelength=1.1, handletextpad=0.4, borderpad=0.0, columnspacing=1.0)
    ax_a = ax

    # ---- (b) 批次重复性 ----
    ax = fig.add_subplot(gs[0, 1])
    subj = rep.sort_values("pearson_log")
    y = np.arange(len(subj))
    ax.hlines(y, 0, subj["pearson_log"], color=PALETTE["neutral_light"], lw=1.2, zorder=1)
    sc = ax.scatter(subj["pearson_log"], y, s=20, c=subj["n_shared"], cmap="cividis",
                    edgecolors="white", lw=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(subj["sample"], fontsize=5.5)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Pearson $r$  (injection 1 vs 2, log$_{10}$)")
    ax.set_title("Injection reproducibility", fontsize=6.5, pad=10)
    m = np.median(subj["pearson_log"])
    ax.axvline(m, color=PALETTE["red_strong"], lw=0.8, ls="--", zorder=2)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(m, 1.005, f"median {m:.2f}", transform=tr, fontsize=5.8, ha="center",
            va="bottom", color=PALETTE["red_strong"])
    # colorbar 必须给独立的 cax，否则 constrained_layout 会把它挤进相邻面板、
    # 与 (c) 的 y 轴标签叠成一团。
    cb = fig.colorbar(sc, ax=ax, location="right", pad=0.02, fraction=0.055, aspect=14)
    cb.set_label("Units shared by\nboth injections", fontsize=5.5)
    cb.ax.tick_params(labelsize=5, width=0.6)
    cb.outline.set_linewidth(0.6)
    ax_b = ax

    # ---- (c) 覆盖度分布 ----
    ax = fig.add_subplot(gs[0, 2])
    n = quant["n_samples_quantified"].to_numpy()
    ax.hist(n, bins=np.arange(-0.5, 25.5, 1), color=PALETTE["blue_secondary"], lw=0)
    ax.set_xlabel("Samples in which a\nposition is quantified")
    ax.set_ylabel("Number of positions")
    ax.set_title("Detection is sparse by design", fontsize=6.5, pad=10)
    ax.set_xlim(-0.6, 24.6)
    med = int(np.median(n))
    ax.axvline(med, color=PALETTE["red_strong"], lw=0.8, ls="--")
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(med, 1.005, f"median {med}/24", transform=tr, fontsize=5.8, ha="center",
            va="bottom", color=PALETTE["red_strong"])
    # 注释落在 x>13 的空白区，避开 x=12 那根最高柱
    ax.text(0.98, 0.42, "DDA fragments only a\nsubset per injection —\nnot a filtering artefact",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
            color=PALETTE["neutral_mid"], style="italic")
    panel_letters(fig, {ax_a: "a", ax_b: "b", ax: "c"})

    save_pub(fig, out)
    plt.close(fig)


def run_pb_plots(pb_output_dir, fig_dir=None):
    if fig_dir is None:
        fig_dir = os.path.join(pb_output_dir, "figures")
    prof, summ, sp, qc, rep, quant = load_results(pb_output_dir)
    figure1(prof, summ, sp, quant, os.path.join(fig_dir, "Fig1_position_dimension"))
    figure2(qc, rep, quant, os.path.join(fig_dir, "Fig2_quality_control"))
    print(f"[PB-fig] 图已输出到 {fig_dir}")


PB_OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "sample", "PB", "pkl")


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass
    run_pb_plots(PB_OUT_DIR)
