# -*- coding: utf-8 -*-
"""
第五阶段：PB 定量结果的质控 + 差异分析。

差异分析刻意分成**两层**，因为这批数据的价值全在第二层：

  第一层 物质总量（PB_species_quant.csv）——**常规脂质组学也能拿到这一层**。
      比的是 "PE 38:4 在 A 组比 B 组多/少"。PB 在这里没有任何新增信息。

  第二层 双键位置组成（PB_dbposition_quant.csv）——**PB 独有、常规看不见**。
      比的不是量，是**同一个物质内部各双键位置的份额**：
          share(物质, 链, n-x, 样本) = area(物质,链,n-x,样本) / area(物质,样本)
      share 已在物质内部归一，**与该物质总量无关**——物质总量翻倍不会动 share。
      所以 share 的组间变化 = 双键位置分布重排，常规方法完全测不到。

  ★ 最终交付的关键分类（quadrant 列）：把每个 (物质,链,位置) 单元按
    "总量是否显著变" × "位置份额是否显著变" 分成四类：
      - 位置独有  ：总量不显著 + 份额显著 -> **只有 PB 能发现的**，本管线的核心产出
      - 量与位置同变：两者都显著 -> 常规能看到量变，PB 额外解释了位置怎么变的
      - 仅总量    ：只有量变 -> 常规方法即可，PB 无新增
      - 无变化
    另出 omega 家族指数（n-3/n-6/n-9 各自的份额之和、n-3/n-6 比值）——脂质生物学常用，
    也是常规方法拿不到的。

实验设计（默认从样本名解析，形如 PB1_A1 / PB2_B119）：
  batch = PB1/PB2 —— **同一批生物样本的两次测量**（PB1_A1 与 PB2_A1 是同一个 A1），
          故先按生物样本对两批取均值再做组间比较，批次效应被平均掉；
          批次自身的重复性在 QC 里单独评估。
  group = A/B —— 比较的两组。不符合该命名的样本（如 PB_MSMS1/2）不进差异分析，
          只在 QC 里列出。

  也可由调用方传入 `design`（{样本名: 组名}）覆盖上面的命名约定，并用 `case`/`control`
  指定比较哪两组——agent 让用户用自然语言说"用哪些文件、跟哪些比"就是这么落地的。
  **不在 design 里的样本 group=None**，即被排除出差异分析（QC 仍会列出）。

  ★ 输出表里的 `_A`/`_B` 后缀是**位置**而非字面组名：A 槽 = case，B 槽 = control。
    默认 case/control 就是 "A"/"B"，所以不传 design 时输出与历史完全一致，pb_plots
    也照常工作；真实组名记在 `PB_design.json` 里，由报告层负责翻译回去。

缺失值：按用户要求**不做覆盖度过滤**——"有就是有，没有就是没有"。NaN 一律参与不了
统计，故每个检验都记录实际用上的 n_A/n_B；两组任一 n<MIN_N 时 p 值给 NaN 并在
`test_note` 里说明，而不是硬算一个没有意义的 p。
"""

import os
import re
import json

import numpy as np
import pandas as pd
from scipy import stats

MIN_N = 3           # 每组至少要有这么多有值样本才做检验，否则 p 给 NaN
ALPHA = 0.05        # BH-FDR 显著性阈值
SAMPLE_RE = re.compile(r"^([A-Za-z]+)(\d+)$")   # 定量表的列已是样品名(A1/B119)，进样已合并
DEFAULT_CASE = "A"      # A 槽的默认组名
DEFAULT_CONTROL = "B"   # B 槽的默认组名


def log(msg):
    print(f"[PB-QC] {msg}")


def parse_design(area_columns, design=None):
    """从 <sample>__area 列名解析实验设计。返回 DataFrame(sample,batch,group,subject)。

    design: {样本名: 组名} 时覆盖命名约定；组名为 None 或样本不在 design 里 -> 排除出
    差异分析（QC 仍列出）。不传则退回按样本名前缀推断，行为与历史一致。
    """
    rows = []
    for col in area_columns:
        s = col[:-len("__area")]
        if design is not None:
            g = design.get(s)
        else:
            m = SAMPLE_RE.match(s)
            # 不符合 <组><编号> 的（PB_MSMS1/2 等）不进差异分析
            g = m.group(1) if m else None
        # 进样已在 pb_quantify 里合并，故一个样品 = 一个生物样本 = 一列
        rows.append({"sample": s, "batch": None,
                     "group": g, "subject": s if g else None})
    return pd.DataFrame(rows)


def bh_fdr(p):
    """Benjamini-Hochberg FDR。NaN 原样返回 NaN。"""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(pv)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


# ============================================================
# 质控
# ============================================================
def _injection_reproducibility(pb_dir):
    """同一样品两针之间的一致性。回到逐针的 *_PB_dbpos.csv（定量表里进样已合并）。"""
    from pb_quantify import pool_dbpos_tables

    pooled = pool_dbpos_tables(pb_dir)
    key = ["total_chain_name", "f", "double_bond_position"]
    rows = []
    for smp, g in pooled.groupby("sample"):
        runs = sorted(g["run"].unique())
        if len(runs) < 2:
            continue
        w = (g.groupby(key + ["run"])["pair_total_intensity"].max()
             .unstack("run"))
        a, b = w[runs[0]].to_numpy(dtype=float), w[runs[1]].to_numpy(dtype=float)
        both = ~np.isnan(a) & ~np.isnan(b) & (a > 0) & (b > 0)
        if both.sum() < 3:
            continue
        la, lb = np.log10(a[both]), np.log10(b[both])
        rows.append({
            "sample": smp, "run1": runs[0], "run2": runs[1],
            "n_run1": int(np.sum(~np.isnan(a) & (a > 0))),
            "n_run2": int(np.sum(~np.isnan(b) & (b > 0))),
            "n_shared": int(both.sum()),
            "n_union": int(np.sum((~np.isnan(a) & (a > 0)) | (~np.isnan(b) & (b > 0)))),
            "pearson_log": stats.pearsonr(la, lb)[0],
            "spearman_log": stats.spearmanr(la, lb)[0],
            "median_rel_diff": float(np.nanmedian(
                np.abs(a[both] - b[both]) / ((a[both] + b[both]) / 2))),
        })
    return pd.DataFrame(rows)


def qc_report(dbpos_csv, species_csv, output_csv=None, design=None):
    db = pd.read_csv(dbpos_csv)
    sp = pd.read_csv(species_csv)
    area_cols = [c for c in db.columns if c.endswith("__area")]
    design = parse_design(area_cols, design)

    log("=" * 56)
    log("质控")
    log("=" * 56)

    rows = []
    for _, d in design.iterrows():
        v = db[f"{d['sample']}__area"].to_numpy(dtype=float)
        sv = sp[f"{d['sample']}__area"].to_numpy(dtype=float)
        rows.append({
            "sample": d["sample"], "batch": d["batch"], "group": d["group"],
            "subject": d["subject"],
            "n_dbpos_quantified": int(np.sum(~np.isnan(v) & (v > 0))),
            "dbpos_missing_rate": float(np.mean(np.isnan(v))),
            "n_species_quantified": int(np.sum(sv > 0)),
            "total_species_area": float(np.nansum(sv)),
            "median_dbpos_area": float(np.nanmedian(v[v > 0])) if np.any(v > 0) else np.nan,
        })
    qc = pd.DataFrame(rows)

    tot = qc["total_species_area"]
    log(f"样本 {len(qc)}；物质级总面积 中位 {tot.median():.3g}，"
        f"最小/最大 {tot.min():.3g}/{tot.max():.3g}（{tot.max()/tot.min():.1f} 倍差）")
    log(f"链×双键单元 定量到的数量：中位 {qc['n_dbpos_quantified'].median():.0f}，"
        f"范围 {qc['n_dbpos_quantified'].min()}~{qc['n_dbpos_quantified'].max()}")
    log(f"缺失率：中位 {qc['dbpos_missing_rate'].median():.1%}"
        "（缺失=该样本 MS2 没打到该单元，非过滤所致）")

    # 进样重复性：同一样品的 PB1 vs PB2 两针。定量表的列已经是合并后的样品，故这里必须
    # 回到未合并的逐针 *_PB_dbpos.csv 去比——这也正是评估“两针一致性”该用的粒度。
    rep = _injection_reproducibility(os.path.dirname(dbpos_csv))
    if len(rep):
        log("")
        log("进样重复性（同一样品 PB1 针 vs PB2 针，特征离子对强度 log10）：")
        log(f"  两针共有的单元数 中位 {rep['n_shared'].median():.0f}"
            f"（单针各测到 ~{rep['n_run1'].median():.0f} 个 -> 合并两针能显著扩覆盖）")
        log(f"  Pearson r  中位 {rep['pearson_log'].median():.3f}，"
            f"范围 {rep['pearson_log'].min():.3f}~{rep['pearson_log'].max():.3f}")
        log(f"  Spearman ρ 中位 {rep['spearman_log'].median():.3f}")
        log(f"  相对差 中位 {rep['median_rel_diff'].median():.1%}")

    if output_csv:
        qc.to_csv(output_csv, index=False, encoding="utf-8-sig")
        rep_csv = output_csv.replace(".csv", "_replicate.csv")
        rep.to_csv(rep_csv, index=False, encoding="utf-8-sig")
        log(f"质控表已保存 {output_csv} / {rep_csv}")
    return qc, rep


# ============================================================
# 差异分析
# ============================================================
def _subject_mean(df, value_cols, design):
    """把同一生物样本的 PB1/PB2 两次测量取均值 -> 每个 subject 一列。"""
    d = design.dropna(subset=["subject"])
    out = {}
    for subj, g in d.groupby("subject"):
        cols = [c for c in value_cols if any(c.startswith(s + "__") for s in g["sample"])]
        out[subj] = df[cols].mean(axis=1, skipna=True)
    return pd.DataFrame(out, index=df.index), d.drop_duplicates("subject").set_index("subject")["group"]


def _test_two_groups(mat, groups, ga=DEFAULT_CASE, gb=DEFAULT_CONTROL):
    """逐行做 ga vs gb（结果列的 _A/_B 后缀是位置：A 槽=ga，B 槽=gb）。

    主检验用 **Welch t**，另给 Mann-Whitney U 作稳健性对照。
    为什么不拿 U 当主检验：本设计 n_A=5 / n_B=6，U 检验的**最小可能 p = 2/C(11,5)
    = 0.0043**；在 ~900 个检验上做 BH 校正后最小 q ≈ 0.0043×900 ≈ 3.9，**数学上永远
    不可能 <0.05**——报 0 个显著是检验的功效上限所致，不是生物学结论。Welch t 的 p
    可以取到任意小，能真正区分。代价是 n=5/6 下正态性无法检验，故两个 p 都给出，
    结论以两者一致者为准。
    """
    ca = [c for c in mat.columns if groups.get(c) == ga]
    cb = [c for c in mat.columns if groups.get(c) == gb]
    res = []
    for i in range(len(mat)):
        a = mat.iloc[i][ca].to_numpy(dtype=float)
        b = mat.iloc[i][cb].to_numpy(dtype=float)
        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]
        na, nb = len(a), len(b)
        ma = float(np.mean(a)) if na else np.nan
        mb = float(np.mean(b)) if nb else np.nan
        if na < MIN_N or nb < MIN_N:
            res.append((na, nb, ma, mb, np.nan, np.nan, np.nan, f"n不足(A={na},B={nb})"))
            continue
        if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.isclose(a[0], b[0]):
            res.append((na, nb, ma, mb, np.nan, np.nan, np.nan, "两组取值完全相同"))
            continue
        p_t = stats.ttest_ind(a, b, equal_var=False).pvalue
        try:
            p_u = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            p_u = np.nan
        lfc = np.log2((ma + 1e-12) / (mb + 1e-12))
        res.append((na, nb, ma, mb, lfc, p_t, p_u, ""))
    return pd.DataFrame(res, index=mat.index,
                        columns=["n_A", "n_B", "mean_A", "mean_B", "log2FC",
                                 "p", "p_mannwhitney", "test_note"])


def differential(dbpos_csv, species_csv, output_prefix, alpha=ALPHA,
                 design=None, case=DEFAULT_CASE, control=DEFAULT_CONTROL):
    db = pd.read_csv(dbpos_csv)
    sp = pd.read_csv(species_csv)
    area_cols = [c for c in db.columns if c.endswith("__area")]
    design = parse_design(area_cols, design)
    used = design[design["group"].isin([case, control])]
    log("=" * 56)
    log("差异分析")
    log("=" * 56)
    log(f"比较 {case}(A 槽) vs {control}(B 槽)")
    log(f"进入分析的样本 {len(used)}（{used['group'].value_counts().to_dict()}），"
        f"同一生物样本的多次进样按均值合并")
    skipped = design[~design["group"].isin([case, control])]["sample"].tolist()
    if skipped:
        log(f"不在 {case}/{control} 设计里、已排除：{skipped}")
    for g in (case, control):
        n = int((used["group"] == g).sum())
        if n < MIN_N:
            log(f"⚠ 组 {g} 只有 {n} 个样本（<{MIN_N}），该组所有检验的 p 都会是 NaN")

    # ---- 第一层：物质总量（常规脂质组学能看到的）----
    sp_key = sp["species_key"]   # 由 pb_quantify 写入，物质峰唯一键
    sp_norm_cols = [c for c in sp.columns if c.endswith("__norm")]
    sp_mat, groups = _subject_mean(sp.set_index(sp_key), sp_norm_cols, design)
    sp_res = _test_two_groups(sp_mat, groups, case, control)
    sp_res["q"] = bh_fdr(sp_res["p"].to_numpy())
    sp_res.insert(0, "species_key", sp_res.index)
    sp_res.insert(1, "total_chain_name", sp.set_index(sp_key)["total_chain_name"].reindex(sp_res.index).values)
    sp_res.insert(2, "subclass", sp.set_index(sp_key)["subclass"].reindex(sp_res.index).values)
    sp_res = sp_res.reset_index(drop=True)
    n_sig_sp = int((sp_res["q"] < alpha).sum())
    log(f"[第一层 物质总量] 检验 {int(sp_res['p'].notna().sum())}/{len(sp_res)} 个物质峰，"
        f"FDR<{alpha} 显著 {n_sig_sp}")

    # ---- 第二层：双键位置组成（PB 独有）----
    # share = 该 (物质,链,位置) 面积 / 该物质该样本的总面积。与物质总量无关。
    sp_area_by_key = sp.set_index(sp_key)
    assert not sp_area_by_key.index.duplicated().any(), "species 键不唯一"
    db_key = db["species_key"]
    share = pd.DataFrame(index=db.index)
    for s in design["sample"]:
        num = db[f"{s}__area"].to_numpy(dtype=float)
        den = sp_area_by_key[f"{s}__area"].reindex(db_key).to_numpy(dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            share[f"{s}__share"] = np.where((den > 0), num / den, np.nan)
    share_cols = list(share.columns)
    sh_mat, _ = _subject_mean(share, share_cols, design)
    sh_res = _test_two_groups(sh_mat, groups, case, control)
    sh_res["q"] = bh_fdr(sh_res["p"].to_numpy())
    n_sig_sh = int((sh_res["q"] < alpha).sum())
    log(f"[第二层 位置份额] 检验 {int(sh_res['p'].notna().sum())}/{len(sh_res)} 个单元，"
        f"FDR<{alpha} 显著 {n_sig_sh}")

    out = db[["species_key", "subclass", "total_chain_name", "double_bond_position",
              "f", "consensus_rt", "n_samples_quantified"]].copy()
    out = out.join(sh_res.add_prefix("share_"))
    sp_lookup = sp_res.set_index("species_key")
    for c in ["log2FC", "p", "q", "n_A", "n_B"]:
        out[f"species_{c}"] = sp_lookup[c].reindex(out["species_key"]).values

    # ★ 四象限：这是本表最重要的一列
    sig_sp = out["species_q"] < alpha
    sig_sh = out["share_q"] < alpha
    out["quadrant"] = np.select(
        [sig_sh & ~sig_sp, sig_sh & sig_sp, ~sig_sh & sig_sp],
        ["位置独有(总量不变,位置份额变)", "量与位置同变", "仅总量变(常规即可)"],
        default="无变化",
    )
    out.loc[out["share_q"].isna() | out["species_q"].isna(), "quadrant"] = "未检验(n不足)"

    dbpos_out = f"{output_prefix}_dbposition_diff.csv"
    out.sort_values(["quadrant", "share_q"], kind="mergesort").to_csv(
        dbpos_out, index=False, encoding="utf-8-sig")
    sp_res.to_csv(f"{output_prefix}_species_diff.csv", index=False, encoding="utf-8-sig")

    log("")
    log("四象限分布（每个 物质×链×双键位置 单元）：")
    for k, v in out["quadrant"].value_counts().items():
        log(f"  {k:26} {v}")
    uniq = out[out["quadrant"] == "位置独有(总量不变,位置份额变)"]
    if len(uniq):
        log("")
        log(f"★ 只有 PB 能发现的 {len(uniq)} 个单元（总量无显著变化，但双键位置份额显著变化）：")
        for _, r in uniq.sort_values("share_q").head(15).iterrows():
            log(f"  {r['total_chain_name']:12} f={r['f']} {r['double_bond_position']:5} "
                f"份额 A={r['share_mean_A']:.3f} B={r['share_mean_B']:.3f} "
                f"q={r['share_q']:.3g} | 总量 log2FC={r['species_log2FC']:.2f} q={r['species_q']:.3g}")

    # ---- omega 家族指数（也是常规拿不到的）----
    omega = _omega_index(db, design, groups, alpha, case, control)
    omega_out = f"{output_prefix}_omega_index.csv"
    omega.to_csv(omega_out, index=False, encoding="utf-8-sig")
    n_sig_om = int((omega["q"] < alpha).sum())
    log("")
    log(f"[omega 家族指数] {len(omega)} 行（每物质的 n-3/n-6/n-9 份额及 n-3/n-6 比），"
        f"FDR<{alpha} 显著 {n_sig_om}")
    log(f"结果已保存 {dbpos_out} / {output_prefix}_species_diff.csv / {omega_out}")
    return out, sp_res, omega


def _omega_index(db, design, groups, alpha, case=DEFAULT_CASE, control=DEFAULT_CONTROL):
    """每物质每样本的 n-3 / n-6 / n-9 份额之和，再做 case vs control。"""
    pos = db["double_bond_position"].str.extract(r"n-(\d+)")[0].astype(float)
    recs = []
    for fam in [3, 6, 9]:
        sel = pos == fam
        if not sel.any():
            continue
        for key, g in db[sel].groupby("total_chain_name"):
            row = {"total_chain_name": key, "metric": f"n-{fam} 份额之和"}
            for s in design["sample"]:
                row[s] = np.nansum(g[f"{s}__area"].to_numpy(dtype=float))
            recs.append(row)
    if not recs:
        return pd.DataFrame()
    W = pd.DataFrame(recs)
    val_cols = [c for c in W.columns if c not in ("total_chain_name", "metric")]
    # 归一到该物质在该样本的 n-3/6/9 总和，得组成比
    denom = W.groupby("total_chain_name")[val_cols].transform("sum")
    frac = W[val_cols] / denom.replace(0, np.nan)
    frac.columns = [f"{c}__share" for c in val_cols]
    mat, _ = _subject_mean(frac, list(frac.columns), design)
    res = _test_two_groups(mat, groups, case, control)
    res["q"] = bh_fdr(res["p"].to_numpy())
    return pd.concat([W[["total_chain_name", "metric"]].reset_index(drop=True),
                      res.reset_index(drop=True)], axis=1)


DESIGN_JSON = "PB_design.json"


def available_samples(pb_output_dir):
    """定量表里有哪些样品（列名去掉 __area）。给 agent 列"能用哪些文件"用。"""
    dbpos_csv = os.path.join(pb_output_dir, "PB_dbposition_quant.csv")
    cols = pd.read_csv(dbpos_csv, nrows=0).columns
    return [c[:-len("__area")] for c in cols if c.endswith("__area")]


def load_design(pb_output_dir):
    """读回上次存的设计；没有就返回 None（调用方退回样本名前缀约定）。"""
    path = os.path.join(pb_output_dir, DESIGN_JSON)
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_design(pb_output_dir, design, case, control, note=""):
    """把设计连同 A/B 槽对应的真实组名存盘，报告层据此把 _A/_B 翻译回去。"""
    path = os.path.join(pb_output_dir, DESIGN_JSON)
    payload = {"groups": design, "case": case, "control": control, "note": note,
               "slot_A": case, "slot_B": control}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    return path


def run_pb_qc_diff(pb_output_dir, design=None, case=None, control=None):
    """质控 + 差异分析。

    design/case/control 为 None 时退回历史行为（样本名前缀，A vs B），故 run.py、
    命令行、老脚本的调用不受影响。传了就按传的来，并把设计存进 PB_design.json。
    """
    dbpos_csv = os.path.join(pb_output_dir, "PB_dbposition_quant.csv")
    species_csv = os.path.join(pb_output_dir, "PB_species_quant.csv")
    prefix = os.path.join(pb_output_dir, "PB")
    case = case or DEFAULT_CASE
    control = control or DEFAULT_CONTROL
    if design:
        save_design(pb_output_dir, design, case, control)
    qc_report(dbpos_csv, species_csv, os.path.join(pb_output_dir, "PB_qc.csv"), design=design)
    diff = differential(dbpos_csv, species_csv, prefix, design=design, case=case, control=control)
    pref = position_preference(dbpos_csv, prefix, design=design, case=case, control=control)
    return diff, pref


def position_preference(dbpos_csv, output_prefix, min_group_n=2,
                        design=None, case=DEFAULT_CASE, control=DEFAULT_CONTROL):
    """同一 (物质, f) 内部各双键位置的构成谱，逐组给出偏向。**份额同 f 之和 = 1**。

    为什么按 (物质, f) 分组，这是本分析成立的前提：
      f = 该双键**甲基侧还有几个双键**。在给定 f 层，**每个分子有且只有一个双键**——
      有的分子它落在 n-6，有的落在 n-9。所以同 f 内各位置的份额是对这批分子的**真实
      划分**，天然加总为 1，“A 偏向 n-9 / B 偏向 n-7”这种话只在这里成立。
      跨 f 比份额没有这个含义(同一批分子被数了多次)，故不跨 f 归一。
      不分链：诊断离子质量只由 (n,f) 与母离子决定，与链如何切分无关(见 pb_localize)。

    份额怎么算才加总为 1：
      在**某样本已经看到该 (物质,f) 层的前提下**，该层没检出的位置意味着其诊断离子对
      在这张 MS2 谱里低于检出——对“构成”而言就是 ~0，不是“未知”。故：
        - 该样本对该 (物质,f) 有任何检出 -> 未检出的位置记 0，该样本各位置和 = 1；
        - 该样本对该物质完全没打到     -> 该样本整组 NaN，不参与。
      先按样本归一再按组求均值时，若把未检出记 NaN，各样本的位置集合不同，组均值就不
      加总为 1(实测出现过 n-6:83% + n-3:52% = 135% 这种读不通的结果)。

    不做假设检验——n=5/6 功效太低。这里给效应量：每组各位置平均份额、各组优势位置、
    优势位置是否在两组间翻转(flip)。
    """
    db = pd.read_csv(dbpos_csv)
    area_cols = [c for c in db.columns if c.endswith("__area")]
    design = parse_design(area_cols, design)
    used = design[design["group"].isin([case, control])]
    gmap = used.drop_duplicates("subject").set_index("subject")["group"]

    rows = []
    for (skey, species, f), g in db.groupby(["species_key", "total_chain_name", "f"], sort=False):
        if len(g) < 2:
            continue                      # 该 f 层只检出 1 个位置 -> 无“偏向”可谈
        prof = {}
        for _, d in used.iterrows():
            v = g[f"{d['sample']}__area"].to_numpy(dtype=float)
            if np.all(np.isnan(v)):
                prof[d["sample"]] = np.full(len(v), np.nan)   # 该样本没打到这一层
                continue
            v = np.nan_to_num(v, nan=0.0)                     # 打到了但该位置没检出 -> ~0
            tot = v.sum()
            prof[d["sample"]] = v / tot if tot > 0 else np.full(len(v), np.nan)
        P = pd.DataFrame(prof, index=g["double_bond_position"].values)
        subj = {sj: P[list(sg["sample"])].mean(axis=1, skipna=True)
                for sj, sg in used.groupby("subject")}
        S = pd.DataFrame(subj)
        ca = [c for c in S.columns if gmap.get(c) == case]
        cb = [c for c in S.columns if gmap.get(c) == control]
        mean_a, mean_b = S[ca].mean(axis=1), S[cb].mean(axis=1)
        n_a, n_b = S[ca].notna().sum(axis=1), S[cb].notna().sum(axis=1)
        if n_a.max() < min_group_n or n_b.max() < min_group_n:
            continue
        dom_a = mean_a.idxmax() if mean_a.notna().any() else None
        dom_b = mean_b.idxmax() if mean_b.notna().any() else None
        for pos in P.index:
            rows.append({
                "species_key": skey, "total_chain_name": species, "f": f,
                "double_bond_position": pos, "n_competing_positions": len(g),
                "share_A": mean_a.get(pos, np.nan), "share_B": mean_b.get(pos, np.nan),
                "delta_share_A_minus_B": mean_a.get(pos, np.nan) - mean_b.get(pos, np.nan),
                "n_A": int(n_a.get(pos, 0)), "n_B": int(n_b.get(pos, 0)),
                "dominant_A": dom_a, "dominant_B": dom_b,
                "preference_flip": (dom_a is not None and dom_b is not None and dom_a != dom_b),
            })
    prof = pd.DataFrame(rows)
    if prof.empty:
        log("没有任何 (物质,f) 组同时检出 ≥2 个竞争位置且两组都有足够样本")
        return prof, pd.DataFrame()

    long_csv = f"{output_prefix}_position_profile.csv"
    prof.to_csv(long_csv, index=False, encoding="utf-8-sig")

    summ = (prof.groupby(["species_key", "total_chain_name", "f"], sort=False)
            .agg(positions=("double_bond_position", lambda x: ",".join(x)),
                 n_positions=("double_bond_position", "count"),
                 dominant_A=("dominant_A", "first"), dominant_B=("dominant_B", "first"),
                 preference_flip=("preference_flip", "first"),
                 share_sum_A=("share_A", "sum"), share_sum_B=("share_B", "sum"),
                 max_abs_delta=("delta_share_A_minus_B", lambda x: np.nanmax(np.abs(x))))
            .reset_index())
    summ = summ.sort_values(["preference_flip", "max_abs_delta"], ascending=[False, False])
    summ_csv = f"{output_prefix}_position_preference.csv"
    summ.to_csv(summ_csv, index=False, encoding="utf-8-sig")

    log("=" * 56)
    log("位置偏好谱（同一物质、同一 f 层内部，各双键位置的构成）")
    log("=" * 56)
    bad = summ[(summ["share_sum_A"] - 1).abs() > 0.01]
    log(f"份额归一自检：{len(summ)} 组中 A 组份额和偏离 1 超过 1% 的 {len(bad)} 组"
        f"（应为 0）")
    log(f"可比较的 (物质,f) 组 {len(summ)} 个（该 f 层检出 ≥2 个竞争位置）")
    n_flip = int(summ["preference_flip"].sum())
    log(f"其中**优势位置在 A/B 间翻转**的 {n_flip} 个 —— 这类最能体现 PB 的价值：")
    log("  (常规脂质组学只看到同一个物质，完全看不到内部的位置构成变了)")
    log("")
    for _, r in summ[summ["preference_flip"]].head(12).iterrows():
        sub = prof[(prof.species_key == r["species_key"]) & (prof.f == r["f"])]
        a_txt = "; ".join(f"{x.double_bond_position}:{x.share_A:.0%}"
                          for x in sub.itertuples() if pd.notna(x.share_A))
        b_txt = "; ".join(f"{x.double_bond_position}:{x.share_B:.0%}"
                          for x in sub.itertuples() if pd.notna(x.share_B))
        log(f"  {r['total_chain_name']:12} f={r['f']}  A 偏向 {r['dominant_A']:5} | "
            f"B 偏向 {r['dominant_B']:5}  (最大份额差 {r['max_abs_delta']:.0%})")
        log(f"      A: {a_txt}")
        log(f"      B: {b_txt}")
    log("")
    log(f"位置构成谱已保存 {long_csv}")
    log(f"偏好摘要已保存 {summ_csv}")
    return prof, summ


PB_OUT_DIR = r"D:\bio_inf\LipidIN_github\PB_demo\sample\PB\pkl"


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass
    run_pb_qc_diff(PB_OUT_DIR)
