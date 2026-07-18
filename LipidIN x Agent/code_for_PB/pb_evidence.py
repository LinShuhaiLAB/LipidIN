# -*- coding: utf-8 -*-
"""
PB 结果的证据层：把盘上的 CSV 读成"能喂给大模型"和"能打分"的东西。

这一层**不调用大模型**，只做确定性的读表、摘要、切片和打分。大模型相关的都在
`pb_llm_analysis.py` / `pb_verify.py` 里。这样分开的好处：假设是否成立最终由数字说了算，
数字怎么算的可以单独测、不受模型输出漂移影响。

三种可信度（最终报告里每条结论都要带）：
  注释可信度 —— 这个脂质到底注释得靠不靠谱。来自 PB_dbposition_quant 的
      master_match_status(是否与常规注释主表对上) / master_evidence_type(MS2 还是仅 MS1)
      / n_samples_ms2(多少样本真打到了 MS2)。
  强度可信度 —— 信号强不强、测得稳不稳。来自 max_pair_intensity(诊断离子对强度)、
      n_samples_quantified(多少样本定量到)、median_pair_share(离子对占比)，
      以及 PB_qc_replicate 的两针相关性(整体重复性)。
  是否够新 —— 文献里有没有人报道过。由 pb_llm_analysis 检索文献后按命中数给分：
      没人报过=新(高分)，报过很多=不新(低分)。注意"新"不等于"对"，故在报告里与
      前两项分开列，不合并成一个总分。
"""

import os
import re

import numpy as np
import pandas as pd

# 四象限的字面值（pb_qc_diff 写出来就是中文字符串，别改，改了对不上）
Q_POSITION_ONLY = "位置独有(总量不变,位置份额变)"
Q_BOTH = "量与位置同变"
Q_SPECIES_ONLY = "仅总量变(常规即可)"
Q_NONE = "无变化"
Q_UNTESTED = "未检验(n不足)"

MASTER_MATCHED = "master匹配"
EVIDENCE_MS2 = "MS2注释"

ALPHA = 0.05

TABLE_FILES = {
    "dbpos_quant": "PB_dbposition_quant.csv",
    "species_quant": "PB_species_quant.csv",
    "dbpos_diff": "PB_dbposition_diff.csv",
    "species_diff": "PB_species_diff.csv",
    "omega": "PB_omega_index.csv",
    "qc": "PB_qc.csv",
    "qc_replicate": "PB_qc_replicate.csv",
    "position_profile": "PB_position_profile.csv",
    "position_preference": "PB_position_preference.csv",
}


def load_tables(pb_dir):
    """读入所有 PB 产物表。缺哪张就没有哪个键——调用方自己判断够不够跑。"""
    out = {}
    for key, fname in TABLE_FILES.items():
        path = os.path.join(pb_dir, fname)
        if os.path.exists(path):
            try:
                out[key] = pd.read_csv(path, encoding="utf-8-sig")
            except Exception:
                pass
    return out


def missing_tables(tables, needed):
    return [TABLE_FILES[k] for k in needed if k not in tables or tables[k] is None]


# ------------------------------------------------------------------ #
# 摘要：喂给大模型的"这批数据发生了什么"
# ------------------------------------------------------------------ #
def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def design_label(design_json):
    """A 槽/B 槽对应的真实组名。没有设计文件就是默认的 A/B。"""
    if not design_json:
        return "A", "B"
    return design_json.get("case") or "A", design_json.get("control") or "B"


def qc_digest(tables, design_json=None):
    """质控摘要：样本数、缺失率、总面积离散、两针重复性。"""
    lines = []
    qc = tables.get("qc")
    if qc is not None and len(qc):
        used = qc[qc["group"].notna()]
        counts = used["group"].value_counts().to_dict()
        lines.append(f"样本 {len(qc)} 个，进入比较的分组：{counts}；"
                     f"不分组(仅质控)：{qc[qc['group'].isna()]['sample'].tolist()}")
        tot = qc["total_species_area"]
        lines.append(f"物质级总面积：中位 {_fmt(tot.median())}，"
                     f"最小/最大 {_fmt(tot.min())}/{_fmt(tot.max())}"
                     f"（{_fmt(tot.max()/tot.min() if tot.min() else np.nan, 2)} 倍差）")
        lines.append(f"链×双键单元定量数：中位 {_fmt(qc['n_dbpos_quantified'].median(), 4)}，"
                     f"范围 {int(qc['n_dbpos_quantified'].min())}~{int(qc['n_dbpos_quantified'].max())}")
        lines.append(f"缺失率中位 {qc['dbpos_missing_rate'].median():.1%}"
                     f"（缺失=该样本 MS2 没打到该单元，不是过滤掉的）")
    rep = tables.get("qc_replicate")
    if rep is not None and len(rep):
        lines.append(f"两针重复性（同一样品 PB1 vs PB2，log10 强度）：Pearson r 中位 "
                     f"{_fmt(rep['pearson_log'].median())}，相对差中位 "
                     f"{rep['median_rel_diff'].median():.1%}，两针共有单元中位 "
                     f"{_fmt(rep['n_shared'].median(), 4)}")
    return "\n".join(lines) if lines else "（无质控表）"


def diff_digest(tables, design_json=None, top_n=18):
    """差异摘要。重点给 PB 独有的那一层，常规能做的那层只给一句话。"""
    case, control = design_label(design_json)
    lines = []
    sp = tables.get("species_diff")
    if sp is not None and len(sp):
        n_sig = int((sp["q"] < ALPHA).sum())
        lines.append(f"[第一层 物质总量：常规脂质组学也能拿到] 检验 "
                     f"{int(sp['p'].notna().sum())}/{len(sp)} 个物质峰，FDR<{ALPHA} 显著 {n_sig}")
    db = tables.get("dbpos_diff")
    if db is not None and len(db):
        n_sig_sh = int((db["share_q"] < ALPHA).sum())
        lines.append(f"[第二层 双键位置份额：PB 独有、常规完全看不见] 检验 "
                     f"{int(db['share_p'].notna().sum())}/{len(db)} 个 物质×链×位置 单元，"
                     f"FDR<{ALPHA} 显著 {n_sig_sh}")
        lines.append("四象限分布：" + str(db["quadrant"].value_counts().to_dict()))

        uniq = db[db["quadrant"] == Q_POSITION_ONLY].sort_values("share_q")
        if len(uniq):
            lines.append("")
            lines.append(f"★ 只有 PB 能发现的单元（总量无显著变化，但双键位置份额显著变化），"
                         f"共 {len(uniq)} 个，最显著的 {min(top_n, len(uniq))} 个：")
            for _, r in uniq.head(top_n).iterrows():
                lines.append(
                    f"  {r['total_chain_name']} (f={r['f']}) {r['double_bond_position']}：份额 "
                    f"{case}={_fmt(r['share_mean_A'])} {control}={_fmt(r['share_mean_B'])}，"
                    f"share_q={_fmt(r['share_q'])}；同物质总量 log2FC={_fmt(r['species_log2FC'])} "
                    f"q={_fmt(r['species_q'])}")
        both = db[db["quadrant"] == Q_BOTH].sort_values("share_q")
        if len(both):
            lines.append("")
            lines.append(f"量与位置同变的单元 {len(both)} 个，最显著的 {min(8, len(both))} 个：")
            for _, r in both.head(8).iterrows():
                lines.append(
                    f"  {r['total_chain_name']} (f={r['f']}) {r['double_bond_position']}："
                    f"份额 {case}={_fmt(r['share_mean_A'])} {control}={_fmt(r['share_mean_B'])} "
                    f"share_q={_fmt(r['share_q'])}；总量 log2FC={_fmt(r['species_log2FC'])} "
                    f"q={_fmt(r['species_q'])}")
    return "\n".join(lines) if lines else "（无差异表）"


def preference_digest(tables, design_json=None, top_n=12):
    """位置偏好翻转：同一物质同一 f 层里，优势双键位置在两组间换人了。"""
    case, control = design_label(design_json)
    pref = tables.get("position_preference")
    prof = tables.get("position_profile")
    if pref is None or not len(pref):
        return "（无位置偏好表）"
    flips = pref[pref["preference_flip"] == True]  # noqa: E712  (csv 读回来是 bool/str)
    lines = [f"可比较的 (物质,f) 层 {len(pref)} 个（该层检出 ≥2 个竞争位置），"
             f"其中优势位置在 {case}/{control} 间**翻转**的 {len(flips)} 个。"]
    if len(flips):
        lines.append(f"翻转最强的 {min(top_n, len(flips))} 个（构成谱，同 f 层各位置份额和=1）：")
        for _, r in flips.head(top_n).iterrows():
            txt = ""
            if prof is not None:
                sub = prof[(prof["species_key"] == r["species_key"]) & (prof["f"] == r["f"])]
                a_txt = "; ".join(f"{x.double_bond_position}:{x.share_A:.0%}"
                                  for x in sub.itertuples() if pd.notna(x.share_A))
                b_txt = "; ".join(f"{x.double_bond_position}:{x.share_B:.0%}"
                                  for x in sub.itertuples() if pd.notna(x.share_B))
                txt = f"\n      {case}: {a_txt}\n      {control}: {b_txt}"
            lines.append(f"  {r['total_chain_name']} f={r['f']}：{case} 偏向 {r['dominant_A']}，"
                         f"{control} 偏向 {r['dominant_B']}（最大份额差 {r['max_abs_delta']:.0%}）{txt}")
    return "\n".join(lines)


def omega_digest(tables, design_json=None, top_n=10):
    """omega 家族指数（n-3/n-6/n-9 份额），脂质生物学常用、常规方法拿不到。"""
    case, control = design_label(design_json)
    om = tables.get("omega")
    if om is None or not len(om):
        return "（无 omega 表）"
    sig = om[om["q"] < ALPHA].sort_values("q")
    lines = [f"omega 家族指数 {len(om)} 行（每物质的 n-3/n-6/n-9 份额之和），"
             f"FDR<{ALPHA} 显著 {len(sig)}。"]
    for _, r in sig.head(top_n).iterrows():
        lines.append(f"  {r['total_chain_name']} {r['metric']}：{case}={_fmt(r['mean_A'])} "
                     f"{control}={_fmt(r['mean_B'])} log2FC={_fmt(r['log2FC'])} q={_fmt(r['q'])}")
    return "\n".join(lines)


def power_reality(tables):
    """这批数据的统计现实——**必须**原样喂给大模型，否则它会编出不存在的显著性。

    实测 n=5/6 时 FDR 校正后常常一个显著都没有（这是功效上限，不是"没有生物学"）。
    这种情况下唯一诚实的做法是：结论只能停在**效应量**层面（位置偏好翻转、份额差），
    并明说是探索性的、需要扩样本验证。
    """
    bits = []
    sp = tables.get("species_diff")
    db = tables.get("dbpos_diff")
    om = tables.get("omega")
    n_sig_sp = int((sp["q"] < ALPHA).sum()) if sp is not None and len(sp) else 0
    n_sig_sh = int((db["share_q"] < ALPHA).sum()) if db is not None and len(db) else 0
    n_sig_om = int((om["q"] < ALPHA).sum()) if om is not None and len(om) else 0
    qc = tables.get("qc")
    if qc is not None and len(qc):
        counts = qc[qc["group"].notna()]["group"].value_counts().to_dict()
        bits.append(f"每组样本数：{counts}。")
    if n_sig_sp == 0 and n_sig_sh == 0 and n_sig_om == 0:
        bits.append(
            "⚠ 关键：**三个层面 FDR<0.05 的显著项都是 0**。这是小样本的功效上限，"
            "不代表没有生物学差异。因此：\n"
            "  - 不许把任何东西说成'显著'、'统计学上确认'；\n"
            "  - 结论只能建立在**效应量**上（位置偏好翻转、份额差、omega 份额差）；\n"
            "  - 每条假设都必须标明是探索性的、需要更大样本验证；\n"
            "  - 若某假设在数据里连效应量都很弱，就不要提这条假设。")
    else:
        bits.append(f"FDR<{ALPHA} 显著项：物质总量层 {n_sig_sp}，位置份额层 {n_sig_sh}，"
                    f"omega 层 {n_sig_om}。未达显著的部分只能当效应量看。")
    pref = tables.get("position_preference")
    if pref is not None and len(pref):
        n_flip = int((pref["preference_flip"] == True).sum())  # noqa: E712
        bits.append(f"位置偏好翻转 {n_flip} 个（这一层按设计就不做假设检验，只给效应量）。")
    return "\n".join(bits)


def full_digest(tables, design_json=None):
    """喂给大模型的完整数据摘要。"""
    case, control = design_label(design_json)
    parts = [
        f"### 比较：{case}（A 槽）vs {control}（B 槽）。log2FC>0 表示在 {case} 中更高。",
        "",
        "### 统计现实（先读这一段，它决定了你能说什么）", power_reality(tables), "",
        "### 质控", qc_digest(tables, design_json), "",
        "### 差异分析（两层）", diff_digest(tables, design_json), "",
        "### 双键位置偏好", preference_digest(tables, design_json), "",
        "### omega 家族指数", omega_digest(tables, design_json),
    ]
    return "\n".join(parts)


# ------------------------------------------------------------------ #
# 切片：某条假设涉及哪些数据行
# ------------------------------------------------------------------ #
def _norm_name(x):
    return re.sub(r"\s+", " ", str(x or "")).strip().lower()


def match_rows(df, lipids=None, classes=None, positions=None, name_col="total_chain_name"):
    """按脂质名 / 类别 / 双键位置挑出相关行。都没给就返回空表——宁可说"没证据"。

    这里的取舍很关键，直接决定证据切片和三个可信度算得准不准：

    具体脂质优先，**不与类别取并集**。一条"PC 40:5 的双键位置重排"的假设，
    involved_classes 通常也带着 "PC"；若并集，就会把全部 PC 的行（实测 768 行）都算进来，
    证据被稀释成整个 PC 类的平均水平，这条假设自己的强弱反而看不见了。所以给了具体
    脂质就只看这些脂质，没给才退到类别。

    位置是**过滤**（AND）不是并集：并集会把所有脂质里出现 n-9 的行都拉进来，那不是这条
    假设的证据。但过滤后若一行不剩（模型给的位置在该脂质上不存在），就退回不带位置的
    结果——保留证据让子 agent 自己看出"位置对不上"，比直接判 0 匹配更接近真相。
    """
    if df is None or not len(df):
        return pd.DataFrame()
    if lipids and name_col in df.columns:
        want = {_norm_name(x) for x in lipids}
        mask = df[name_col].map(lambda v: _norm_name(v) in want)
    elif classes and "subclass" in df.columns:
        want = {_norm_name(x) for x in classes}
        mask = df["subclass"].map(lambda v: _norm_name(v) in want)
    else:
        return pd.DataFrame()
    base = df[mask]
    if positions and "double_bond_position" in base.columns and len(base):
        want = {_norm_name(x) for x in positions}
        narrowed = base[base["double_bond_position"].map(lambda v: _norm_name(v) in want)]
        if len(narrowed):
            return narrowed
    return base


def evidence_slice(tables, hypothesis, design_json=None, max_rows=25):
    """把某条假设涉及的**真实数据行**渲染成文本，交给验证子 agent 自己判断。

    这是"数据驱动验证"的关键：子 agent 看到的是这条假设名下的具体数字（份额、q 值、
    象限、omega），而不是假设的复述，所以它能发现"假设说的和数据不符"。
    """
    case, control = design_label(design_json)
    lipids = hypothesis.get("involved_lipids") or []
    classes = hypothesis.get("involved_classes") or []
    positions = hypothesis.get("involved_positions") or []
    lines = []

    db = match_rows(tables.get("dbpos_diff"), lipids, classes, positions)
    lines.append(f"【双键位置层 差异行】匹配到 {len(db)} 行"
                 + ("（这条假设在数据里找不到对应的位置层证据行）" if not len(db) else ""))
    if len(db):
        cnt = db["quadrant"].value_counts().to_dict()
        lines.append(f"  象限分布：{cnt}")
        for _, r in db.sort_values("share_q").head(max_rows).iterrows():
            lines.append(
                f"  {r['total_chain_name']} (f={r['f']}) {r['double_bond_position']}｜"
                f"份额 {case}={_fmt(r['share_mean_A'])} {control}={_fmt(r['share_mean_B'])}｜"
                f"share_log2FC={_fmt(r['share_log2FC'])} share_q={_fmt(r['share_q'])}"
                f"（n={_fmt(r['share_n_A'],2)}/{_fmt(r['share_n_B'],2)}）｜"
                f"总量 log2FC={_fmt(r['species_log2FC'])} q={_fmt(r['species_q'])}｜{r['quadrant']}")

    sp = match_rows(tables.get("species_diff"), lipids, classes)
    lines.append(f"\n【物质总量层 差异行】匹配到 {len(sp)} 行")
    for _, r in sp.sort_values("q").head(12).iterrows():
        lines.append(f"  {r['total_chain_name']}｜{case}={_fmt(r['mean_A'])} "
                     f"{control}={_fmt(r['mean_B'])}｜log2FC={_fmt(r['log2FC'])} "
                     f"q={_fmt(r['q'])} p={_fmt(r['p'])}｜{r['test_note'] or '已检验'}")

    pref = match_rows(tables.get("position_preference"), lipids, classes)
    if len(pref):
        lines.append(f"\n【位置偏好层】匹配到 {len(pref)} 个 (物质,f) 层，"
                     f"其中翻转 {int((pref['preference_flip'] == True).sum())} 个")  # noqa: E712
        for _, r in pref.head(10).iterrows():
            lines.append(f"  {r['total_chain_name']} f={r['f']}：{case} 偏向 {r['dominant_A']}，"
                         f"{control} 偏向 {r['dominant_B']}，最大份额差 "
                         f"{_fmt(r['max_abs_delta'])}，翻转={r['preference_flip']}")

    om = match_rows(tables.get("omega"), lipids, classes)
    if len(om):
        lines.append(f"\n【omega 指数层】匹配到 {len(om)} 行")
        for _, r in om.sort_values("q").head(8).iterrows():
            lines.append(f"  {r['total_chain_name']} {r['metric']}：{case}={_fmt(r['mean_A'])} "
                         f"{control}={_fmt(r['mean_B'])} q={_fmt(r['q'])}")
    return "\n".join(lines)


# ------------------------------------------------------------------ #
# 三种可信度
# ------------------------------------------------------------------ #
def _level(score, lang="zh"):
    if score >= 0.7:
        return "高" if lang == "zh" else "high"
    if score >= 0.4:
        return "中" if lang == "zh" else "medium"
    return "低" if lang == "zh" else "low"


def annotation_confidence(tables, hypothesis, lang="zh"):
    """注释可信度：这些脂质注释得靠不靠谱。

    master匹配 = 与常规注释主表对上了（不是 PB 这边孤零零冒出来的峰）；
    MS2注释 = 有二级证据，不是只有一级质量数；
    n_samples_ms2 = 多少样本真打到了 MS2（打到的样本越多越不像偶然）。
    """
    q = tables.get("dbpos_quant")
    rows = match_rows(q, hypothesis.get("involved_lipids"), hypothesis.get("involved_classes"),
                      hypothesis.get("involved_positions"))
    if not len(rows):
        return {"n": 0, "score": 0.0, "level": _level(0, lang),
                "detail": "没有匹配到定量表里的行 / no matching rows"}
    n_samples = len([c for c in q.columns if c.endswith("__area")]) or 1
    master_frac = float((rows["master_match_status"] == MASTER_MATCHED).mean())
    ms2_frac = float((rows["master_evidence_type"] == EVIDENCE_MS2).mean())
    ms2_cov = float(np.clip(rows["n_samples_ms2"].mean() / n_samples, 0, 1))
    score = round(0.40 * master_frac + 0.35 * ms2_frac + 0.25 * ms2_cov, 3)
    return {
        "n": int(len(rows)), "score": score, "level": _level(score, lang),
        "master_matched_fraction": round(master_frac, 3),
        "ms2_annotated_fraction": round(ms2_frac, 3),
        "ms2_sample_coverage": round(ms2_cov, 3),
        "detail": (f"{len(rows)} 个单元：master 匹配 {master_frac:.0%}，MS2 注释 {ms2_frac:.0%}，"
                   f"MS2 样本覆盖 {ms2_cov:.0%}"),
    }


def intensity_confidence(tables, hypothesis, lang="zh"):
    """强度可信度：信号强不强、测得稳不稳。

    强度用**分位数**而不是绝对值——绝对强度没有可比的尺子，但"在本批数据里排前 20%"
    是有意义的。重复性用两针 Pearson r 的中位数（整批一个值，不分假设）。
    """
    q = tables.get("dbpos_quant")
    rows = match_rows(q, hypothesis.get("involved_lipids"), hypothesis.get("involved_classes"),
                      hypothesis.get("involved_positions"))
    if not len(rows) or q is None:
        return {"n": 0, "score": 0.0, "level": _level(0, lang),
                "detail": "没有匹配到定量表里的行 / no matching rows"}
    n_samples = len([c for c in q.columns if c.endswith("__area")]) or 1
    # 强度分位：这些单元的诊断离子对强度在全体单元里排第几
    all_int = q["max_pair_intensity"].to_numpy(dtype=float)
    all_int = all_int[np.isfinite(all_int)]
    med_int = float(np.nanmedian(rows["max_pair_intensity"].to_numpy(dtype=float)))
    pct = float(np.mean(all_int <= med_int)) if len(all_int) and np.isfinite(med_int) else 0.0
    quant_cov = float(np.clip(rows["n_samples_quantified"].mean() / n_samples, 0, 1))
    pair_share = float(np.clip(np.nanmedian(rows["median_pair_share"].to_numpy(dtype=float)), 0, 1))
    rep = tables.get("qc_replicate")
    rep_r = float(rep["pearson_log"].median()) if rep is not None and len(rep) else np.nan
    rep_term = float(np.clip(rep_r, 0, 1)) if np.isfinite(rep_r) else 0.5
    score = round(0.35 * pct + 0.30 * quant_cov + 0.15 * pair_share + 0.20 * rep_term, 3)
    return {
        "n": int(len(rows)), "score": score, "level": _level(score, lang),
        "intensity_percentile": round(pct, 3),
        "median_max_pair_intensity": med_int,
        "quantified_sample_coverage": round(quant_cov, 3),
        "median_pair_share": round(pair_share, 3),
        "replicate_pearson_log": None if not np.isfinite(rep_r) else round(rep_r, 3),
        "detail": (f"强度处于全体单元的 {pct:.0%} 分位，定量样本覆盖 {quant_cov:.0%}，"
                   f"离子对占比中位 {pair_share:.0%}，两针相关 r="
                   f"{'NA' if not np.isfinite(rep_r) else f'{rep_r:.2f}'}"),
    }


def novelty_score(references, lang="zh"):
    """是否够新：文献里有没有人报道过。

    没检索到 = 可能新，但也可能是检索词没写好——所以 detail 里必须把检索到几篇讲清楚，
    让人自己判断，不要把"没搜到"直接当成"世界首报"。
    """
    n = len(references or [])
    if n == 0:
        score, verdict = 0.85, ("未检索到直接相关报道，可能是新发现（也可能是检索词太窄）"
                                if lang == "zh" else "no direct report found; possibly novel")
    elif n <= 2:
        score, verdict = 0.6, (f"检索到 {n} 篇相关报道，属于少有人碰的方向"
                               if lang == "zh" else f"{n} related report(s); lightly studied")
    elif n <= 5:
        score, verdict = 0.4, (f"检索到 {n} 篇相关报道，已有一定研究"
                               if lang == "zh" else f"{n} related reports; some prior work")
    else:
        score, verdict = 0.2, (f"检索到 {n} 篇相关报道，是研究较多的方向"
                               if lang == "zh" else f"{n} related reports; well studied")
    return {"n_references": n, "score": score, "is_novel": n <= 2, "detail": verdict,
            "references": list(references or [])[:5]}


def pathway_frame(tables, layer="pb"):
    """把 PB 的差异结果整理成 lipid_pathways 认识的形状。

    lipid_pathways.pathway_enrichment 要 subclass / significant / log2fc / direction 四列。
    layer='pb'      -> 用位置份额层的显著性（PB 独有的那层，本管线的核心）
    layer='species' -> 用物质总量层的显著性（常规脂质组学那层，用来对照）
    """
    if layer == "species":
        sp = tables.get("species_diff")
        if sp is None or not len(sp):
            return pd.DataFrame()
        out = pd.DataFrame({
            "subclass": sp["subclass"],
            "log2fc": sp["log2FC"],
            "significant": sp["q"] < ALPHA,
        })
    else:
        db = tables.get("dbpos_diff")
        if db is None or not len(db):
            return pd.DataFrame()
        out = pd.DataFrame({
            "subclass": db["subclass"],
            "log2fc": db["share_log2FC"],
            "significant": db["quadrant"].isin([Q_POSITION_ONLY, Q_BOTH]),
        })
    out["direction"] = np.where(out["significant"] & (out["log2fc"] > 0), "up",
                                np.where(out["significant"] & (out["log2fc"] <= 0), "down", "ns"))
    return out.dropna(subset=["subclass"])
