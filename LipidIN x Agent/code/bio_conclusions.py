"""Generate biological conclusions from the differential analysis (LLM-driven, offline fallback).

从差异分析生成生物学结论（大模型驱动，含离线回退）。

If an LLM is configured, it proposes N conclusions as structured JSON from a compact evidence digest.
Otherwise a deterministic rule engine derives conclusions from the most-changed lipid classes, the
unsaturation/chain trend and the random-forest separability.

若配置了大模型，则从紧凑的证据摘要中产出 N 条结构化 JSON 结论；否则用确定性规则引擎，从变化最大的
脂质类别、不饱和度/链长趋势与随机森林可分性中推导结论。
"""

import numpy as np

import lipid_pathways


def _digest(analysis, annotation_overview, sample_context=None):
    """Build a compact, model-friendly evidence digest string from the analysis result.

    构建紧凑、适合模型阅读的证据摘要字符串。
    """
    lines = []
    if sample_context:
        lines.append(f"Study context: {sample_context}")
    lines.append(f"Comparison: {analysis['case']} (case) vs {analysis['control']} (control).")
    lines.append(f"Features tested: {analysis['n_features']}, significant: {analysis['n_significant']} "
                 f"(up {analysis['n_up']}, down {analysis['n_down']}).")
    up = ", ".join(f"{c['name']}({c['subclass']},log2FC {c['log2fc']})" for c in analysis["top_up"][:8])
    down = ", ".join(f"{c['name']}({c['subclass']},log2FC {c['log2fc']})" for c in analysis["top_down"][:8])
    lines.append(f"Top up: {up or 'none'}")
    lines.append(f"Top down: {down or 'none'}")
    classes = sorted(analysis["class_summary"], key=lambda r: abs(r["mean_log2fc"]), reverse=True)[:10]
    lines.append("Class changes (subclass: mean log2FC, n, n_sig): " + "; ".join(
        f"{r['subclass']}: {round(r['mean_log2fc'],2)}, {r['n']}, {int(r['n_sig'])}" for r in classes))
    rf = analysis.get("random_forest") or {}
    if rf:
        lines.append(f"Random forest CV accuracy {rf.get('cv_accuracy')}, AUC {rf.get('cv_auc')}.")
    pathway_lines = lipid_pathways.pathway_digest(analysis.get("pathways"))
    if pathway_lines:
        lines.append(pathway_lines)
    if annotation_overview:
        lines.append(f"Total annotations {annotation_overview.get('total_annotations')}, "
                     f"grades {annotation_overview.get('grade_counts')}.")
    return "\n".join(lines)


def _search_query(involved_classes, involved_lipids, direction):
    terms = list(dict.fromkeys([c for c in involved_classes if c]))[:3]
    base = " ".join(terms) if terms else " ".join(involved_lipids[:2])
    return f"{base} lipidomics dysregulation biomarker".strip()


def generate_conclusions_llm(analysis, annotation_overview, llm, language, max_conclusions=5,
                             sample_context=None):
    lang_word = "Chinese" if language == "zh" else "English"
    system = (
        "You are a lipidomics domain expert. From the provided differential-analysis evidence, the LIPID "
        "MAPS pathway enrichment, and the study context, propose concise, testable biological conclusions. "
        "Interpret the changes in terms of LIPID MAPS categories and canonical lipid metabolic pathways "
        "(e.g. fatty-acid de novo synthesis/beta-oxidation, the Kennedy pathway, the Lands remodeling "
        "cycle, cardiolipin synthesis, sphingolipid de novo synthesis, sterol-ester metabolism). Base "
        f"every claim ONLY on the evidence. Write the 'statement' and 'rationale' fields in {lang_word}."
    )
    user = (
        f"{_digest(analysis, annotation_overview, sample_context)}\n\n"
        f"Propose up to {max_conclusions} distinct biological conclusions. Prefer conclusions that connect "
        "the dysregulated lipid classes to a specific LIPID MAPS metabolic pathway when the evidence supports it. "
        "Return a JSON object of the form {\"conclusions\": [ ... ]}, where each element has: "
        "id (int), statement (string), direction ('up'|'down'|'mixed'), "
        "involved_classes (array of lipid subclass strings), "
        "involved_lipids (array of specific lipid names from the evidence), "
        "involved_pathways (array of LIPID MAPS pathway names from the evidence, may be empty), "
        "rationale (string), search_query (short English literature query string)."
    )
    data = llm.chat_json(system, user)
    if isinstance(data, dict):
        data = data.get("conclusions") or data.get("results") or [data]
    conclusions = []
    for i, item in enumerate(data[:max_conclusions], start=1):
        conclusions.append({
            "id": int(item.get("id", i)),
            "statement": str(item.get("statement", "")).strip(),
            "direction": str(item.get("direction", "mixed")),
            "involved_classes": list(item.get("involved_classes", []) or []),
            "involved_lipids": list(item.get("involved_lipids", []) or []),
            "involved_pathways": list(item.get("involved_pathways", []) or []),
            "rationale": str(item.get("rationale", "")).strip(),
            "search_query": str(item.get("search_query", "")).strip()
                            or _search_query(item.get("involved_classes", []), item.get("involved_lipids", []),
                                             item.get("direction", "mixed")),
            "source": "llm",
        })
    return [c for c in conclusions if c["statement"]]


def generate_conclusions_rule(analysis, annotation_overview, language, max_conclusions=5):
    zh = language == "zh"
    conclusions = []
    class_summary = [r for r in analysis["class_summary"] if r["n_sig"] >= 1]
    class_summary = sorted(class_summary, key=lambda r: abs(r["mean_log2fc"]), reverse=True)

    diff = analysis.get("differential_table")
    for rank, row in enumerate(class_summary[:max_conclusions - 1], start=1):
        subclass = row["subclass"]
        direction = "up" if row["mean_log2fc"] > 0 else "down"
        dir_word = ("升高" if direction == "up" else "降低") if zh else ("increased" if direction == "up" else "decreased")
        lipids = []
        if diff is not None:
            sub = diff[(diff["subclass"].astype(str) == subclass) & (diff["significant"])]
            sub = sub.reindex(sub["neg_log10_fdr"].sort_values(ascending=False).index).head(5)
            lipids = [str(x) for x in sub.get("sub_chain_name", sub.get("total_chain_name", [])).tolist()]
        if zh:
            statement = f"在 {analysis['case']} 组相比 {analysis['control']} 组，{subclass} 类脂质整体{dir_word}"
            rationale = f"{subclass} 类平均 log2FC 为 {round(row['mean_log2fc'],2)}，其中 {int(row['n_sig'])} 个特征显著。"
        else:
            statement = f"{subclass} lipids are overall {dir_word} in {analysis['case']} vs {analysis['control']}"
            rationale = f"{subclass} mean log2FC {round(row['mean_log2fc'],2)} with {int(row['n_sig'])} significant features."
        conclusions.append({
            "id": rank,
            "statement": statement,
            "direction": direction,
            "involved_classes": [subclass],
            "involved_lipids": lipids,
            "involved_pathways": [p["name_en"] for p in (analysis.get("pathways") or {}).get("enrichment", [])
                                  if subclass in (p.get("subclasses") or "")][:2],
            "rationale": rationale,
            "search_query": _search_query([subclass], lipids, direction),
            "source": "rule",
        })

    # LIPID MAPS pathway enrichment conclusion / LIPID MAPS 代谢通路富集结论
    pathways = (analysis.get("pathways") or {}).get("enrichment", [])
    enriched = [p for p in pathways if p.get("n_sig", 0) >= 1]
    if enriched and len(conclusions) < max_conclusions:
        top = enriched[0]
        pname = top["name_zh"] if zh else top["name_en"]
        p_dir = "up" if (top.get("mean_log2fc") or 0) > 0 else "down"
        dir_word = ("升高" if p_dir == "up" else "降低") if zh else ("increased" if p_dir == "up" else "decreased")
        if zh:
            statement = (f"在 {analysis['case']} 组，{pname} 通路相关脂质整体{dir_word}"
                         f"（{int(top['n_sig'])}/{int(top['n_features'])} 个特征显著）")
            rationale = (f"该通路涉及亚类 {top['subclasses']}，平均 log2FC {top['mean_log2fc']}；"
                         f"超几何富集 p={('%.1e' % top['enrich_p']) if top.get('enrich_p')==top.get('enrich_p') else 'NA'}。"
                         f"通路说明：{top['desc_zh']}")
        else:
            statement = (f"{pname} pathway lipids are overall {dir_word} in {analysis['case']} vs "
                         f"{analysis['control']} ({int(top['n_sig'])}/{int(top['n_features'])} significant)")
            rationale = (f"Subclasses {top['subclasses']}, mean log2FC {top['mean_log2fc']}; hypergeometric "
                         f"enrichment p={('%.1e' % top['enrich_p']) if top.get('enrich_p')==top.get('enrich_p') else 'NA'}. "
                         f"Pathway: {top['desc_en']}")
        conclusions.append({
            "id": len(conclusions) + 1,
            "statement": statement, "direction": p_dir,
            "involved_classes": [s.strip() for s in top["subclasses"].split(",") if s.strip()],
            "involved_lipids": [],
            "involved_pathways": [top["name_en"]],
            "rationale": rationale,
            "search_query": f"{top['name_en']} lipidomics dysregulation",
            "source": "rule",
        })

    # Unsaturation trend conclusion / 不饱和度趋势结论
    if diff is not None:
        sub = diff.dropna(subset=["total_unsaturation", "log2fc"])
        if len(sub) >= 8 and sub["total_unsaturation"].nunique() >= 3:
            corr = float(np.corrcoef(sub["total_unsaturation"], sub["log2fc"])[0, 1])
            if abs(corr) >= 0.2:
                trend = ("更不饱和" if corr > 0 else "更饱和") if zh else ("more unsaturated" if corr > 0 else "more saturated")
                if zh:
                    statement = f"{analysis['case']} 组倾向于富集{trend}的脂质（不饱和度与 log2FC 相关 r={round(corr,2)}）"
                    rationale = f"不饱和度与 log2FC 的 Pearson 相关为 {round(corr,2)}。"
                else:
                    statement = f"{analysis['case']} tends to enrich {trend} lipids (unsaturation vs log2FC r={round(corr,2)})"
                    rationale = f"Pearson r between double-bond number and log2FC is {round(corr,2)}."
                conclusions.append({
                    "id": len(conclusions) + 1,
                    "statement": statement, "direction": "mixed",
                    "involved_classes": [], "involved_lipids": [], "involved_pathways": [],
                    "rationale": rationale,
                    "search_query": "lipid unsaturation remodeling desaturase lipidomics",
                    "source": "rule",
                })

    # Random-forest separability conclusion / 随机森林可分性结论
    rf = analysis.get("random_forest") or {}
    if rf and rf.get("cv_auc") == rf.get("cv_auc"):  # not NaN
        top_feats = ", ".join(f[0] for f in rf.get("selected_features", [])[:5])
        if zh:
            statement = (f"基于 Top 脂质特征的随机森林可区分 {analysis['case']} 与 {analysis['control']}"
                         f"（交叉验证 AUC={round(rf['cv_auc'],2)}，准确率={round(rf['cv_accuracy'],2)}）")
            rationale = f"关键特征包括 {top_feats}。"
        else:
            statement = (f"A random forest on top lipid features separates {analysis['case']} from "
                         f"{analysis['control']} (CV AUC={round(rf['cv_auc'],2)}, ACC={round(rf['cv_accuracy'],2)})")
            rationale = f"Key features include {top_feats}."
        conclusions.append({
            "id": len(conclusions) + 1,
            "statement": statement, "direction": "mixed",
            "involved_classes": [], "involved_lipids": [f[0] for f in rf.get("selected_features", [])[:5]],
            "involved_pathways": [],
            "rationale": rationale,
            "search_query": "lipidomics random forest classifier biomarker panel",
            "source": "rule",
        })

    return conclusions[:max_conclusions]


def generate_conclusions(analysis, annotation_overview, llm, language, max_conclusions=5, log=print,
                         sample_context=None):
    """Try the LLM engine first; fall back to the rule engine on any failure or when disabled.

    优先用大模型引擎；不可用或失败时回退到规则引擎。
    """
    if llm is not None and llm.is_enabled():
        try:
            conclusions = generate_conclusions_llm(analysis, annotation_overview, llm, language,
                                                   max_conclusions, sample_context=sample_context)
            if conclusions:
                return conclusions
            log("LLM returned no conclusions; using rule engine / 大模型未返回结论，改用规则引擎")
        except Exception as e:
            log(f"LLM conclusion generation failed ({type(e).__name__}: {e}); using rule engine "
                f"/ 大模型生成结论失败，改用规则引擎")
    return generate_conclusions_rule(analysis, annotation_overview, language, max_conclusions)
