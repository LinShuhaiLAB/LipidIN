"""One verification sub-agent per biological conclusion (n conclusions -> n sub-agents).

每条生物学结论一个验证子智能体（n 条结论 -> n 个子智能体）。

Each sub-agent evaluates its conclusion along four axes:
  1. Annotation confidence  - evidence type, fragment completeness, cosine / IoU of the involved lipids
  2. RT confidence          - RT-fit score and on-line fraction of the involved lipids
  3. Biological support     - LLM (or rule) judgement of plausibility given the evidence
  4. Similar reports        - online literature retrieval + LLM (or rule) similarity judgement
The sub-agents run concurrently; each returns a structured report and an overall verdict.

每个子智能体从四个维度评估其结论：注释可信度、RT 可信度、生物结论支撑度、是否有相似报道。子智能体
并发运行，各自返回结构化报告与总体判定。
"""

import concurrent.futures

import numpy as np
import pandas as pd

from literature_search import find_similar_reports, format_reports_for_prompt
from annotation_stats import grade_row


RT_ON_LINE = "在拟合线上"


def _match_rows(conclusion, diff_table):
    if diff_table is None or len(diff_table) == 0:
        return diff_table
    lipids = [str(x).strip() for x in conclusion.get("involved_lipids", []) if str(x).strip()]
    classes = [str(x).strip() for x in conclusion.get("involved_classes", []) if str(x).strip()]

    name_cols = [c for c in ("sub_chain_name", "total_chain_name") if c in diff_table.columns]
    mask = pd.Series(False, index=diff_table.index)
    if lipids and name_cols:
        for col in name_cols:
            mask |= diff_table[col].astype(str).str.strip().isin(lipids)
    if classes and "subclass" in diff_table.columns:
        mask |= diff_table["subclass"].astype(str).str.strip().isin(classes)
    subset = diff_table[mask]
    if subset.empty and "significant" in diff_table.columns:
        subset = diff_table[diff_table["significant"]]
    return subset


def _safe_mean(series):
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.mean()) if not values.empty else None


def assess_annotation(subset):
    if subset is None or len(subset) == 0:
        return {"n": 0, "score": 0.0, "level": "low", "detail": "no matching annotations / 无匹配注释"}
    ev = subset.get("evidence_type", pd.Series(["" for _ in range(len(subset))]))
    frag = subset.get("fragment_check", pd.Series(["" for _ in range(len(subset))]))
    grades = [grade_row(e, f) for e, f in zip(ev, frag)]
    ms2_fraction = float(np.mean([str(e) == "MS2注释" for e in ev]))
    mean_cos = _safe_mean(subset.get("cosine_similarity"))
    mean_iou = _safe_mean(subset.get("intersection_over_union"))

    score = 0.5 * ms2_fraction + 0.3 * (mean_cos or 0.0) + 0.2 * (mean_iou or 0.0)
    level = "high" if score >= 0.7 else ("medium" if score >= 0.4 else "low")
    grade_counts = {g: grades.count(g) for g in set(grades)}
    return {
        "n": int(len(subset)),
        "score": round(score, 3),
        "level": level,
        "ms2_fraction": round(ms2_fraction, 3),
        "mean_cosine": None if mean_cos is None else round(mean_cos, 3),
        "mean_iou": None if mean_iou is None else round(mean_iou, 3),
        "grade_counts": grade_counts,
        "detail": f"{len(subset)} lipids, MS2 fraction {ms2_fraction:.2f}, "
                  f"mean cosine {('%.3f' % mean_cos) if mean_cos is not None else 'NA'}, "
                  f"mean IoU {('%.3f' % mean_iou) if mean_iou is not None else 'NA'}",
    }


def assess_rt(subset):
    if subset is None or len(subset) == 0:
        return {"n": 0, "score": 0.0, "level": "low", "detail": "no matching lipids / 无匹配脂质"}
    mean_rt_fit = _safe_mean(subset.get("rt_fit_score"))
    status = subset.get("rt_fit_status", pd.Series(["" for _ in range(len(subset))])).astype(str)
    has_status = bool((status.str.strip() != "").any())
    on_line_fraction = float(np.mean([s == RT_ON_LINE for s in status])) if has_status else None
    score = mean_rt_fit if mean_rt_fit is not None else (on_line_fraction or 0.0)
    level = "high" if score >= 0.7 else ("medium" if score >= 0.4 else "low")
    detail = f"mean RT-fit score {('%.3f' % mean_rt_fit) if mean_rt_fit is not None else 'NA'}"
    if on_line_fraction is not None:
        detail += f", on-line fraction {on_line_fraction:.2f}"
    return {
        "n": int(len(subset)),
        "score": round(float(score), 3),
        "level": level,
        "mean_rt_fit_score": None if mean_rt_fit is None else round(mean_rt_fit, 3),
        "on_line_fraction": None if on_line_fraction is None else round(on_line_fraction, 3),
        "detail": detail,
    }


def assess_biology(conclusion, annotation, rt, llm, language, sample_context=None):
    if llm is not None and llm.is_enabled():
        lang_word = "Chinese" if language == "zh" else "English"
        system = ("You are a rigorous lipidomics reviewer. Judge how well the stated conclusion is "
                  "supported by the analytical evidence and is plausible in the study context. "
                  f"Be skeptical and concise. Write 'comment' in {lang_word}.")
        user = (
            f"Study context: {sample_context or 'not provided'}\n"
            f"Conclusion: {conclusion['statement']}\n"
            f"Rationale: {conclusion.get('rationale','')}\n"
            f"Annotation evidence: {annotation}\n"
            f"RT evidence: {rt}\n\n"
            "Return JSON: {\"support_score\": 0..1, \"comment\": string}."
        )
        try:
            data = llm.chat_json(system, user)
            return {"score": float(data.get("support_score", 0.5)),
                    "comment": str(data.get("comment", "")).strip(), "source": "llm"}
        except Exception:
            pass
    # Rule fallback: support tracks analytical confidence.
    # 规则回退：支撑度跟随分析置信度。
    score = round(0.6 * annotation["score"] + 0.4 * rt["score"], 3)
    if language == "zh":
        comment = f"依据注释可信度({annotation['level']})与 RT 可信度({rt['level']})综合评估支撑度。"
    else:
        comment = f"Support estimated from annotation confidence ({annotation['level']}) and RT confidence ({rt['level']})."
    return {"score": score, "comment": comment, "source": "rule"}


def assess_literature(conclusion, llm, language):
    reports = find_similar_reports(conclusion.get("search_query", conclusion["statement"]), limit=5)
    has_similar = len(reports) > 0
    if llm is not None and llm.is_enabled() and reports:
        lang_word = "Chinese" if language == "zh" else "English"
        system = ("You check whether retrieved literature reports findings similar to the conclusion. "
                  f"Write 'comment' in {lang_word}.")
        user = (
            f"Conclusion: {conclusion['statement']}\n\n"
            f"Retrieved references:\n{format_reports_for_prompt(reports)}\n\n"
            "Return JSON: {\"has_similar_report\": true/false, \"comment\": string}."
        )
        try:
            data = llm.chat_json(system, user)
            has_similar = bool(data.get("has_similar_report", has_similar))
            comment = str(data.get("comment", "")).strip()
        except Exception:
            comment = ""
    else:
        if language == "zh":
            comment = ("存在主题相关的公开文献。" if reports else "未检索到明显相关的文献。")
        else:
            comment = ("Topically related literature exists." if reports
                       else "No clearly related references retrieved.")
    return {
        "has_similar_report": has_similar,
        "n_references": len(reports),
        "references": reports,
        "comment": comment,
        "score": 0.8 if has_similar else 0.3,
    }


def _overall(annotation, rt, biology, literature, language):
    score = (0.30 * annotation["score"] + 0.25 * rt["score"]
             + 0.30 * biology["score"] + 0.15 * literature["score"])
    if score >= 0.7:
        verdict = "strongly supported" if language == "en" else "强支撑"
    elif score >= 0.45:
        verdict = "moderately supported" if language == "en" else "中等支撑"
    else:
        verdict = "weakly supported" if language == "en" else "弱支撑"
    return round(score, 3), verdict


def verify_one(conclusion, diff_table, llm, language, log=print, sample_context=None):
    log(f"Sub-agent verifying conclusion #{conclusion['id']} / 子智能体正在验证结论 #{conclusion['id']}")
    subset = _match_rows(conclusion, diff_table)
    annotation = assess_annotation(subset)
    rt = assess_rt(subset)
    biology = assess_biology(conclusion, annotation, rt, llm, language, sample_context=sample_context)
    literature = assess_literature(conclusion, llm, language)
    overall_score, verdict = _overall(annotation, rt, biology, literature, language)
    return {
        "conclusion": conclusion,
        "annotation_confidence": annotation,
        "rt_confidence": rt,
        "biological_support": biology,
        "similar_reports": literature,
        "overall_score": overall_score,
        "verdict": verdict,
    }


def verify_conclusions(conclusions, diff_table, llm, language, log=print, max_workers=None,
                       sample_context=None):
    """Launch one sub-agent per conclusion, concurrently, and collect their reports in order.

    每条结论并发启动一个子智能体，按顺序收集其报告。
    """
    if not conclusions:
        return []
    workers = max_workers or min(len(conclusions), 6)
    results = [None] * len(conclusions)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_index = {
            executor.submit(verify_one, c, diff_table, llm, language, log, sample_context): i
            for i, c in enumerate(conclusions)
        }
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            try:
                results[index] = future.result()
            except Exception as e:
                log(f"Verification sub-agent failed for #{conclusions[index]['id']}: {e} "
                    f"/ 结论 #{conclusions[index]['id']} 的验证子智能体失败")
                results[index] = {"conclusion": conclusions[index], "error": str(e),
                                  "overall_score": 0.0, "verdict": "error"}
    return results
