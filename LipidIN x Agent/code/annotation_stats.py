"""Annotation counting, confidence grading, and per-class/per-polarity summaries.

注释计数、置信度分级，以及按分类/极性的汇总。

Reads the master_table.csv (or master_table_quant.csv) of one polarity and derives: total
annotations, evidence-type breakdown, a confidence grade per annotation, score distributions and
per-subclass counts. These feed both the QC/report sections and the conclusion verification.

读取单极性的 master_table.csv（或 master_table_quant.csv），得到：总注释数、证据类型分布、每条注释的
置信度分级、打分分布与各 subclass 计数。既用于质控/报告，也用于结论验证。
"""

import numpy as np
import pandas as pd


EVIDENCE_MS2 = "MS2注释"
MISSING_TAG = "缺失"

# Confidence grades / 置信度分级
GRADE_1 = "Level 1 (MS2, complete fragments)"
GRADE_2 = "Level 2 (MS2, partial fragments)"
GRADE_3 = "Level 3 (MS1 inference + RT)"
GRADE_UNKNOWN = "Ungraded"


def grade_row(evidence_type, fragment_check):
    evidence = str(evidence_type)
    fragment = str(fragment_check)
    if evidence == EVIDENCE_MS2:
        if MISSING_TAG in fragment:
            return GRADE_2
        return GRADE_1
    if evidence.startswith("MS1"):
        return GRADE_3
    return GRADE_UNKNOWN


def annotation_confidence_grade(df):
    """Add a 'confidence_grade' column derived from evidence type + fragment completeness.

    根据证据类型与碎片完整度，新增 'confidence_grade' 列。
    """
    df = df.copy()
    ev = df.get("evidence_type", pd.Series([""] * len(df)))
    frag = df.get("fragment_check", pd.Series([""] * len(df)))
    df["confidence_grade"] = [grade_row(e, f) for e, f in zip(ev, frag)]
    return df


def _score_stats(series):
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return {"n": 0, "median": None, "mean": None}
    return {"n": int(values.size), "median": round(float(values.median()), 4),
            "mean": round(float(values.mean()), 4)}


def summarize_annotations(master_csv, polarity):
    """Compute the annotation statistics dict for one polarity's master table.

    计算单极性注释总表的统计字典。
    """
    df = pd.read_csv(master_csv, encoding="utf-8-sig")
    df = annotation_confidence_grade(df)

    evidence_counts = df.get("evidence_type", pd.Series(dtype=str)).value_counts().to_dict()
    grade_counts = df["confidence_grade"].value_counts().to_dict()
    subclass_counts = df.get("subclass", pd.Series(dtype=str)).astype(str).value_counts().head(25).to_dict()

    return {
        "polarity": polarity,
        "total_annotations": int(len(df)),
        "evidence_counts": {str(k): int(v) for k, v in evidence_counts.items()},
        "grade_counts": {str(k): int(v) for k, v in grade_counts.items()},
        "n_subclasses": int(df.get("subclass", pd.Series(dtype=str)).astype(str).nunique()),
        "top_subclass_counts": {str(k): int(v) for k, v in subclass_counts.items()},
        "cosine_similarity": _score_stats(df.get("cosine_similarity")),
        "intersection_over_union": _score_stats(df.get("intersection_over_union")),
        "rt_fit_score": _score_stats(df.get("rt_fit_score")),
        "master_csv": master_csv,
    }


def combine_polarity_stats(stats_list):
    """Merge per-polarity annotation stats into a combined overview.

    把各极性的注释统计合并成总体概览。
    """
    stats_list = [s for s in stats_list if s]
    total = sum(s["total_annotations"] for s in stats_list)
    combined_grades = {}
    combined_evidence = {}
    for s in stats_list:
        for k, v in s["grade_counts"].items():
            combined_grades[k] = combined_grades.get(k, 0) + v
        for k, v in s["evidence_counts"].items():
            combined_evidence[k] = combined_evidence.get(k, 0) + v
    return {
        "total_annotations": total,
        "per_polarity": {s["polarity"]: s["total_annotations"] for s in stats_list},
        "grade_counts": combined_grades,
        "evidence_counts": combined_evidence,
        "polarity_details": stats_list,
    }
