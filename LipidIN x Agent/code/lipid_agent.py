"""Interactive LipidIN agent: runs the whole pipeline and downstream analysis, then writes a report.

交互式 LipidIN 智能体：运行全流程与下游分析，并生成报告。

Seven steps:
  1. Ask the session language (Chinese or English) and use it for everything afterwards.
  2. Ask for parameters, detect pos/neg data folders and confirm what to run.
  3. Run the qualitative + quantitative pipeline (run.py steps).
  4. Ask for sample information, then run quality control.
  5. Auto-run lipid differential analysis (volcano, heatmap, FA-chain change, class bubble, random forest).
  6. Derive biological conclusions and verify each with one sub-agent (annotation / RT / support / literature).
  7. Write the full Markdown report and print a short summary.

七个步骤：选语言 -> 参数与文件夹检测 -> 定性定量流程 -> 样本信息与质控 -> 脂质差异分析 ->
生物学结论与多智能体验证 -> 生成完整 Markdown 报告并给出简要总结。

The interactive entry is `interactive_main()`. The reusable engine is `run_agent(config)`, which takes a
fully-resolved config dict and can be called programmatically (e.g. for testing) without any prompts.
交互入口是 interactive_main()；可复用引擎是 run_agent(config)，接收已解析好的配置，可无提示地程序化调用。
"""

import os
import sys
from pathlib import Path


def configure_utf8_stdio():
    for stream_name in ("stdout", "stderr", "stdin"):
        stream = getattr(sys, stream_name, None)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


configure_utf8_stdio()

import pandas as pd

from i18n import Messages
from llm_client import LLMClient
import pipeline_runner
from qc_analysis import run_qc, group_of, is_qc_group
from annotation_stats import summarize_annotations, combine_polarity_stats, annotation_confidence_grade
from master_table import EVIDENCE_MS2, EVIDENCE_MS1_ONLY, EVIDENCE_MS1_UNANNOTATED
import lipid_analysis
from bio_conclusions import generate_conclusions
from conclusion_verification import verify_conclusions
import report_writer


CODE_DIR = Path(__file__).resolve().parent
DEFAULT_LIBRARY_DIR = str(CODE_DIR.parent / "library")


# ------------------------------------------------------------------ #
# Group / comparison helpers
# 分组 / 比较辅助
# ------------------------------------------------------------------ #
def infer_groups(sample_names):
    return {s: group_of(s) for s in sample_names}


def pick_case_control(groups):
    """Pick the two most populous non-QC groups as (case, control).

    选择样本数最多的两个非 QC 分组作为 (病例, 对照)。
    """
    counts = {}
    for sample, group in groups.items():
        if is_qc_group(group):
            continue
        counts[group] = counts.get(group, 0) + 1
    ordered = sorted(counts, key=lambda g: (-counts[g], g))
    if len(ordered) < 2:
        return None, None
    # Convention: the alphabetically later group is treated as case (e.g. B vs A).
    # 约定：字母靠后的组作为病例组（例如 B vs A）。
    top_two = sorted(ordered[:2])
    return top_two[1], top_two[0]


def build_combined_corrected_table(qc_results):
    """Concatenate each polarity's total-area-corrected table into one lipidome (tagging polarity).

    把各极性的总峰面积校正表拼成一个总脂质组（标注极性）。
    """
    frames = []
    for polarity, qc in qc_results.items():
        corrected_csv = qc.get("corrected_csv")
        if corrected_csv and os.path.exists(corrected_csv):
            df = pd.read_csv(corrected_csv, encoding="utf-8-sig")
            df.insert(0, "polarity", polarity)
            frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


from annotation_stats import GRADE_1, GRADE_2, GRADE_3

GRADE_BY_DIGIT = {"1": GRADE_1, "2": GRADE_2, "3": GRADE_3}

DEFAULT_ANNOTATION_FILTER = {
    "include_ms2": True,
    "include_ms1_only": True,
    "include_unannotated": True,
    "allowed_grades": None,   # None = keep all grades; else a set of grade strings
}


def apply_annotation_filter(df, filter_cfg, msgs, log):
    """Keep only the annotation rows the user allowed (by evidence type and confidence grade).

    仅保留用户允许参与的注释行（按证据类型与置信度等级）。
    """
    if not filter_cfg:
        return df
    df = annotation_confidence_grade(df)
    n0 = len(df)
    evidence = df.get("evidence_type", pd.Series([""] * len(df))).astype(str)
    keep = pd.Series(True, index=df.index)
    if not filter_cfg.get("include_ms1_only", True):
        keep &= evidence != EVIDENCE_MS1_ONLY
    if not filter_cfg.get("include_unannotated", True):
        keep &= evidence != EVIDENCE_MS1_UNANNOTATED
    if not filter_cfg.get("include_ms2", True):
        keep &= evidence != EVIDENCE_MS2
    allowed = filter_cfg.get("allowed_grades")
    if allowed:
        keep &= df["confidence_grade"].isin(set(allowed))
    out = df[keep].copy()
    log(msgs.t("filter_applied", stage="annotation", n0=n0, n1=len(out), dropped=n0 - len(out)))
    return out


def apply_class_rt_rules(df, class_rules, log=print):
    """Apply per-lipid-class RT-strictness rules, e.g. {"FA": {"rt_min": 0.85}, "*": {"rt_min": 0.3}}.

    应用按脂质类别的 RT 严格度规则，例如对 FA 类要求 rt_fit_score>=0.85，其它类别更宽松。

    Keys are matched against the 'subclass' column (case-insensitive); "*"/"other"/"default" is the
    fallback for any class not named. A row is dropped only when its class has rt_min>0 and its
    rt_fit_score is missing or below rt_min. Classes with rt_min==0 keep every row.
    键与 subclass 列匹配（不区分大小写）；"*"/"other"/"default" 为未命名类别的兜底规则。仅当某类
    rt_min>0 且该行 rt_fit_score 缺失或低于 rt_min 时才剔除；rt_min==0 的类别保留全部行。
    """
    if not class_rules or "subclass" not in df.columns:
        return df
    subclass = df["subclass"].astype(str)
    rt_score = pd.to_numeric(df.get("rt_fit_score"), errors="coerce")

    def _fallback(rules):
        for key in ("*", "other", "others", "default", "其它", "其他"):
            if key in rules:
                return float(rules[key].get("rt_min", 0.0))
        return 0.0

    rt_min = pd.Series(_fallback(class_rules), index=df.index)
    for cls, rule in class_rules.items():
        if str(cls).lower() in ("*", "other", "others", "default", "其它", "其他"):
            continue
        mask = subclass.str.upper() == str(cls).upper()
        rt_min[mask] = float(rule.get("rt_min", 0.0))

    strict = rt_min > 0
    keep = (~strict) | (rt_score.notna() & (rt_score >= rt_min))
    n0 = len(df)
    out = df[keep].copy()
    if n0 != len(out):
        log(f"Per-class RT rules: {n0} -> {len(out)} features ({n0 - len(out)} removed)"
            f" / 按类别 RT 规则：{n0} -> {len(out)} 个特征（移除 {n0 - len(out)} 个）")
    return out


def apply_cv_filter(df, sample_names, groups, cv_cfg, msgs, log):
    """Drop features whose CV across the QC samples exceeds the user threshold.

    剔除 QC 样本间 CV 超过用户阈值的特征。
    """
    if not cv_cfg or not cv_cfg.get("enabled"):
        return df
    qc_samples = [s for s in sample_names if is_qc_group(groups.get(s)) and s in df.columns]
    if len(qc_samples) < 2:
        return df
    max_cv = float(cv_cfg.get("max_cv_percent", 30.0))
    values = df[qc_samples].apply(pd.to_numeric, errors="coerce")
    mean = values.mean(axis=1)
    std = values.std(axis=1, ddof=1)
    cv = 100.0 * std / mean.replace(0, pd.NA)
    n0 = len(df)
    # Drop only where CV is computable and too high; keep QC-undetected features.
    # 仅在 CV 可计算且过高时剔除；QC 未检出的特征保留。
    drop = (mean > 0) & (cv > max_cv)
    out = df[~drop.fillna(False)].copy()
    log(msgs.t("filter_applied", stage=f"CV>{max_cv:.0f}%", n0=n0, n1=len(out), dropped=n0 - len(out)))
    return out


# ------------------------------------------------------------------ #
# The reusable engine (no prompts)
# 可复用引擎（无交互提示）
# ------------------------------------------------------------------ #
def run_agent(config):
    configure_utf8_stdio()
    language = config.get("language", "en")
    msgs = Messages(language)

    def say(key, **fmt):
        print(msgs.t(key, **fmt))

    def log(message):
        print(message)

    llm = LLMClient.from_dict(config.get("llm"))
    data_root = config["data_root"]
    library_dir = config.get("library_dir", DEFAULT_LIBRARY_DIR)
    polarities = config["polarities"]
    params = dict(pipeline_runner.DEFAULT_PARAMS)
    params.update(config.get("params", {}) or {})
    output_dir = config.get("output_dir", os.path.join(data_root, "agent_report"))
    os.makedirs(output_dir, exist_ok=True)
    max_conclusions = config.get("max_conclusions", 5)

    if not llm.is_enabled():
        say("llm_none")

    # ---- Step 3 | pipeline ---- #
    say("banner_pipeline")
    polarity_paths = {}
    sample_names = None
    for polarity in polarities:
        say("pipeline_polarity_start", polarity=polarity)
        if config.get("run_pipeline", True):
            paths = pipeline_runner.run_pipeline_for_polarity(
                data_root, polarity, library_dir, params=params, log=log)
        else:
            paths = pipeline_runner.polarity_paths(data_root, polarity, library_dir)
            _, mzml_files = pipeline_runner.validate_polarity_inputs(data_root, polarity, library_dir)
            paths["sample_names"] = [os.path.splitext(os.path.basename(p))[0] for p in mzml_files]
        polarity_paths[polarity] = paths
        sample_names = sample_names or paths["sample_names"]
    say("pipeline_done", polarities=", ".join(polarities))

    # ---- Step 4 | sample info + QC ---- #
    say("banner_sampleinfo")
    groups = config.get("groups") or infer_groups(sample_names)
    say("detected_samples", n=len(sample_names), samples=", ".join(sample_names))
    say("auto_groups", groups=groups)

    case = config.get("case")
    control = config.get("control")
    if not case or not control:
        case, control = pick_case_control(groups)

    say("qc_running")
    qc_results = {}
    for polarity, paths in polarity_paths.items():
        qc_dir = os.path.join(data_root, f"QC_{polarity}")
        qc_results[polarity] = run_qc(
            quant_csv=paths["quant_table_csv"], output_dir=qc_dir,
            sample_names=sample_names, polarity=polarity)

    # ---- Annotation statistics ---- #
    per_polarity_stats = []
    for polarity, paths in polarity_paths.items():
        if os.path.exists(paths["master_table_csv"]):
            per_polarity_stats.append(summarize_annotations(paths["master_table_csv"], polarity))
    annotation_overview = combine_polarity_stats(per_polarity_stats)

    # ---- Annotation & CV filtering (applied before differential analysis) ---- #
    sample_context = config.get("sample_context")
    combined = build_combined_corrected_table(qc_results)
    if combined is not None:
        say("banner_filter")
        combined = apply_annotation_filter(combined, config.get("annotation_filter"), msgs, log)
        combined = apply_cv_filter(combined, sample_names, groups, config.get("cv_filter"), msgs, log)
        combined = apply_class_rt_rules(combined, config.get("class_rules"), log)

    # ---- Step 5 | differential analysis ---- #
    say("banner_analysis")
    say("analysis_running")
    analysis_dir = os.path.join(output_dir, "analysis")
    analyses = []
    if combined is not None and len(combined) and case and control:
        analysis = lipid_analysis.run_analysis(
            combined, sample_names, groups, case, control, analysis_dir,
            tag="_".join(polarities))
        analyses.append(analysis)

    # ---- Step 6 | conclusions + verification ---- #
    say("banner_conclusions")
    conclusions = []
    verifications = []
    if analyses:
        confirm_callback = config.get("confirm_callback")
        while True:
            conclusions = generate_conclusions(
                analyses[0], annotation_overview, llm, language,
                max_conclusions=max_conclusions, log=log, sample_context=sample_context)
            say("conclusions_generated", n=len(conclusions))
            if confirm_callback is None:
                break
            decision = confirm_callback(conclusions, msgs)
            if decision == "REGENERATE":
                continue
            conclusions = decision
            break
        if conclusions:
            verifications = verify_conclusions(
                conclusions, analyses[0].get("differential_table"), llm, language, log=log,
                sample_context=sample_context)
        else:
            say("no_conclusions_kept")

    # ---- Step 7 | report ---- #
    say("banner_report")
    md_path = os.path.join(output_dir, "LipidIN_report.md")
    report_params = {
        "language": language, "polarities": ",".join(polarities),
        "case": case, "control": control, "llm_provider": llm.provider,
        "annotation_filter": config.get("annotation_filter"),
        "cv_filter": config.get("cv_filter"),
        **{k: params[k] for k in ("MS2_filter", "ppm_err1", "da_err1", "ppm_err2", "da_err2")},
    }
    report_writer.build_report(
        md_path=md_path, language=language, params=report_params,
        annotation_overview=annotation_overview, qc_results=qc_results,
        analyses=analyses, conclusions=conclusions, verifications=verifications,
        extra_paths={"output_dir": output_dir}, sample_context=sample_context)
    say("report_saved", path=md_path)

    summary = report_writer.build_short_summary(
        language, annotation_overview, qc_results, analyses, verifications)
    print("\n" + msgs.t("summary_header"))
    print(summary)

    return {
        "md_path": md_path,
        "annotation_overview": annotation_overview,
        "qc_results": qc_results,
        "analyses": analyses,
        "conclusions": conclusions,
        "verifications": verifications,
        "summary": summary,
    }


# ------------------------------------------------------------------ #
# Interactive front-end
# 交互式前端
# ------------------------------------------------------------------ #
def _ask(prompt):
    try:
        return input(prompt).strip()
    except EOFError:
        return ""


def _detect_polarities(data_root):
    found = {}
    for polarity in ("pos", "neg"):
        folder = os.path.join(data_root, polarity)
        files = pipeline_runner.list_mzml_files(folder) if os.path.isdir(folder) else []
        if files:
            found[polarity] = files
    return found


def interactive_main():
    configure_utf8_stdio()

    # Step 1 | language
    lang = ""
    while lang not in ("en", "zh"):
        lang = _ask(Messages("en").t("ask_language")).lower() or "en"
    msgs = Messages(lang)
    print(msgs.t("lang_set"))

    # LLM configuration
    provider = _ask(msgs.t("ask_llm_provider")).lower() or "none"
    llm_config = {"provider": provider}
    if provider not in ("none", ""):
        llm_config["api_key"] = _ask(msgs.t("ask_llm_key")) or None
        llm_config["model"] = _ask(msgs.t("ask_llm_model")) or None
        if provider == "custom":
            llm_config["base_url"] = _ask(msgs.t("ask_llm_baseurl")) or None

    # Step 2 | parameters + folder detection
    print(msgs.t("banner_params"))
    data_root = _ask(msgs.t("ask_data_root"))
    library_dir = _ask(msgs.t("ask_library_dir")) or DEFAULT_LIBRARY_DIR

    found = _detect_polarities(data_root)
    if not found:
        print(msgs.t("detect_found_none"))
        return
    if "pos" in found and "neg" in found:
        print(msgs.t("detect_found_both", n_pos=len(found["pos"]), n_neg=len(found["neg"])))
        choice = _ask(msgs.t("ask_run_both")).lower()
        if choice == "both":
            polarities = ["pos", "neg"]
        elif choice in ("pos", "neg"):
            polarities = [choice]
        else:
            polarities = ["pos", "neg"]
    else:
        only = next(iter(found))
        print(msgs.t("detect_found_one", polarity=only, n=len(found[only])))
        if _ask(msgs.t("ask_run_one", polarity=only)).lower().startswith("y"):
            polarities = [only]
        else:
            print(msgs.t("aborted"))
            return

    params = dict(pipeline_runner.DEFAULT_PARAMS)
    print(msgs.t("params_summary", **params))
    if not _ask(msgs.t("ask_use_default_params")).lower().startswith("y"):
        for name in ("MS2_filter", "ppm_err1", "da_err1", "ppm_err2", "da_err2", "eq3_workers"):
            value = _ask(msgs.t("ask_param_value", name=name, current=params[name]))
            if value:
                try:
                    params[name] = type(params[name])(value) if params[name] is not None else float(value)
                except ValueError:
                    params[name] = value

    # Peek at samples to ask for grouping (uses the first polarity's folder).
    # 先看一眼样本，便于询问分组（使用第一个极性的目录）。
    first_files = found[polarities[0]]
    sample_names = [os.path.splitext(os.path.basename(p))[0] for p in first_files]
    print(msgs.t("banner_sampleinfo"))
    print(msgs.t("detected_samples", n=len(sample_names), samples=", ".join(sample_names)))
    groups = infer_groups(sample_names)
    print(msgs.t("auto_groups", groups=groups))
    if not _ask(msgs.t("ask_edit_groups")).lower().startswith("y"):
        for sample in sample_names:
            value = _ask(msgs.t("ask_group_for", sample=sample))
            if value:
                groups[sample] = value
    case, control = pick_case_control(groups)
    cc = _ask(msgs.t("ask_case_control"))
    if "," in cc.replace("，", ","):
        parts = [p.strip() for p in cc.replace("，", ",").split(",")]
        if len(parts) == 2 and all(parts):
            case, control = parts[0], parts[1]

    # Detailed sample / study context (disease, tissue, design ...) to guide the conclusions.
    # 详细样本/研究背景（疾病、组织、设计等），用于指导结论生成。
    sample_context = _ask(msgs.t("ask_sample_context")) or None

    # Annotation & quantification filtering, applied before the differential analysis.
    # 注释与定量筛选，在差异分析之前应用。
    print(msgs.t("banner_filter"))
    annotation_filter = dict(DEFAULT_ANNOTATION_FILTER)
    annotation_filter["include_ms1_only"] = not _ask(msgs.t("ask_include_ms1only")).lower().startswith("n")
    annotation_filter["include_unannotated"] = not _ask(msgs.t("ask_include_unannotated")).lower().startswith("n")
    grade_answer = _ask(msgs.t("ask_custom_grades")).replace("，", ",").strip()
    if grade_answer:
        chosen = [GRADE_BY_DIGIT[d.strip()] for d in grade_answer.split(",") if d.strip() in GRADE_BY_DIGIT]
        if chosen:
            annotation_filter["allowed_grades"] = chosen

    cv_filter = {"enabled": False, "max_cv_percent": 30.0}
    if _ask(msgs.t("ask_cv_filter")).lower().startswith("y"):
        cv_filter["enabled"] = True
        threshold = _ask(msgs.t("ask_cv_threshold")).strip()
        if threshold:
            try:
                cv_filter["max_cv_percent"] = float(threshold)
            except ValueError:
                pass

    def confirm_callback(conclusions, msgs):
        while True:
            print(msgs.t("confirm_conclusions_header"))
            for c in conclusions:
                print(f"  [{c['id']}] {c['statement']}")
            answer = _ask(msgs.t("ask_confirm_conclusions")).strip().lower()
            if answer in ("", "y", "yes"):
                return conclusions
            if answer == "r":
                return "REGENERATE"
            if answer == "q":
                return []
            try:
                drop_ids = {int(x) for x in answer.replace("，", ",").split(",") if x.strip()}
            except ValueError:
                continue
            kept = [c for c in conclusions if c["id"] not in drop_ids]
            print(msgs.t("conclusions_updated", n=len(kept)))
            return kept

    config = {
        "language": lang,
        "llm": llm_config,
        "data_root": data_root,
        "library_dir": library_dir,
        "polarities": polarities,
        "params": params,
        "groups": groups,
        "case": case,
        "control": control,
        "sample_context": sample_context,
        "annotation_filter": annotation_filter,
        "cv_filter": cv_filter,
        "confirm_callback": confirm_callback,
        "output_dir": os.path.join(data_root, "agent_report"),
        "run_pipeline": True,
        "max_conclusions": 5,
    }
    run_agent(config)


if __name__ == "__main__":
    interactive_main()
