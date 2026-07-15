"""Local web UI (chat + progress bar) for the LipidIN agent.

LipidIN 智能体的本地网页界面（对话框 + 进度条）。

A minimalist single-page app served by Flask. The user talks to the agent in a chat box; progress
streams live over Server-Sent Events. Every parameter can be changed at any point in the conversation
(annotation grade filter, CV cutoff, case/control, per-class RT strictness, etc.) and a downstream
stage can be re-run with the new settings without repeating the heavy pipeline. Free-form personalized
requests (e.g. "keep RT strict for FA but loose for other classes") are parsed by the configured LLM;
common structured commands work without any LLM.

由 Flask 提供的极简单页应用。用户在对话框里与智能体交流，进度通过 SSE 实时推送。对话过程中任何参数
都可随时修改（注释等级筛选、CV 阈值、病例/对照、按类别的 RT 严格度等），改完可只重跑下游步骤而无需
重复耗时的流程。自然语言的 个性化需求（如"FA 类 RT 严格、其它类别宽松"）由所配置的大模型解析；常见的
结构化命令无需大模型即可使用。

Run:  python agent_server.py           # then open http://127.0.0.1:7860
运行： python agent_server.py           # 然后打开 http://127.0.0.1:7860
"""

import os
import re
import sys
import json
import queue
import threading
import webbrowser
from pathlib import Path


def configure_utf8_stdio():
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


configure_utf8_stdio()

import pandas as pd
from flask import Flask, request, Response, jsonify

import pipeline_runner
import raw_converter
from llm_client import LLMClient
from qc_analysis import run_qc, is_qc_group, group_of
from annotation_stats import summarize_annotations, combine_polarity_stats
import lipid_analysis
from bio_conclusions import generate_conclusions
from conclusion_verification import verify_conclusions
import report_writer
from i18n import Messages
from lipid_agent import (
    infer_groups, pick_case_control, build_combined_corrected_table,
    apply_annotation_filter, apply_cv_filter, apply_class_rt_rules,
    DEFAULT_ANNOTATION_FILTER, GRADE_BY_DIGIT,
)

CODE_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = str(CODE_DIR.parent / "common_demo")
DEFAULT_LIBRARY_DIR = str(CODE_DIR.parent / "library")


def _load_logo_data_uri():
    """Read the LipidIN icon and return a base64 data URI (empty string if unavailable).

    读取 LipidIN 图标并返回 base64 data URI（读不到则返回空串）。
    """
    import base64
    icon_path = CODE_DIR / "icon.png"
    try:
        data = icon_path.read_bytes()
        return "data:image/png;base64," + base64.b64encode(data).decode("ascii")
    except Exception:
        return ""


LOGO_DATA_URI = _load_logo_data_uri()


# ------------------------------------------------------------------ #
def deep_merge(base, patch):
    for key, value in (patch or {}).items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base


class AgentUISession:
    """Holds the editable config, cached intermediate results, and the event stream.

    保存可编辑配置、缓存的中间结果，以及事件流。
    """

    def __init__(self):
        self.lock = threading.Lock()
        self.subscribers = []
        self.busy = False
        self.pending = None  # multi-turn wizard state, e.g. "llm_provider" / 多轮向导状态
        self.pending_action = None  # a run action awaiting preflight confirmation / 待确认的运行动作
        self.confirmed = False  # parameters have been confirmed for running / 参数是否已确认可运行
        self.config = {
            "language": "zh",
            "llm": {"provider": "none", "api_key": None, "model": None, "base_url": None},
            "data_root": DEFAULT_DATA_ROOT,
            "library_dir": DEFAULT_LIBRARY_DIR,
            "polarities": ["pos", "neg"],
            "params": dict(pipeline_runner.DEFAULT_PARAMS),
            "groups": None,
            "case": None,
            "control": None,
            "comparisons": None,  # list of {"case","control"} for multi-group differential / 多组对照
            "sample_context": None,
            "annotation_filter": dict(DEFAULT_ANNOTATION_FILTER),
            "cv_filter": {"enabled": False, "max_cv_percent": 30.0},
            "class_rules": None,
            "output_dir": None,
            "max_conclusions": 5,
        }
        self.cache = {}
        self.jobs = queue.Queue()
        self.worker = threading.Thread(target=self._worker_loop, daemon=True)
        self.worker.start()

    # -------- language / i18n helpers -------- #
    @property
    def lang(self):
        return self.config.get("language", "zh")

    def t(self, zh, en):
        return zh if self.lang == "zh" else en

    def msgs(self):
        return Messages(self.lang)

    def llm(self):
        return LLMClient.from_dict(self.config.get("llm"))

    # -------- event stream -------- #
    def subscribe(self):
        q = queue.Queue()
        self.subscribers.append(q)
        q.put(self._state_event())
        return q

    def unsubscribe(self, q):
        if q in self.subscribers:
            self.subscribers.remove(q)

    def push(self, event):
        for q in list(self.subscribers):
            q.put(event)

    def assistant(self, zh, en=None):
        self.push({"type": "assistant", "text": self.t(zh, en if en is not None else zh)})

    def log(self, text):
        self.push({"type": "log", "text": str(text)})
        # Nudge the progress bar when a stage prints "Step i/N".
        # 当某一步打印 "Step i/N" 时推动进度条。
        m = re.search(r"Step\s+(\d)\s*/\s*(\d)", str(text))
        if m and self.cache.get("_stage_span"):
            lo, hi = self.cache["_stage_span"]
            frac = int(m.group(1)) / max(1, int(m.group(2)))
            self.progress(lo + (hi - lo) * frac)

    def progress(self, value, label=None):
        self.push({"type": "progress", "value": max(0, min(100, round(value)))
                   , **({"label": label} if label else {})})

    def set_busy(self, busy):
        self.busy = busy
        self.push({"type": "busy", "busy": busy})

    def _state_event(self):
        cfg = self.config
        chips = {
            "data_root": cfg["data_root"],
            "polarities": ",".join(cfg["polarities"]),
            "case_control": f"{cfg.get('case') or '?'} vs {cfg.get('control') or '?'}",
            "provider": cfg["llm"].get("provider", "none"),
            "grades": _grades_label(cfg["annotation_filter"].get("allowed_grades")),
            "ms1_only": "on" if cfg["annotation_filter"].get("include_ms1_only", True) else "off",
            "cv": (f"<{cfg['cv_filter']['max_cv_percent']:.0f}%" if cfg["cv_filter"].get("enabled") else "off"),
            "class_rules": _class_rules_label(cfg.get("class_rules")),
            "context": (cfg.get("sample_context") or "")[:40],
        }
        try:
            panel = _panel_data(self)
        except Exception:
            panel = None
        return {"type": "config", "chips": chips, "busy": self.busy, "panel": panel}

    def emit_state(self):
        self.push(self._state_event())

    # -------- job queue -------- #
    def enqueue(self, action):
        self.jobs.put(action)

    def _worker_loop(self):
        while True:
            action = self.jobs.get()
            try:
                self.set_busy(True)
                STAGES[action](self)
            except Exception as e:
                import traceback
                self.log(traceback.format_exc(limit=6))
                self.assistant(f"运行 {action} 时出错：{e}", f"Error while running {action}: {e}")
            finally:
                self.set_busy(False)
                self.progress(0)
                self.emit_state()


SESSION = AgentUISession()


def _grades_label(allowed):
    if not allowed:
        return "all"
    digits = [d for d, g in GRADE_BY_DIGIT.items() if g in allowed]
    return ",".join(sorted(digits)) or "custom"


def _class_rules_label(rules):
    if not rules:
        return "none"
    return ", ".join(f"{k}:rt>={v.get('rt_min', 0)}" for k, v in rules.items())


# ------------------------------------------------------------------ #
# Multi-group comparison helpers / 多组对照辅助
# ------------------------------------------------------------------ #
def enumerate_pairs(groups):
    """All unordered pairs of non-QC groups as {case, control} (alphabetically later = case).

    所有非 QC 组的两两组合（字母靠后者作为病例）。
    """
    names = sorted({g for g in (groups or {}).values() if g and not is_qc_group(g)})
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            pairs.append({"case": b, "control": a})
    return pairs


def resolve_comparisons(s):
    """The list of {case, control} pairs to run. Explicit config wins; else the single case/control.

    需要运行的 {病例, 对照} 对列表。优先用显式配置的多组对照；否则回退到单组 case/control。
    """
    comps = s.config.get("comparisons")
    if comps:
        seen, out = set(), []
        for c in comps:
            case, control = c.get("case"), c.get("control")
            if case and control and case != control and (case, control) not in seen:
                seen.add((case, control))
                out.append({"case": case, "control": control})
        if out:
            return out
    if s.config.get("case") and s.config.get("control"):
        return [{"case": s.config["case"], "control": s.config["control"]}]
    return []


def _comparisons_label(comps):
    if not comps:
        return "auto"
    return "; ".join(f"{c['case']} vs {c['control']}" for c in comps)


def _detect_samples(s):
    """Detect per-polarity input file counts (mzML or vendor raw) and sample names.

    检测每个极性的输入文件数（mzML 或厂商原始文件）与样本名。返回 (detected, sample_names)。
    """
    detected = {}
    sample_names = None
    for polarity in s.config["polarities"]:
        # Polarity-aware: single-polarity raws are attributed to their real polarity, not both.
        # 极性感知：单极性 raw 归到其真实极性，而不是同时计入正负两个极性。
        kind, files, _src = raw_converter.list_polarity_inputs(
            s.config["data_root"], polarity, log=(lambda *a, **k: None))
        kind_label = {"mzml": "mzML", "raw": "raw", "wiff": "wiff", "none": "none"}.get(kind, kind)
        detected[polarity] = {"count": len(files), "kind": kind_label}
        if files and sample_names is None:
            if kind == "mzml":
                sample_names = [os.path.splitext(os.path.basename(f))[0] for f in files]
            else:
                sample_names = [raw_converter._base_name(f) for f in files]
    return detected, sample_names


def _runnable_polarities(s):
    """Split the configured polarities into (have-files, empty). Empty ones must be skipped, not crash.

    把配置的极性分为（有文件的、没文件的）。没文件的应跳过而不是让流程崩溃。
    """
    runnable, empty = [], []
    for pol in s.config["polarities"]:
        _k, files, _src = raw_converter.list_polarity_inputs(s.config["data_root"], pol, log=(lambda *a, **k: None))
        (runnable if files else empty).append(pol)
    return runnable, empty


def _resolve_run_polarities(s):
    """Return the polarities to actually run, warning about (and skipping) any with no input files.

    返回真正要运行的极性，并对没有输入文件的极性给出提示并跳过。全为空时返回 []。
    """
    runnable, empty = _runnable_polarities(s)
    if empty and runnable:
        s.assistant(f"极性 {','.join(empty)} 下没有数据文件，已自动跳过，仅运行 {','.join(runnable)}。",
                    f"No data for {','.join(empty)}; skipping it and running only {','.join(runnable)}.")
    return runnable


def _panel_data(s):
    """Structured data for the left UI panel: run params, groups, comparisons, biological question.

    左侧界面面板的结构化数据：运行参数、组别、对照、生物学问题。
    """
    cfg = s.config
    detected, sample_names = _detect_samples(s)
    # groups + counts (prefer already-inferred groups, else infer from detected samples)
    groups = cfg.get("groups") or (
        {sn: group_of(sn) for sn in sample_names} if sample_names else {})
    counts = {}
    for g in groups.values():
        counts[g] = counts.get(g, 0) + 1
    comps = resolve_comparisons(s)
    used = set()
    for c in comps:
        used.add(c["case"]); used.add(c["control"])
    group_list = []
    for g, n in sorted(counts.items()):
        role = "qc" if is_qc_group(g) else ("case/control" if g in used else "group")
        group_list.append({"name": g, "n": n, "role": role})
    af = cfg["annotation_filter"]
    params = cfg["params"]
    return {
        "data_root": cfg["data_root"],
        "library_dir": cfg["library_dir"],
        "polarities": cfg["polarities"],
        "inputs": detected,
        "grades": _grades_label(af.get("allowed_grades")),
        "ms1_only": "on" if af.get("include_ms1_only", True) else "off",
        "cv": (f"<{cfg['cv_filter']['max_cv_percent']:.0f}%" if cfg["cv_filter"].get("enabled") else "off"),
        "class_rules": _class_rules_label(cfg.get("class_rules")),
        "ppm": {"ms1": params.get("ppm_err1"), "ms2": params.get("ppm_err2"),
                "ms2_filter": params.get("MS2_filter")},
        "provider": cfg["llm"].get("provider", "none"),
        "model": cfg["llm"].get("model") or "",
        "groups": group_list,
        "comparisons": comps,
        "context": cfg.get("sample_context") or "",
        "busy": s.busy,
    }


# ------------------------------------------------------------------ #
# Stage runners (reuse the pipeline building blocks, with caching)
# 阶段执行器（复用流程构件，并带缓存）
# ------------------------------------------------------------------ #
def _resolve_paths(s):
    paths = {}
    for polarity in s.config["polarities"]:
        paths[polarity] = pipeline_runner.polarity_paths(
            s.config["data_root"], polarity, s.config["library_dir"])
    return paths


def stage_annotation(s):
    polarities = _resolve_run_polarities(s)
    if not polarities:
        s.assistant("所有极性下都没有数据文件，请检查数据目录。", "No data files for any polarity; check the data folder.")
        return
    s.cache["_stage_span"] = (5, 38)
    s.progress(5, "annotation")
    s.assistant("开始注释（仅定性）……", "Running annotation (qualitative only) ...")
    paths = dict(s.cache.get("paths") or {})
    sample_names = None
    for polarity in polarities:
        p = pipeline_runner.run_annotation_for_polarity(
            s.config["data_root"], polarity, s.config["library_dir"],
            params=s.config["params"], log=s.log)
        paths[polarity] = p
        sample_names = sample_names or p["sample_names"]
    s.cache["paths"] = paths
    s.cache["sample_names"] = sample_names
    s.progress(38)
    s.assistant("注释完成。可继续 `定量`。", "Annotation done. Next: `quantify`.")


def stage_quantify(s):
    polarities = _resolve_run_polarities(s)
    if not polarities:
        s.assistant("所有极性下都没有数据文件，请检查数据目录。", "No data files for any polarity; check the data folder.")
        return
    for polarity in polarities:
        p = pipeline_runner.polarity_paths(s.config["data_root"], polarity, s.config["library_dir"])
        if not os.path.exists(p["master_table_csv"]):
            s.assistant(f"极性 {polarity} 还没有注释结果，请先运行 `注释`。",
                        f"No annotation for {polarity}; run `annotation` first.")
            return
    s.cache["_stage_span"] = (38, 55)
    s.progress(39, "quantify")
    s.assistant("开始定量（仅定量）……", "Running quantification (only) ...")
    paths = dict(s.cache.get("paths") or {})
    sample_names = None
    for polarity in polarities:
        p = pipeline_runner.run_quantification_for_polarity(
            s.config["data_root"], polarity, s.config["library_dir"],
            params=s.config["params"], log=s.log)
        paths[polarity] = p
        sample_names = sample_names or p["sample_names"]
    s.cache["paths"] = paths
    s.cache["sample_names"] = sample_names
    s.progress(55)
    s.assistant("定量完成。可继续 `质控`。", "Quantification done. Next: `qc`.")


def _ensure_paths(s):
    if "paths" not in s.cache:
        paths = _resolve_paths(s)
        missing = [pol for pol, p in paths.items() if not os.path.exists(p["quant_table_csv"])]
        if missing:
            s.assistant(f"这些极性还没有定量结果 {missing}，请先运行 `流程`。",
                        f"No quantification yet for {missing}; run `pipeline` first.")
            return None
        s.cache["paths"] = paths
        first = paths[s.config["polarities"][0]]
        _, mzml = pipeline_runner.validate_polarity_inputs(
            s.config["data_root"], s.config["polarities"][0], s.config["library_dir"])
        s.cache["sample_names"] = [os.path.splitext(os.path.basename(m))[0] for m in mzml]
    return s.cache["paths"]


def stage_qc(s):
    paths = _ensure_paths(s)
    if not paths:
        return
    s.cache["_stage_span"] = (45, 62)
    s.progress(46, "qc")
    s.assistant("开始质控（总峰面积校正 + CV + PCA + PLS-DA）……",
                "Running QC (area correction + CV + PCA + PLS-DA) ...")
    sample_names = s.cache["sample_names"]
    groups = s.config.get("groups") or infer_groups(sample_names)
    s.config["groups"] = groups
    if not s.config.get("case") or not s.config.get("control"):
        s.config["case"], s.config["control"] = pick_case_control(groups)
    if not resolve_comparisons(s):
        s.assistant("未能自动确定病例/对照，请在左侧设置对照，如：对照 B vs A。",
                    "Could not auto-pick case/control; set a comparison, e.g. 'B vs A'.")

    qc_results = {}
    for polarity, p in paths.items():
        qc_dir = os.path.join(s.config["data_root"], f"QC_{polarity}")
        qc_results[polarity] = run_qc(quant_csv=p["quant_table_csv"], output_dir=qc_dir,
                                      sample_names=sample_names, polarity=polarity, )
    s.cache["qc_results"] = qc_results
    s.cache["combined_raw"] = build_combined_corrected_table(qc_results)

    stats = [summarize_annotations(p["master_table_csv"], pol)
             for pol, p in paths.items() if os.path.exists(p["master_table_csv"])]
    s.cache["annotation_overview"] = combine_polarity_stats(stats)

    s.progress(62)
    s.emit_state()
    s.assistant("质控完成。可以调整筛选后运行 `差异分析` / run analysis。",
                "QC done. Adjust filters then `run analysis`.")


def stage_analysis(s):
    if s.cache.get("combined_raw") is None:
        s.assistant("请先运行 `质控` / run qc。", "Run `qc` first.")
        return
    comparisons = resolve_comparisons(s)
    if not comparisons:
        s.assistant("请先设置至少一组对照（如：对照 B vs A）。", "Set at least one comparison first (e.g. 'B vs A').")
        return
    s.cache["_stage_span"] = (62, 80)
    s.progress(63, "analysis")
    n = len(comparisons)
    s.assistant(f"正在对 {n} 组对照做差异分析（火山图/热图/脂肪酸链/气泡图/随机森林/LIPID MAPS 通路）……",
                f"Differential analysis over {n} comparison(s) (volcano/heatmap/FA/bubble/RF/LIPID MAPS pathways) ...")
    msgs = s.msgs()
    base = s.cache["combined_raw"].copy()
    base = apply_annotation_filter(base, s.config.get("annotation_filter"), msgs, s.log)
    base = apply_cv_filter(base, s.cache["sample_names"], s.config["groups"], s.config.get("cv_filter"), msgs, s.log)
    base = apply_class_rt_rules(base, s.config.get("class_rules"), s.log)

    output_dir = s.config.get("output_dir") or os.path.join(s.config["data_root"], "agent_report")
    analysis_dir = os.path.join(output_dir, "analysis")
    if not len(base):
        s.assistant("筛选后没有剩余特征，请放宽条件。", "No features left after filtering; relax the filters.")
        return

    pol_tag = "_".join(s.config["polarities"])
    analyses = []
    for i, comp in enumerate(comparisons):
        case, control = comp["case"], comp["control"]
        lo = 63 + (80 - 63) * i / max(1, n)
        s.progress(lo, f"analysis {case} vs {control}")
        s.assistant(f"[{i+1}/{n}] 差异分析 {case} vs {control} ……", f"[{i+1}/{n}] {case} vs {control} ...")
        try:
            analysis = lipid_analysis.run_analysis(
                base.copy(), s.cache["sample_names"], s.config["groups"], case, control,
                analysis_dir, tag=f"{case}_vs_{control}_{pol_tag}", language=s.lang)
            analyses.append(analysis)
        except Exception as e:
            s.assistant(f"对照 {case} vs {control} 分析失败：{e}", f"Comparison {case} vs {control} failed: {e}")

    if not analyses:
        s.assistant("所有对照分析均失败，请检查分组。", "All comparisons failed; check the grouping.")
        return
    s.cache["analyses"] = analyses
    s.cache["analysis"] = analyses[0]  # primary, kept for backward compatibility / 主对照，保留兼容
    s.cache["output_dir"] = output_dir
    s.progress(80)
    summary = "；".join(f"{a['case']}vs{a['control']}: {a['n_significant']}显著" for a in analyses)
    summary_en = "; ".join(f"{a['case']}vs{a['control']}: {a['n_significant']} sig" for a in analyses)
    s.assistant(f"差异分析完成（{n} 组对照）。{summary}。继续 `生成结论` / run conclusions。",
                f"Done ({n} comparison(s)). {summary_en}. Next: `run conclusions`.")


def stage_conclusions(s):
    analyses = s.cache.get("analyses") or ([s.cache["analysis"]] if s.cache.get("analysis") else [])
    if not analyses:
        s.assistant("请先运行 `差异分析` / run analysis。", "Run `analysis` first.")
        return
    s.progress(82, "conclusions")
    s.assistant("正在生成生物学结论（含 LIPID MAPS 代谢通路）……",
                "Generating biological conclusions (incl. LIPID MAPS pathways) ...")
    # Spread the conclusion budget across comparisons; keep the total bounded for verification.
    # 把结论额度分摊到各组对照；总数设上限以控制验证子智能体数量。
    n = len(analyses)
    max_total = max(s.config.get("max_conclusions", 5), 4)
    per = s.config.get("max_conclusions", 5) if n == 1 else max(2, -(-max_total // n))
    all_cons = []
    for a in analyses:
        cons = generate_conclusions(
            a, s.cache.get("annotation_overview"), s.llm(), s.lang,
            max_conclusions=per, log=s.log, sample_context=s.config.get("sample_context"))
        for c in cons:
            c["comparison"] = f"{a['case']} vs {a['control']}"
            c["case"], c["control"] = a["case"], a["control"]
        all_cons.extend(cons)
    all_cons = all_cons[:max_total] if n > 1 else all_cons
    for i, c in enumerate(all_cons, start=1):
        c["id"] = i
    s.cache["conclusions"] = all_cons
    s.progress(88)
    tag = (lambda c: f"({c['comparison']}) " if n > 1 and c.get("comparison") else "")
    listing = "\n".join(f"[{c['id']}] {tag(c)}{c['statement']}" for c in all_cons)
    s.assistant("已生成以下结论，请确认：\n" + listing +
                "\n\n回复 `确认` 进入验证，或 `删除 2,4` 去掉某几条，或 `重新生成`。",
                "Generated conclusions, please confirm:\n" + listing +
                "\n\nReply `confirm` to verify, `drop 2,4` to remove some, or `regenerate`.")


def stage_verify(s):
    conclusions = s.cache.get("conclusions")
    if not conclusions:
        s.assistant("还没有结论可验证，请先 `生成结论`。", "No conclusions to verify; run `conclusions` first.")
        return
    s.progress(89, "verify")
    s.assistant(f"启动 {len(conclusions)} 个验证子智能体……", f"Launching {len(conclusions)} verification sub-agents ...")
    analyses = s.cache.get("analyses") or [s.cache["analysis"]]
    # Union of every comparison's differential table so each conclusion finds its lipids' evidence.
    # 合并各组对照的差异表，使每条结论都能检索到相关脂质的证据。
    diff_tables = [a.get("differential_table") for a in analyses if a.get("differential_table") is not None]
    combined_diff = pd.concat(diff_tables, ignore_index=True) if diff_tables else None
    verifications = verify_conclusions(
        conclusions, combined_diff, s.llm(), s.lang, log=s.log,
        sample_context=s.config.get("sample_context"))
    s.cache["verifications"] = verifications

    output_dir = s.cache.get("output_dir") or os.path.join(s.config["data_root"], "agent_report")
    md_path = os.path.join(output_dir, "LipidIN_report.md")
    report_params = {
        "language": s.lang, "polarities": ",".join(s.config["polarities"]),
        "case": s.config["case"], "control": s.config["control"],
        "llm_provider": s.config["llm"].get("provider"),
        "annotation_filter": s.config.get("annotation_filter"),
        "cv_filter": s.config.get("cv_filter"), "class_rules": s.config.get("class_rules"),
        **{k: s.config["params"][k] for k in ("MS2_filter", "ppm_err1", "da_err1", "ppm_err2", "da_err2")},
    }
    report_params["comparisons"] = _comparisons_label(resolve_comparisons(s))
    report_writer.build_report(
        md_path=md_path, language=s.lang, params=report_params,
        annotation_overview=s.cache.get("annotation_overview"), qc_results=s.cache.get("qc_results"),
        analyses=analyses, conclusions=conclusions, verifications=verifications,
        extra_paths={"output_dir": output_dir}, sample_context=s.config.get("sample_context"))
    summary = report_writer.build_short_summary(
        s.lang, s.cache.get("annotation_overview"), s.cache.get("qc_results"),
        analyses, verifications)
    s.progress(100)
    s.assistant("报告已生成：\n" + md_path + "\n\n" + summary,
                "Report written:\n" + md_path + "\n\n" + summary)


def stage_full(s):
    """The whole workflow: annotation -> quantification -> QC -> differential -> conclusions.

    全流程：注释 -> 定量 -> 质控 -> 差异分析 -> 生成结论。（验证并出报告需你确认结论后再进行）
    """
    polarities = _resolve_run_polarities(s)
    if not polarities:
        s.assistant("所有极性下都没有数据文件，请检查数据目录。", "No data files for any polarity; check the data folder.")
        return
    # Only run the polarities that actually have data (so a neg-only dataset won't fail on pos).
    # 只运行真正有数据的极性（这样仅负离子的数据集不会在 pos 上失败）。
    s.config["polarities"] = polarities
    for polarity in polarities:
        p = pipeline_runner.polarity_paths(s.config["data_root"], polarity, s.config["library_dir"])
        if not os.path.exists(p["quant_table_csv"]):
            s.cache["_stage_span"] = (5, 50)
            pipeline_runner.run_pipeline_for_polarity(
                s.config["data_root"], polarity, s.config["library_dir"],
                params=s.config["params"], log=s.log)
    _ensure_paths(s)
    stage_qc(s)
    stage_analysis(s)
    stage_conclusions(s)


STAGES = {
    "annotation": stage_annotation,
    "quantify": stage_quantify,
    "qc": stage_qc,
    "analysis": stage_analysis,
    "conclusions": stage_conclusions,
    "verify": stage_verify,
    "full": stage_full,
}


# ------------------------------------------------------------------ #
# Message interpretation: deterministic commands first, then LLM for free text
# 消息解析：先走确定性命令，自由文本再交给大模型
# ------------------------------------------------------------------ #
def interpret(s, text):
    raw = text.strip()
    low = raw.lower()

    # language toggle
    if low in ("lang en", "english", "英文"):
        s.config["language"] = "en"; s.assistant("语言已切换为英文。", "Switched to English."); s.emit_state(); return
    if low in ("lang zh", "chinese", "中文"):
        s.config["language"] = "zh"; s.assistant("已切换为中文。", "Switched to Chinese."); s.emit_state(); return

    # An active multi-turn wizard (LLM config or run confirmation) takes priority.
    # 激活中的多轮向导（大模型配置或运行确认）优先处理。
    if s.pending:
        if s.pending == "confirm_run":
            _confirm_run_step(s, raw, low)
        else:
            _llm_wizard_step(s, raw)
        return
    if re.search(r"配置大模型|配置模型|设置大模型|设置模型|接入大模型|配置\s*api|设置\s*api|"
                 r"configure\s+(the\s+)?(llm|model)|llm\s*config|set\s*up\s*(the\s+)?(llm|model)", low):
        _start_llm_wizard(s)
        return

    if low in ("help", "帮助", "?"):
        _reply_help(s); return
    if low in ("show", "config", "配置", "状态"):
        s.push(s._state_event()); s.assistant(_config_text(s), _config_text(s)); return

    # "set data/library folder" intent but no path supplied -> tell the user the exact command.
    # 表达了“设置数据/谱库目录”的意图但没给路径 -> 直接告诉用户确切的命令格式。
    _has_path = bool(re.search(r"[A-Za-z]:[\\/]|[\\/].+[\\/]|[:：]\s*\S", raw))
    if not _has_path and re.search(r"(数据(目录|地址|路径|文件夹)|原始数据|data\s*(root|dir|folder|path))", low):
        s.assistant("请把路径写在同一句里，例如：\n数据目录 D:\\bio_inf\\LipidIN_github\\common_demo_raw",
                    "Include the path in the same message, e.g.:\ndata root D:\\bio_inf\\LipidIN_github\\common_demo_raw")
        return
    if not _has_path and re.search(r"(谱库(目录|地址|路径|文件夹)|谱图库|library\s*(dir|folder|path))", low):
        s.assistant("请把路径写在同一句里，例如：\n谱库目录 D:\\bio_inf\\LipidIN_github\\library",
                    "Include the path in the same message, e.g.:\nlibrary dir D:\\bio_inf\\LipidIN_github\\library")
        return

    # confirmation flow
    m = re.search(r"(drop|删除|去掉|去除)\s*([\d,，\s]+)", low)
    if m:
        ids = {int(x) for x in re.split(r"[,\s，]+", m.group(2)) if x.strip().isdigit()}
        conclusions = [c for c in (s.cache.get("conclusions") or []) if c["id"] not in ids]
        s.cache["conclusions"] = conclusions
        s.assistant(f"已删除 {sorted(ids)}，剩余 {len(conclusions)} 条。回复 `确认` 继续验证。",
                    f"Dropped {sorted(ids)}, {len(conclusions)} left. Reply `confirm` to verify.")
        return
    if low in ("regenerate", "重新生成", "再生成"):
        s.enqueue("conclusions"); return
    if low in ("confirm", "确认", "接受", "yes", "y"):
        s.enqueue("verify"); return

    # run actions
    action = _match_action(low)
    patch, notes = _parse_patch(s, raw, low)

    if patch:
        deep_merge(s.config, patch)
        if any(k in patch for k in CONFIRM_SENSITIVE_KEYS):
            s.confirmed = False  # directory/polarity/library/grouping changed -> re-confirm before running
        s.emit_state()
        s.assistant("已更新：" + "; ".join(notes), "Updated: " + "; ".join(notes))

    if action:
        if s.busy:
            s.assistant("正在运行中，请等当前步骤完成。", "A stage is running; please wait.")
            return
        # Every run must be confirmed against the full parameter list first.
        # 每次运行前都要先对完整参数清单进行确认。
        if not s.confirmed:
            _send_preflight(s, action)
        else:
            s.enqueue(action)
        return

    if patch:
        return

    # nothing matched -> LLM free-form parse
    llm = s.llm()
    if llm.is_enabled():
        _llm_interpret(s, raw)
    else:
        s.assistant("我需要先接入大模型才能理解这类自然语言。请回复『配置大模型』，我会一步步带你配置（提供方/密钥/模型）；"
                    "配置好后就能用自然语言了（如“FA类RT严格，其它宽松”）。也可以直接用命令：运行流程/质控/差异分析/生成结论/确认。",
                    "I need an LLM to understand free-form requests. Reply 'configure model' and I'll walk you through it "
                    "(provider/key/model); afterwards you can use natural language (e.g. 'FA strict RT, others loose'). "
                    "Or use commands: run pipeline/qc/analysis/conclusions/confirm.")


def _start_llm_wizard(s):
    s.pending = "llm_provider"
    s.assistant("好的，先配置大模型。请选择提供方：anthropic / openai / deepseek / custom（或回复 none 关闭大模型）。",
                "Let's configure the LLM. Choose a provider: anthropic / openai / deepseek / custom (or 'none' to disable).")


def _is_skip(text):
    return text.strip().lower() in ("默认", "default", "跳过", "skip", "-", "回车", "空")


def _llm_wizard_step(s, raw):
    step = s.pending
    if step == "llm_provider":
        provider = raw.strip().lower()
        if provider in ("none", "关闭", "取消", "cancel"):
            s.config["llm"] = {"provider": "none", "api_key": None, "model": None, "base_url": None}
            s.pending = None
            s.assistant("已关闭大模型，将使用离线规则引擎。", "LLM disabled; using the offline rule engine.")
            s.emit_state(); return
        if provider not in ("anthropic", "openai", "codex", "deepseek", "custom"):
            s.assistant("请从 anthropic / openai / deepseek / custom / none 中选择。",
                        "Please choose one of anthropic / openai / deepseek / custom / none.")
            return
        s.config["llm"] = {"provider": provider, "api_key": None, "model": None, "base_url": None}
        s.pending = "llm_key"
        s.assistant("请输入 API 密钥（回复『跳过』则从环境变量读取）。",
                    "Enter the API key (reply 'skip' to read it from the environment).")
        return
    if step == "llm_key":
        s.config["llm"]["api_key"] = None if _is_skip(raw) else raw.strip()
        s.pending = "llm_model"
        s.assistant("请输入模型名称（回复『默认』使用该提供方的默认模型）。",
                    "Enter the model name (reply 'default' to use the provider default).")
        return
    if step == "llm_model":
        s.config["llm"]["model"] = None if _is_skip(raw) else raw.strip()
        if s.config["llm"]["provider"] == "custom":
            s.pending = "llm_baseurl"
            s.assistant("请输入 Base URL（自定义网关地址，需兼容 OpenAI /chat/completions）。",
                        "Enter the Base URL (custom gateway, OpenAI /chat/completions compatible).")
            return
        _finish_llm_wizard(s)
        return
    if step == "llm_baseurl":
        s.config["llm"]["base_url"] = None if _is_skip(raw) else raw.strip()
        _finish_llm_wizard(s)
        return


def _finish_llm_wizard(s):
    s.pending = None
    llm = s.llm()
    s.emit_state()
    if llm.is_enabled():
        s.assistant(f"大模型已就绪（{llm.provider} · {llm.model}）。现在你可以用自然语言下达要求，"
                    f"例如“FA类RT严格，其它宽松”“只保留MS2且等级1的注释，病例B对照A，然后一键运行”。",
                    f"LLM ready ({llm.provider} · {llm.model}). You can now use natural language, e.g. "
                    f"'FA strict RT, others loose' or 'keep only MS2 grade-1 annotations, case B control A, then run all'.")
    else:
        s.assistant("大模型尚未启用：缺少密钥或模型。可回复『配置大模型』重新配置，或用命令继续。",
                    "LLM not enabled: missing key or model. Reply 'configure model' to retry, or use commands.")


# Short, command-like messages (suggestion buttons, quick typing) map directly to an action.
# 简短的命令式消息（建议按钮、快速输入）直接映射到某个动作。
_EXACT_ACTIONS = {
    "full": {"运行全流程", "全流程", "运行流程", "完整流程", "整个流程", "full workflow", "full", "run all"},
    "annotation": {"注释", "仅注释", "定性", "注释（仅定性注释）", "annotation", "annotate"},
    "quantify": {"定量", "仅定量", "定量（仅定量）", "quantify", "quantification"},
    "qc": {"质控", "仅质控", "质量控制", "质控（仅质控）", "qc"},
    "analysis": {"差异分析", "仅差异分析", "差异分析（仅差异分析）", "analysis", "differential"},
    "conclusions": {"生成结论", "结论", "conclusions"},
    "verify": {"验证", "验证并出报告", "verify", "report"},
}

# A run intent needs an explicit run verb next to the keyword, so a sentence merely *mentioning*
# a stage (e.g. "对于注释结果我只看…") does NOT trigger that stage.
# 运行意图要求关键词旁边有明确的动词，因此只是"提到"某步骤的句子（如"对于注释结果我只看…"）不会触发。
_RUN_VERB = r"(运行|重新|重跑|开始|执行|跑|做|生成|进行|re-?run|rerun|run|start|redo)"
_VERB_ACTIONS = [
    ("full", r"(运行|开始|执行|跑)\s*全流程|full\s*workflow"),
    ("annotation", _RUN_VERB + r"\s*(注释|定性)|re-?run\s+annotation|\bannotate\b|\brun\s+annotation\b"),
    ("quantify", _RUN_VERB + r"\s*定量|run\s+quantif"),
    ("qc", _RUN_VERB + r"\s*(质控|质量控制)|run\s+qc"),
    ("analysis", _RUN_VERB + r"\s*(差异分析|分析)|重新分析|重跑分析|run\s+(the\s+)?analysis|run\s+differential"),
    ("conclusions", _RUN_VERB + r"\s*结论|run\s+conclusion"),
    ("verify", r"(运行|开始|执行)?\s*验证|run\s+verify|生成\s*报告|出\s*报告"),
]


def _match_action(low):
    t = low.strip()
    for action, words in _EXACT_ACTIONS.items():
        if t in words:
            return action
    for action, pattern in _VERB_ACTIONS:
        if re.search(pattern, low):
            return action
    return None


def _grab_abs_path(text):
    """Extract the first absolute path (Windows drive, UNC, or POSIX) from a message, or None.

    从消息中抽取第一个绝对路径（Windows 盘符、UNC 或 POSIX），无则 None。支持带引号与无空格粘贴。
    """
    if not text:
        return None
    # Quoted path wins (may contain spaces).
    qm = re.search(r'["\']([A-Za-z]:[\\/][^"\']+|\\\\[^"\']+|/[^"\']+)["\']', text)
    if qm:
        return qm.group(1).strip()
    # Unquoted: drive path, UNC, or a rooted POSIX path; stop at whitespace.
    m = re.search(r'([A-Za-z]:[\\/]\S+|\\\\\S+|/(?:[^\s/]+/)+[^\s/]*)', text)
    if m:
        return m.group(1).strip().strip('"').strip("'").rstrip("。.,，；;、)）」』】")
    return None


def _parse_patch(s, raw, low):
    patch = {}
    notes = []

    # --- Case/control & multi-group comparisons / 病例对照与多组对照 ---
    m = re.search(r"(case|病例)\s*[:：]?\s*([A-Za-z0-9_]+).*?(control|对照)\s*[:：]?\s*([A-Za-z0-9_]+)", raw, re.I)
    if re.search(r"all\s*pairs|两两|全部对照|所有对照|全部两两|所有两两", low):
        groups = s.config.get("groups")
        if not groups:
            _, sn = _detect_samples(s)
            groups = {x: group_of(x) for x in (sn or [])}
        pairs = enumerate_pairs(groups)
        if pairs:
            patch["comparisons"] = pairs
            patch["case"], patch["control"] = pairs[0]["case"], pairs[0]["control"]
            notes.append(f"comparisons all pairs ({len(pairs)})")
        else:
            notes.append("all-pairs requested but no groups detected yet")
    elif m:
        patch["case"], patch["control"] = m.group(2), m.group(4)
        patch["comparisons"] = [{"case": m.group(2), "control": m.group(4)}]
        notes.append(f"case/control {m.group(2)}/{m.group(4)}")
    else:
        # One or more "X vs Y" pairs (semicolon/comma separated) -> multi-group comparisons.
        # 一个或多个 "X vs Y"（分号/逗号分隔）-> 多组对照。
        vs_pairs = re.findall(r"\b([A-Za-z0-9_]+)\s*(?:vs|VS)\s*([A-Za-z0-9_]+)\b", raw)
        if len(vs_pairs) >= 2:
            comps = [{"case": a, "control": b} for a, b in vs_pairs]
            patch["comparisons"] = comps
            patch["case"], patch["control"] = comps[0]["case"], comps[0]["control"]
            notes.append("comparisons " + _comparisons_label(comps))
        elif len(vs_pairs) == 1:
            a, b = vs_pairs[0]
            patch["case"], patch["control"] = a, b
            patch["comparisons"] = [{"case": a, "control": b}]
            notes.append(f"case/control {a}/{b}")

    m = re.search(r"(grade|等级)\s*[:：]?\s*([\d,，\s]+)", low)
    if m:
        digits = [d for d in re.split(r"[,\s，]+", m.group(2)) if d in GRADE_BY_DIGIT]
        if digits:
            patch.setdefault("annotation_filter", {})["allowed_grades"] = [GRADE_BY_DIGIT[d] for d in digits]
            notes.append("grades " + ",".join(digits))

    if re.search(r"(exclude|不要|排除|去掉).*(ms1[- ]?only|仅ms1|只有ms1)", low) or "ms1only off" in low:
        patch.setdefault("annotation_filter", {})["include_ms1_only"] = False; notes.append("ms1-only off")
    if re.search(r"(include|保留|加入).*(ms1[- ]?only|仅ms1)", low) or "ms1only on" in low:
        patch.setdefault("annotation_filter", {})["include_ms1_only"] = True; notes.append("ms1-only on")
    if re.search(r"(exclude|不要|排除).*(unannotated|未注释)", low):
        patch.setdefault("annotation_filter", {})["include_unannotated"] = False; notes.append("unannotated off")

    m = re.search(r"(cv)\s*[:：]?\s*(\d+(?:\.\d+)?)", low)
    if m:
        patch["cv_filter"] = {"enabled": True, "max_cv_percent": float(m.group(2))}
        notes.append(f"CV<{m.group(2)}%")
    if re.search(r"(no cv|不做cv|关闭cv|取消cv)", low):
        patch["cv_filter"] = {"enabled": False, "max_cv_percent": 30.0}; notes.append("CV filter off")

    for label, key in (("ms2 ppm", ("ppm_err2",)), ("ms1 ppm", ("ppm_err1",)),
                       ("ms2 da", ("da_err2",)), ("ms1 da", ("da_err1",)), ("ms2_filter", ("MS2_filter",))):
        m = re.search(re.escape(label).replace(r"\ ", r"\s*") + r"\s*[:：]?\s*(\d+(?:\.\d+)?)", low)
        if m:
            patch.setdefault("params", {})[key[0]] = float(m.group(1))
            notes.append(f"{label}={m.group(1)}")

    m = re.search(r"(context|背景|疾病)\s*[:：]\s*(.+)", raw)
    if m:
        patch["sample_context"] = m.group(2).strip(); notes.append("context set")

    m = re.search(r"(provider)\s+([A-Za-z]+)", low)
    if m:
        patch.setdefault("llm", {})["provider"] = m.group(2); notes.append(f"provider {m.group(2)}")
    m = re.search(r"(key|密钥)\s+(\S+)", raw)
    if m:
        patch.setdefault("llm", {})["api_key"] = m.group(2); notes.append("api key set")
    m = re.search(r"(model|模型)\s+(\S+)", raw)
    if m:
        patch.setdefault("llm", {})["model"] = m.group(2); notes.append(f"model {m.group(2)}")
    # base_url only for the LLM gateway; keep it specific so it does NOT swallow "数据地址"/"谱库地址".
    # base_url 仅指大模型网关；写得具体些，避免误吞“数据地址”/“谱库地址”。
    m = re.search(r"(base\s*url|网关地址|接口地址|网关)\s*[:：]?\s+(\S+)", raw, re.I)
    if m:
        patch.setdefault("llm", {})["base_url"] = m.group(2); notes.append("base url set")

    # Data / library folder. Extract any absolute path from the message so ALL phrasings work
    # deterministically (e.g. "地址更新为D:\...", "路径:D:\...", quoted, or a bare pasted path) —
    # never rely on the LLM for a path, which may acknowledge without returning the value.
    # 数据/谱库目录：从消息中抽取绝对路径，使各种说法都能确定性生效（如“地址更新为D:\…”、“路径:D:\…”、
    # 带引号，或直接粘贴的裸路径）——绝不依赖大模型识别路径（它可能只嘴上确认却不返回该值）。
    lib_kw = (r"library\s*dir(?:ectory)?|library\s*folder|library\s*path|lib\s*dir|"
              r"谱库目录|谱库地址|谱库路径|谱库文件夹|谱图库目录|谱图库|谱库")
    data_kw = (r"data\s*root|data\s*dir(?:ectory)?|data\s*folder|data\s*path|"
               r"数据目录|数据地址|数据路径|数据文件夹|数据根|原始数据目录|原始数据|数据位置|"
               r"数据|地址|路径|目录|文件夹|位置")
    grabbed = _grab_abs_path(raw)
    if grabbed:
        bare = raw.strip().strip('"').strip("'").strip()
        if re.search(lib_kw, raw, re.I):
            patch["library_dir"] = grabbed; notes.append(f"library dir -> {grabbed}")
        elif re.search(data_kw, raw, re.I) or bare == grabbed:
            patch["data_root"] = grabbed; notes.append(f"data root -> {grabbed}")

    # Auto ion-mode from detected files: set polarities to whichever actually have data.
    # 根据检测到的文件自动设定离子模式：把 polarities 设为真正有数据的那些极性。
    if re.search(r"自动.*(离子|极性|模式)|根据.*(检测|文件|数据).*(离子|极性|模式)|"
                 r"auto.*(polarit|ion\s*mode)|detect.*(polarit|ion\s*mode)", low):
        found = []
        for pol in ("pos", "neg"):
            _k, _files, _s = raw_converter.list_polarity_inputs(s.config["data_root"], pol, log=(lambda *a, **k: None))
            if _files:
                found.append(pol)
        if found:
            patch["polarities"] = found; notes.append("polarities auto -> " + ",".join(found))
    # Explicit polarity words. Accept many verbs (only/仅/只运行/运行/跑/用) and both en/zh tokens.
    # 显式极性词。接受多种动词与中英文写法。
    wants_pos = bool(re.search(r"(?<![a-z])pos(?![a-z])|正离子|正模式|阳离子|阳性", low))
    wants_neg = bool(re.search(r"(?<![a-z])neg(?![a-z])|负离子|负模式|阴离子|阴性", low))
    if "polarities" in patch:
        pass  # auto already decided
    elif re.search(r"both|都跑|正负|两个极性|正负都", low) or (wants_pos and wants_neg):
        patch["polarities"] = ["pos", "neg"]; notes.append("polarities pos+neg")
    elif wants_pos:
        patch["polarities"] = ["pos"]; notes.append("polarities pos")
    elif wants_neg:
        patch["polarities"] = ["neg"]; notes.append("polarities neg")

    # simple per-class RT strictness heuristic
    cls_rules = _parse_class_rules_heuristic(raw, low)
    if cls_rules:
        patch["class_rules"] = cls_rules
        notes.append("class RT rules: " + _class_rules_label(cls_rules))

    return patch, notes


KNOWN_CLASSES = ["FA", "PC", "PE", "PG", "PI", "PS", "PA", "TG", "DG", "MG", "CE", "SM",
                 "Cer", "LPC", "LPE", "CL", "SPB", "ST"]


def _parse_class_rules_heuristic(raw, low):
    if not re.search(r"strict|严格|宽松|loose|放宽", low):
        return None
    rules = {}
    for cls in KNOWN_CLASSES:
        # ASCII-letter boundaries only, so "FA类" / "FA " both match but "LPA" won't match "PA".
        # 仅用 ASCII 字母边界，使 "FA类"/"FA " 都能匹配，但 "LPA" 不会误配 "PA"。
        if re.search(rf"(?<![a-z]){re.escape(cls.lower())}(?![a-z])", low) and re.search(r"strict|严格", low):
            rules[cls] = {"rt_min": 0.85}
    if re.search(r"(other|其它|其他).*(loose|宽松|放宽)|(loose|宽松|放宽).*(other|其它|其他)", low):
        rules["*"] = {"rt_min": 0.3}
    return rules or None


def _llm_interpret(s, raw):
    llm = s.llm()
    schema_hint = (
        "You convert a lipidomics analyst's request into a JSON config patch for an analysis agent. "
        "Output ONLY JSON with optional keys:\n"
        "set (object; allowed fields: data_root (absolute path to the data folder, e.g. "
        "'D:\\\\bio_inf\\\\LipidIN_github\\\\common_demo_raw'), library_dir (absolute path to the spectral "
        "library folder), case, control, "
        "comparisons (list of {case, control} for multiple group-vs-group differentials; use this when the "
        "user wants several comparisons run together, e.g. B vs A and C vs A), "
        "sample_context, max_conclusions, polarities[list], "
        "annotation_filter{include_ms1_only(bool), include_unannotated(bool), allowed_grades[list of "
        "'Level 1 (MS2, complete fragments)'|'Level 2 (MS2, partial fragments)'|'Level 3 (MS1 inference + RT)']}, "
        "cv_filter{enabled(bool), max_cv_percent(number)}, "
        "params{MS2_filter, ppm_err1, da_err1, ppm_err2, da_err2}),\n"
        "class_rules (per-lipid-class RT strictness, e.g. {\"FA\":{\"rt_min\":0.85},\"*\":{\"rt_min\":0.9}}; "
        "rt_min is the minimum rt_fit_score in 0..1, higher = stricter RT match),\n"
        "action (one of: annotation, quantify, qc, analysis, conclusions, verify, full, none),\n"
        "reply (a short confirmation message in the user's language).\n\n"
        "IMPORTANT ROUTING RULES:\n"
        "1) Re-running 'annotation' re-does qualitative identification and is expensive. Choose "
        "action='annotation' ONLY when the user changes annotation hyperparameters "
        "(params.ppm_err1/ppm_err2/da_err1/da_err2/MS2_filter) or the spectral library. "
        "NEVER re-annotate just because the user talks about the annotation RESULTS.\n"
        "2) Keeping or discarding annotation results by reliability is FILTERING, not annotation: "
        "'reliable MS2 / 可靠的MS2' -> annotation_filter.allowed_grades must include "
        "'Level 1 (MS2, complete fragments)' (add 'Level 2 (MS2, partial fragments)' only if the user accepts "
        "partial MS2), and set include_ms1_only=false; drop the rest. "
        "'RT fully matches / RT完全符合 / RT可靠' -> set class_rules {\"*\":{\"rt_min\":0.9}} (strict RT for all classes).\n"
        "3) The scientific question / disease and which group is disease vs control -> set sample_context, "
        "and set case = the disease group, control = the control group. For multiple comparisons set "
        "'comparisons'; 'all pairs / 两两对照' means every pairwise combination of the non-QC groups.\n"
        "4) If the user wants conclusions redone AND filters or case/control changed, use action='analysis' "
        "(recomputes the differential with the new settings; conclusions follow). If ONLY the context changed "
        "and a differential analysis already exists, use action='conclusions'. Otherwise action='none'."
    )
    current = json.dumps({k: s.config[k] for k in ("data_root", "library_dir", "case", "control",
                          "comparisons", "polarities", "params", "annotation_filter", "cv_filter",
                          "class_rules", "sample_context")},
                         ensure_ascii=False)
    try:
        data = llm.chat_json(schema_hint, f"Current config: {current}\nUser request: {raw}")
    except Exception as e:
        s.assistant(f"自然语言解析失败：{e}", f"Free-form parsing failed: {e}")
        return
    patch = dict(data.get("set") or {})
    if data.get("class_rules"):
        patch["class_rules"] = data["class_rules"]
    if patch:
        deep_merge(s.config, patch)
        if any(k in patch for k in CONFIRM_SENSITIVE_KEYS):
            s.confirmed = False
        s.emit_state()
    reply = data.get("reply") or "已根据你的要求更新配置。"
    s.assistant(reply, reply)
    action = data.get("action")
    if action and action in STAGES:
        if s.busy:
            s.assistant("正在运行中，请稍候。", "A stage is running; please wait.")
        elif not s.confirmed:
            _send_preflight(s, action)
        else:
            s.enqueue(action)


ACTION_LABELS = {
    "annotation": ("注释（仅定性注释）", "annotation only"),
    "quantify": ("定量（仅定量）", "quantification only"),
    "qc": ("质控（仅质控）", "QC only"),
    "analysis": ("差异分析（仅差异分析）", "differential analysis only"),
    "conclusions": ("生成结论", "generate conclusions"),
    "verify": ("验证并出报告", "verify & report"),
    "full": ("运行全流程（定性→定量→质控→差异分析→生物学结论）", "the full workflow"),
}

# Editing any of these invalidates a prior confirmation and forces a re-confirm before running.
# 修改这些项会使之前的确认失效，运行前需要重新确认。
CONFIRM_SENSITIVE_KEYS = ("data_root", "library_dir", "polarities", "groups", "case", "control",
                          "comparisons", "params")


def _populate_preflight(s):
    """Detect samples and fill in groups / case / control so the user can confirm real values.

    检测样本并补全分组/病例/对照，让用户确认的是真实生效的值。
    """
    detected, sample_names = _detect_samples(s)  # counts mzML or vendor raw (Thermo/SCIEX) per polarity
    if sample_names:
        if not s.config.get("groups"):
            s.config["groups"] = infer_groups(sample_names)
        if not s.config.get("case") or not s.config.get("control"):
            s.config["case"], s.config["control"] = pick_case_control(s.config["groups"])
    return detected, sample_names


def _preflight_text(s, action):
    detected, sample_names = _populate_preflight(s)
    zh = s.lang == "zh"
    label = ACTION_LABELS.get(action, (action, action))[0 if zh else 1]
    cfg = s.config
    lines = []
    lines.append((f"运行前请确认参数（即将：{label}）：" if zh else f"Please confirm the parameters before {label}:"))

    lines.append((f"• 数据目录：{cfg['data_root']}" if zh else f"• Data folder: {cfg['data_root']}"))

    pol_bits = []
    warnings = []
    kinds = set()
    for polarity in cfg["polarities"]:
        info = detected.get(polarity, {"count": 0, "kind": "none"})
        n, kind = info["count"], info["kind"]
        kinds.add(kind)
        pol_bits.append(f"{polarity}={n}" + (f"({kind})" if kind in ("raw", "wiff") else ""))
        if n == 0:
            warnings.append((f"极性 {polarity} 下没有找到 mzML/.raw/.wiff 文件" if zh else f"no mzML/.raw/.wiff for {polarity}"))
    lines.append((f"• 极性与文件数：{', '.join(pol_bits)}" if zh else f"• Polarities & files: {', '.join(pol_bits)}"))
    if "raw" in kinds:
        lines.append(("• 检测到 Thermo .raw：将直接转为 pkl（不生成中间 mzML）。"
                      if zh else "• Thermo .raw detected: converted directly to pkl (no intermediate mzML)."))
    if "wiff" in kinds:
        lines.append(("• 检测到 SCIEX .wiff：将用 msconvert 转为精简的 centroid mzML 再处理。"
                      if zh else "• SCIEX .wiff detected: converted to a compact centroided mzML via msconvert."))

    lines.append((f"• 谱库目录：{cfg['library_dir']}" if zh else f"• Library folder: {cfg['library_dir']}"))
    lib_bits = []
    for polarity in cfg["polarities"]:
        lib = os.path.join(cfg["library_dir"], f"{polarity}_common.pkl")
        ok = os.path.exists(lib)
        lib_bits.append(f"{polarity}_common.pkl {'✓' if ok else '✗'}")
        if not ok:
            warnings.append((f"缺少谱库 {polarity}_common.pkl" if zh else f"missing library {polarity}_common.pkl"))
    lines.append((f"• 谱库文件：{', '.join(lib_bits)}" if zh else f"• Library files: {', '.join(lib_bits)}"))

    if sample_names:
        lines.append((f"• 样本（{len(sample_names)}）：{', '.join(sample_names)}"
                      if zh else f"• Samples ({len(sample_names)}): {', '.join(sample_names)}"))
        groups = cfg.get("groups") or {}
        counts = {}
        for g in groups.values():
            counts[g] = counts.get(g, 0) + 1
        grp = ", ".join(f"{g}×{n}" for g, n in sorted(counts.items()))
        lines.append((f"• 分组：{grp}" if zh else f"• Groups: {grp}"))
        comps = resolve_comparisons(s)
        comp_txt = _comparisons_label(comps)
        lines.append((f"• 对照（{len(comps)} 组）：{comp_txt}" if zh else f"• Comparisons ({len(comps)}): {comp_txt}"))
    else:
        warnings.append(("未检测到样本，请检查数据目录" if zh else "no samples detected; check the data folder"))

    af = cfg["annotation_filter"]
    filt = (f"等级={_grades_label(af.get('allowed_grades'))}, MS1-only={'on' if af.get('include_ms1_only', True) else 'off'}, "
            f"CV={'<%.0f%%' % cfg['cv_filter']['max_cv_percent'] if cfg['cv_filter'].get('enabled') else 'off'}, "
            f"class_rules={_class_rules_label(cfg.get('class_rules'))}" )
    lines.append((f"• 注释/定量筛选：{filt}" if zh else f"• Filters: {filt}"))
    lines.append((f"• 大模型：{cfg['llm'].get('provider')}" if zh else f"• LLM: {cfg['llm'].get('provider')}"))

    if warnings:
        lines.append(("⚠ 注意：" + "；".join(warnings)) if zh else ("⚠ Warnings: " + "; ".join(warnings)))

    lines.append("")
    lines.append(("回复『确认』开始；或直接修改后我会再让你确认，例如："
                  "only pos / 病例 B 对照 A / 对照 B vs A, C vs A / 两两对照 / 数据目录 D:\\path / 谱库目录 D:\\lib / 等级 1,2 / 排除 ms1-only。回复『取消』放弃。"
                  if zh else
                  "Reply 'confirm' to start; or edit first (e.g. only pos / case B control A / data root D:\\path / "
                  "library dir D:\\lib / grades 1,2 / exclude ms1-only) and I'll ask again. Reply 'cancel' to abort."))
    return "\n".join(lines)


def _send_preflight(s, action):
    s.pending = "confirm_run"
    s.pending_action = action
    s.emit_state()
    text = _preflight_text(s, action)
    s.assistant(text, text)


def _confirm_run_step(s, raw, low):
    if low in ("确认", "confirm", "yes", "y", "开始", "start", "go", "ok", "好", "好的", "可以"):
        action = s.pending_action or "all"
        s.pending = None
        s.confirmed = True
        if s.busy:
            s.assistant("正在运行中，请稍候。", "A stage is running; please wait.")
            return
        s.assistant(f"参数已确认，开始：{ACTION_LABELS.get(action, (action,))[0]}", f"Confirmed, starting: {action}")
        s.enqueue(action)
        return
    if low in ("取消", "cancel", "算了", "no", "n"):
        s.pending = None
        s.pending_action = None
        s.assistant("已取消。你可以修改参数后重新发起运行。", "Cancelled. Edit parameters and start again.")
        return
    if low in ("help", "帮助", "?"):
        _reply_help(s)
        _send_preflight(s, s.pending_action)
        return

    # Otherwise treat the message as edits (structured, or LLM free-form) and re-show preflight.
    # 否则把消息当作修改（结构化或大模型自由文本），再重新展示确认信息。
    new_action = _match_action(low)
    patch, notes = _parse_patch(s, raw, low)
    handled = False
    if patch:
        deep_merge(s.config, patch)
        handled = True
        s.assistant("已更新：" + "; ".join(notes), "Updated: " + "; ".join(notes))
    elif s.llm().is_enabled():
        _llm_interpret(s, raw)
        handled = True
    if new_action and new_action in STAGES:
        s.pending_action = new_action
    if not handled:
        s.assistant("没识别到修改。回复『确认』开始，或『取消』放弃，或直接说明要改的参数。",
                    "No change recognized. Reply 'confirm' to start, 'cancel' to abort, or state the parameter to change.")
    _send_preflight(s, s.pending_action)


def _config_text(s):
    c = s.config
    af = c["annotation_filter"]
    return (f"极性={','.join(c['polarities'])} | 病例/对照={c.get('case')}/{c.get('control')} | "
            f"等级={_grades_label(af.get('allowed_grades'))} | MS1-only={'on' if af.get('include_ms1_only', True) else 'off'} | "
            f"CV={'<%.0f%%' % c['cv_filter']['max_cv_percent'] if c['cv_filter'].get('enabled') else 'off'} | "
            f"class_rules={_class_rules_label(c.get('class_rules'))} | provider={c['llm'].get('provider')}")


def _reply_help(s):
    s.assistant(
        "我可以运行整个脂质分析流程。建议第一步先『配置大模型』，配置后即可用自然语言下达要求。\n"
        "• 配置大模型：回复“配置大模型”，我会带你填写 提供方/密钥/模型\n"
        "• 分步运行：注释（仅定性）/ 定量（仅定量）/ 质控（仅质控）/ 差异分析（仅差异分析）/ 生成结论 / 确认（验证并出报告）\n"
        "• 运行全流程：定性→定量→质控→差异分析→生物学结论 一次跑完\n"
        "• 改参数：病例 B 对照 A、等级 1,2、排除 ms1-only、cv 25、ms2 ppm 20、背景：肝癌 vs 癌旁\n"
        "• 个性化（配置大模型后）：如“FA类RT严格，其它类别宽松”\n"
        "• 每次运行前我都会让你确认目录/极性/谱库/分组等参数；随时可改参数并重跑某一步",
        "I run the whole lipid analysis. First, it's best to 'configure model'; then you can use natural language.\n"
        "• configure model: reply 'configure model' and I'll ask for provider/key/model\n"
        "• step by step: annotation (qualitative) / quantify / qc / analysis / conclusions / confirm (verify & report)\n"
        "• full workflow: annotation -> quantify -> qc -> analysis -> conclusions in one go\n"
        "• edit: case B control A, grades 1,2, exclude ms1-only, cv 25, ms2 ppm 20, context: HCC vs adjacent\n"
        "• personalized (after configuring an LLM): e.g. 'FA strict RT, others loose'\n"
        "• before every run I ask you to confirm folder/polarity/library/grouping; change any setting anytime and re-run a stage")


# ------------------------------------------------------------------ #
# Flask app
# ------------------------------------------------------------------ #
app = Flask(__name__)


@app.route("/")
def index():
    if LOGO_DATA_URI:
        logo_html = f'<img class="logo" src="{LOGO_DATA_URI}" alt="LipidIN"/>'
    else:
        logo_html = '<div class="logo logo-fallback"></div>'
    html = INDEX_HTML.replace("__LOGO__", logo_html).replace("__FAVICON__", LOGO_DATA_URI)
    return Response(html, mimetype="text/html")


@app.route("/api/message", methods=["POST"])
def api_message():
    text = (request.get_json(force=True) or {}).get("text", "")
    if text.strip():
        threading.Thread(target=interpret, args=(SESSION, text), daemon=True).start()
    return jsonify({"ok": True})


@app.route("/api/stream")
def api_stream():
    def gen():
        q = SESSION.subscribe()
        try:
            yield "retry: 3000\n\n"
            while True:
                try:
                    item = q.get(timeout=15)
                    yield f"data: {json.dumps(item, ensure_ascii=False)}\n\n"
                except queue.Empty:
                    yield ": keepalive\n\n"
        finally:
            SESSION.unsubscribe(q)
    return Response(gen(), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no", "Connection": "keep-alive"})


@app.route("/api/state")
def api_state():
    return jsonify(SESSION._state_event())


INDEX_HTML = r"""<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<link rel="icon" href="__FAVICON__"/>
<title>LipidIN Agent</title>
<style>
:root{
  --bg:#eef0f4; --panel:#ffffff; --ink:#0f1220; --muted:#8a90a2; --line:#e6e8ef;
  --accent:#6d5efc; --accent2:#9b8bff; --user:#111427; --ok:#18c29c; --warn:#e8a13a;
  --shadow:0 10px 34px rgba(20,22,40,.10);
}
@media (prefers-color-scheme: dark){
  :root{ --bg:#08090d; --panel:#14161e; --ink:#eef0f6; --muted:#8b90a3; --line:#242736;
    --user:#eef0f6; --shadow:0 12px 44px rgba(0,0,0,.5); }
}
*{box-sizing:border-box}
html,body{height:100%;margin:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",Roboto,Helvetica,Arial,sans-serif;
  background:var(--bg);color:var(--ink);-webkit-font-smoothing:antialiased;
  display:flex;align-items:center;justify-content:center;overflow:hidden}
/* 16:9 stage / 16:9 舞台 */
.stage{width:min(100vw,calc(100vh * 16 / 9));height:min(100vh,calc(100vw * 9 / 16));
  aspect-ratio:16/9;display:flex;flex-direction:column;padding:14px;gap:10px}
.head{display:flex;align-items:center;gap:12px;flex:none}
/* Full wide logo on a transparent background, 51px tall, auto width (whole mark shows, no crop). */
.logo{height:51px;width:auto;max-width:260px;object-fit:contain;display:block;background:transparent}
.logo-fallback{width:51px;height:51px;border-radius:12px;
  background:linear-gradient(135deg,var(--accent),var(--accent2));border:none}
.brand{display:flex;flex-direction:column;align-items:flex-start;gap:5px}
.subline{display:flex;align-items:center;gap:7px}
.title{font-weight:650;font-size:16px;letter-spacing:.2px}
.subtitle{color:var(--muted);font-size:11.5px;margin-top:1px}
.spacer{flex:1}
.toggle{border:2px solid var(--line);background:var(--panel);border-radius:999px;padding:5px 11px;font-size:12px;
  color:var(--muted);cursor:pointer;transition:.2s}
.toggle:hover{color:var(--ink)}
.dot{width:7px;height:7px;border-radius:50%;background:var(--ok);display:inline-block;margin-left:6px}
.dot.busy{background:var(--accent);animation:pulse 1.2s infinite}
@keyframes pulse{0%{box-shadow:0 0 0 0 rgba(109,94,252,.5)}70%{box-shadow:0 0 0 7px rgba(109,94,252,0)}100%{box-shadow:0 0 0 0 rgba(109,94,252,0)}}
/* three columns / 三栏 */
.cols{flex:1;min-height:0;display:grid;grid-template-columns:minmax(220px,1fr) minmax(0,1.7fr) minmax(230px,1.05fr);gap:10px}
.col{background:var(--panel);border:2px solid var(--line);border-radius:16px;box-shadow:var(--shadow);
  display:flex;flex-direction:column;min-height:0;overflow:hidden}
.col-title{font-size:12px;font-weight:650;color:var(--ink);padding:11px 14px 7px;letter-spacing:.3px;
  display:flex;align-items:center;gap:6px}
.col-title .en{color:var(--muted);font-weight:500;font-size:11px}
.col-body{flex:1;min-height:0;overflow-y:auto;padding:2px 12px 12px}
.sec{margin-top:8px}
.sec h4{margin:6px 2px;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.6px;font-weight:700}
.kv{display:flex;justify-content:space-between;gap:8px;padding:5px 8px;border-radius:9px;font-size:12px}
.kv:nth-child(odd){background:rgba(128,128,140,.06)}
.kv .k{color:var(--muted);flex:none}
.kv .v{color:var(--ink);text-align:right;word-break:break-all;font-weight:550}
.badge{font-size:10px;padding:1px 6px;border-radius:999px;border:1.5px solid var(--line);color:var(--muted);margin-left:5px}
.badge.raw{color:var(--warn);border-color:var(--warn)}
.badge.mzml{color:var(--ok);border-color:var(--ok)}
.grp{display:flex;align-items:center;gap:7px;padding:5px 8px;border-radius:9px;font-size:12px}
.grp:nth-child(odd){background:rgba(128,128,140,.06)}
.grp .n{margin-left:auto;color:var(--muted);font-variant-numeric:tabular-nums}
.role{font-size:10px;padding:1px 7px;border-radius:999px;font-weight:650}
.role.cc{background:rgba(109,94,252,.14);color:var(--accent)}
.role.qc{background:rgba(24,194,156,.14);color:var(--ok)}
.role.group{background:rgba(128,128,140,.12);color:var(--muted)}
.comp{display:flex;align-items:center;gap:6px;padding:5px 9px;margin:4px 0;border-radius:9px;font-size:12px;
  background:rgba(109,94,252,.08);border:1.5px solid rgba(109,94,252,.18)}
.comp b{color:var(--accent)}
.note{font-size:11px;color:var(--muted);line-height:1.5;margin:6px 2px}
.bioq{font-size:12.5px;line-height:1.55;color:var(--ink);background:rgba(128,128,140,.06);border:1.5px dashed var(--line);
  border-radius:11px;padding:9px 11px;min-height:44px;white-space:pre-wrap}
.bioq.empty{color:var(--muted)}
.mini{font-size:11px;color:var(--accent);cursor:pointer;border:1.5px solid var(--line);border-radius:999px;
  padding:3px 9px;background:var(--panel);transition:.15s;margin-top:6px;display:inline-block}
.mini:hover{border-color:var(--accent)}
/* middle chat / 中间对话 */
.mid{padding:0}
.progress{height:7px;border-radius:999px;background:var(--line);overflow:hidden;margin:12px 14px 0}
.bar{height:100%;width:0%;border-radius:999px;background:linear-gradient(90deg,var(--accent),var(--accent2));transition:width .4s ease}
.bar.idle{opacity:.35}
.chat{flex:1;overflow-y:auto;padding:14px;display:flex;flex-direction:column;gap:10px}
.msg{max-width:86%;padding:10px 13px;border-radius:15px;font-size:13.5px;line-height:1.5;white-space:pre-wrap;
  word-break:break-word;animation:pop .25s ease}
@keyframes pop{from{opacity:0;transform:translateY(6px)}to{opacity:1;transform:none}}
.assistant{align-self:flex-start;background:var(--bg);border:2px solid var(--line);border-bottom-left-radius:5px}
.user{align-self:flex-end;background:linear-gradient(135deg,var(--accent),var(--accent2));color:#fff;border-bottom-right-radius:5px}
.log{align-self:flex-start;max-width:92%;color:var(--muted);font-size:11px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  background:transparent;padding:0 4px;white-space:pre-wrap}
.suggest{display:flex;gap:6px;flex-wrap:wrap;padding:0 12px 8px}
.s{font-size:11.5px;color:var(--ink);background:var(--panel);border:2px solid var(--line);border-radius:999px;
  padding:5px 10px;cursor:pointer;transition:.15s}
.s:hover{border-color:var(--accent);color:var(--accent)}
.inputbar{display:flex;gap:9px;align-items:flex-end;background:var(--panel);border-top:2px solid var(--line);
  padding:10px 12px}
.inputbar .wrap{flex:1;display:flex;align-items:flex-end;background:var(--bg);border:2px solid var(--line);
  border-radius:15px;padding:8px 8px 8px 13px;transition:border-color .2s}
.inputbar .wrap:focus-within{border-color:var(--accent)}
textarea{flex:1;border:none;outline:none;resize:none;background:transparent;color:var(--ink);font-size:13.5px;
  font-family:inherit;max-height:100px;line-height:1.5}
.send{width:38px;height:38px;border:none;border-radius:12px;background:linear-gradient(135deg,var(--accent),var(--accent2));
  color:#fff;font-size:16px;cursor:pointer;flex:none;transition:.2s;display:flex;align-items:center;justify-content:center}
.send:disabled{opacity:.4;cursor:not-allowed}
/* right citations / 右侧引用 */
.cite-head{font-size:12px;line-height:1.5;color:var(--ink);background:rgba(109,94,252,.07);
  border:1.5px solid rgba(109,94,252,.2);border-radius:12px;padding:10px 12px;margin:6px 2px}
.cite-head .t{font-weight:650}
.cite-head .j{color:var(--muted);font-size:11.5px;margin-top:3px}
.link-btn{display:inline-flex;align-items:center;gap:5px;margin-top:8px;font-size:11.5px;color:#fff;
  background:linear-gradient(135deg,var(--accent),var(--accent2));border-radius:999px;padding:5px 12px;
  text-decoration:none}
.fmt{border:2px solid var(--line);border-radius:12px;padding:9px 11px;margin:8px 2px;position:relative}
.fmt .lbl{font-size:10.5px;font-weight:700;color:var(--accent);text-transform:uppercase;letter-spacing:.5px}
.fmt .txt{font-size:11.5px;line-height:1.55;color:var(--ink);margin-top:5px;word-break:break-word;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace;white-space:pre-wrap}
.copy{position:absolute;top:8px;right:8px;font-size:10.5px;border:1.5px solid var(--line);border-radius:8px;
  background:var(--panel);color:var(--muted);cursor:pointer;padding:2px 8px;transition:.15s}
.copy:hover{border-color:var(--accent);color:var(--accent)}
.copy.ok{color:var(--ok);border-color:var(--ok)}
.contact{font-size:12px;color:var(--ink);background:rgba(128,128,140,.06);border-radius:11px;padding:9px 11px;margin:8px 2px}
.contact a{color:var(--accent);text-decoration:none}
::-webkit-scrollbar{width:8px}::-webkit-scrollbar-thumb{background:var(--line);border-radius:8px}
@media (max-aspect-ratio:1/1){ .stage{width:100vw;height:100vh;aspect-ratio:auto} .cols{grid-template-columns:1fr} }
</style>
</head>
<body>
<div class="stage">
  <div class="head">
    <div class="brand">
      __LOGO__
      <div class="subline"><span class="subtitle" id="sub">脂质组学分析 · 对话式 · 大语言模型</span><span class="dot" id="dot"></span></div>
    </div>
    <div class="spacer"></div>
    <div class="toggle" id="langbtn">中 / EN</div>
  </div>
  <div class="cols">
    <!-- LEFT: parameters / groups / biological question -->
    <div class="col">
      <div class="col-title"><span id="lt-param">运行状态</span> <span class="en" id="lt-param-en">Live status</span></div>
      <div class="col-body" id="leftbody">
        <div class="note" id="left-empty">尚未检测到数据。请在中间对话框配置。<br/>No data yet — configure in the chat.</div>
      </div>
    </div>
    <!-- MIDDLE: chat -->
    <div class="col mid">
      <div class="progress"><div class="bar idle" id="bar"></div></div>
      <div class="chat" id="chat"></div>
      <div class="suggest" id="suggest"></div>
      <div class="inputbar">
        <div class="wrap"><textarea id="box" rows="1" placeholder="输入指令，如：运行全流程 / 病例 B 对照 A / 两两对照 / FA类RT严格"></textarea></div>
        <button class="send" id="send">➤</button>
      </div>
    </div>
    <!-- RIGHT: citations -->
    <div class="col">
      <div class="col-title"><span id="lt-cite">引用格式</span> <span class="en" id="lt-cite-en">How to cite</span></div>
      <div class="col-body" id="citebody"></div>
    </div>
  </div>
</div>
<script>
const chat=document.getElementById('chat'), bar=document.getElementById('bar'),
  box=document.getElementById('box'), send=document.getElementById('send'), dot=document.getElementById('dot'),
  suggest=document.getElementById('suggest'), leftbody=document.getElementById('leftbody'),
  citebody=document.getElementById('citebody');
let lang='zh';
const SUGG={zh:['配置大模型','运行全流程','注释','定量','质控','差异分析','两两对照','生成结论','确认','帮助'],
             en:['configure model','full workflow','annotation','quantify','qc','analysis','all pairs','conclusions','confirm','help']};
const L={zh:{param:'运行状态',cite:'引用格式',rp:'运行参数',gc:'组别与对照',bq:'生物学问题',
    data:'数据目录',lib:'谱库目录',pol:'极性/文件',filt:'注释筛选',ppm:'质量精度',llm:'大模型',
    grades:'置信等级',ms1:'仅MS1',cv:'CV筛选',cls:'类别RT规则',model:'模型',
    groups:'检测到的组别',comps:'对照组合',noComp:'尚未设置对照',
    bqEmpty:'尚未设置生物学问题。可在对话框输入：背景：肝癌 vs 癌旁',
    editBq:'✎ 设置生物学问题',setComp:'✎ 设置对照',
    multiNote:'支持多组对照：可输入“对照 B vs A, C vs A”或“两两对照”一次跑多组差异。',
    sub:'脂质组学分析 · 对话式 · 参数/对照/引用三栏',
    ph:'输入指令，如：运行全流程 / 病例 B 对照 A / 两两对照 / FA类RT严格',
    copy:'复制',copied:'已复制',openArt:'查看文章'},
  en:{param:'Live status',cite:'How to cite',rp:'Run parameters',gc:'Groups & comparisons',bq:'Biological question',
    data:'Data folder',lib:'Library',pol:'Polarity/files',filt:'Filters',ppm:'Mass accuracy',llm:'LLM',
    grades:'Grades',ms1:'MS1-only',cv:'CV filter',cls:'Class RT rules',model:'Model',
    groups:'Detected groups',comps:'Comparisons',noComp:'No comparison set yet',
    bqEmpty:'No biological question yet. In the chat type: context: HCC vs adjacent',
    editBq:'✎ Set biological question',setComp:'✎ Set comparison',
    multiNote:'Multi-group supported: type "B vs A, C vs A" or "all pairs" to run several differentials at once.',
    sub:'Lipidomics · conversational · LLM',
    ph:'Type a command, e.g. full workflow / case B control A / all pairs / FA strict RT',
    copy:'Copy',copied:'Copied',openArt:'Open article'}};

// ---- citations (real metadata: Nat. Commun. 2025, 16:4566) ----
const ART_URL='https://www.nature.com/articles/s41467-025-59683-5';
const EMAIL='852152344@qq.com';
const CITES=[
  ['Nature Style','Xu, H., Jiang, T., Lin, Y. et al. LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics. Nat. Commun. 16, 4566 (2025).'],
  ['Science Style','Xu, H., Jiang, T., Lin, Y. et al., LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics. Nat. Commun. 16, 4566 (2025).'],
  ['Cell Style','Xu, H. (2025) LipidIN: a comprehensive repository for flash platform independent annotation and reverse lipidomics. Nat. Commun. 16, 4566 '],
  ['ACS','Xu, H.; Jiang, T.; Lin, Y.; Zhang, L.; Yang, H.; Huang, X.; Mao, R.; Yang, Z.; Zeng, C.; Zhao, S.; Di, L.; Zhang, W.; Zeng, J.; Cai, Z.; Lin, S.-H. LipidIN: A Comprehensive Repository for Flash Platform-Independent Annotation and Reverse Lipidomics. Nat. Commun. 2025, 16, 4566.'],
  ['BibTeX','@article{Xu2025LipidIN,\n  title={LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics},\n  author={Xu, Hao and Jiang, Tianhang and Lin, Yuxiang and Zhang, Lei and Yang, Huan and Huang, Xiaoyun and Mao, Ridong and Yang, Zhu and Zeng, Changchun and Zhao, Shuang and Di, Lijun and Zhang, Wenbin and Zeng, Jun and Cai, Zongwei and Lin, Shu-Hai},\n  journal={Nature Communications},\n  year={2025},\n  volume={16},\n  number={1},\n  pages={4566},\n  doi={10.1038/s41467-025-59683-5}\n}']
];

function esc(s){return (s==null?'':String(s)).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
function copyText(t,btn){ const done=()=>{btn.textContent=L[lang].copied;btn.classList.add('ok');
    setTimeout(()=>{btn.textContent=L[lang].copy;btn.classList.remove('ok');},1400);};
  if(navigator.clipboard&&navigator.clipboard.writeText){navigator.clipboard.writeText(t).then(done,done);}
  else{const ta=document.createElement('textarea');ta.value=t;document.body.appendChild(ta);ta.select();
    try{document.execCommand('copy');}catch(e){} ta.remove();done();} }

function renderCitations(){ const T=L[lang];
  let h=`<div class="cite-head"><div class="t">LipidIN</div>`+
    `<div class="j">Nature Communications 2025, 16:4566 · Xu H. et al.</div>`+
    `<a class="link-btn" href="${ART_URL}" target="_blank" rel="noopener">🔗 ${T.openArt}</a></div>`;
  CITES.forEach((c,i)=>{ h+=`<div class="fmt"><button class="copy" data-i="${i}">${T.copy}</button>`+
    `<div class="lbl">${c[0]}</div><div class="txt">${esc(c[1])}</div></div>`; });
  h+=`<div class="contact">${lang==='zh'?'联系方式':'Contact'}: <a href="mailto:${EMAIL}">${EMAIL}</a>`+
     ` <button class="copy" style="position:static;margin-left:6px" id="cpmail">${T.copy}</button></div>`;
  citebody.innerHTML=h;
  citebody.querySelectorAll('.copy[data-i]').forEach(b=>b.onclick=()=>copyText(CITES[+b.dataset.i][1],b));
  const cm=document.getElementById('cpmail'); if(cm) cm.onclick=()=>copyText(EMAIL,cm);
  document.getElementById('lt-cite').textContent=T.cite;
}

function kv(k,v,extra){ return `<div class="kv"><span class="k">${esc(k)}</span><span class="v">${v}${extra||''}</span></div>`; }
function renderPanel(p){ const T=L[lang];
  document.getElementById('lt-param').textContent=T.param;
  if(!p){ return; }
  const emptyGuess=!p.groups.length && !(p.inputs && Object.values(p.inputs).some(x=>x.count));
  let h='';
  // Run parameters
  h+='<div class="sec"><h4>'+T.rp+'</h4>';
  h+=kv(T.data, esc(p.data_root));
  h+=kv(T.lib, esc(p.library_dir));
  let polbits=(p.polarities||[]).map(pol=>{const inf=(p.inputs||{})[pol]||{count:0,kind:'none'};
    const badge=inf.kind==='raw'?'<span class="badge raw">raw</span>':(inf.kind==='mzML'?'<span class="badge mzml">mzML</span>':'');
    return esc(pol)+'='+inf.count+badge;}).join(' &nbsp;');
  h+=kv(T.pol, polbits||'-');
  h+=kv(T.grades, esc(p.grades)) + kv(T.ms1, esc(p.ms1_only)) + kv(T.cv, esc(p.cv));
  h+=kv(T.cls, esc(p.class_rules));
  h+=kv(T.ppm, 'MS1 '+esc(p.ppm.ms1)+' / MS2 '+esc(p.ppm.ms2)+' ppm');
  h+=kv(T.llm, esc(p.provider)+(p.model?' · '+esc(p.model):''));
  h+='</div>';
  // Groups & comparisons
  h+='<div class="sec"><h4>'+T.gc+'</h4>';
  if(p.groups.length){ h+='<div class="note" style="margin-top:2px">'+T.groups+'</div>';
    p.groups.forEach(g=>{ const rc=g.role==='qc'?'qc':(g.role==='case/control'?'cc':'group');
      h+=`<div class="grp"><span class="role ${rc}">${esc(g.role)}</span><b>${esc(g.name)}</b><span class="n">×${g.n}</span></div>`; }); }
  h+='<div class="note" style="margin-top:8px">'+T.comps+'</div>';
  if(p.comparisons && p.comparisons.length){
    p.comparisons.forEach(c=>{ h+=`<div class="comp">⚖ <b>${esc(c.case)}</b> vs ${esc(c.control)}</div>`; }); }
  else{ h+='<div class="note">'+T.noComp+'</div>'; }
  h+='<div class="mini" onclick="hint(\''+(lang==='zh'?'对照 B vs A, C vs A':'B vs A, C vs A')+'\')">'+T.setComp+'</div>';
  h+='<div class="note">'+T.multiNote+'</div>';
  h+='</div>';
  // Biological question
  h+='<div class="sec"><h4>'+T.bq+'</h4>';
  if(p.context){ h+='<div class="bioq">'+esc(p.context)+'</div>'; }
  else{ h+='<div class="bioq empty">'+T.bqEmpty+'</div>'; }
  h+='<div class="mini" onclick="hint(\''+(lang==='zh'?'背景：':'context: ')+'\')">'+T.editBq+'</div>';
  h+='</div>';
  leftbody.innerHTML=h;
}
function hint(prefix){ box.value=prefix; box.focus(); auto(); }
window.hint=hint;

function renderSuggest(){ suggest.innerHTML=''; SUGG[lang].forEach(t=>{const b=document.createElement('div');
  b.className='s'; b.textContent=t; b.onclick=()=>sendMsg(t); suggest.appendChild(b);}); }
function add(cls,text){ const d=document.createElement('div'); d.className='msg '+cls; d.textContent=text;
  chat.appendChild(d); chat.scrollTop=chat.scrollHeight; }
function addLog(text){ const d=document.createElement('div'); d.className='log'; d.textContent=text;
  chat.appendChild(d); chat.scrollTop=chat.scrollHeight;
  while(chat.querySelectorAll('.log').length>400){ chat.querySelector('.log').remove(); } }
function setBusy(b){ send.disabled=b; dot.className='dot'+(b?' busy':''); bar.className='bar'+(b?'':' idle'); }
let lastPanel=null;
const es=new EventSource('/api/stream');
es.onmessage=(e)=>{ const m=JSON.parse(e.data);
  if(m.type==='assistant') add('assistant',m.text);
  else if(m.type==='log') addLog(m.text);
  else if(m.type==='progress'){ bar.style.width=(m.value||0)+'%'; }
  else if(m.type==='busy') setBusy(m.busy);
  else if(m.type==='config'){ lastPanel=m.panel; renderPanel(m.panel); }
};
function sendMsg(t){ t=(t||box.value).trim(); if(!t)return; add('user',t); box.value=''; auto();
  fetch('/api/message',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({text:t})}); }
send.onclick=()=>sendMsg();
box.addEventListener('keydown',e=>{ if(e.key==='Enter'&&!e.shiftKey){e.preventDefault();sendMsg();} });
function auto(){ box.style.height='auto'; box.style.height=Math.min(box.scrollHeight,100)+'px'; }
box.addEventListener('input',auto);
function applyLang(){ const T=L[lang];
  document.getElementById('sub').textContent=T.sub; box.placeholder=T.ph;
  document.getElementById('lt-param').textContent=T.param;
  document.getElementById('lt-param-en').textContent=(lang==='zh'?'Live status':'实时状态');
  document.getElementById('lt-cite').textContent=T.cite;
  document.getElementById('lt-cite-en').textContent=(lang==='zh'?'How to cite':'如何引用');
  renderSuggest(); renderCitations(); renderPanel(lastPanel); }
document.getElementById('langbtn').onclick=()=>{ lang=(lang==='zh')?'en':'zh'; applyLang();
  sendMsg(lang==='zh'?'中文':'english'); };
renderSuggest(); renderCitations();
add('assistant', '你好，我是 LipidIN 分析助手。左栏实时显示运行参数、组别与对照（支持多组对照）、生物学问题；右栏是多期刊引用格式。\n'
  +'建议先“配置大模型”，然后“运行全流程”，或分步：注释→定量→质控→差异分析→生成结论→确认。\n'
  +'支持 mzML，以及 Thermo .raw SCIEX .wiff 原始文件。\n'
  +'Hi, I am the LipidIN assistant. Left column: live parameters, groups & comparisons (multi-group), biological question. Right column: citation formats. Also accepts Thermo .raw / SCIEX .wiff.');
</script>
</body>
</html>
"""


def main():
    host, port = "127.0.0.1", 7860
    url = f"http://{host}:{port}"
    print(f"[INFO] LipidIN Agent UI at {url}\n[信息] LipidIN 智能体界面： {url}")
    try:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    except Exception:
        pass
    app.run(host=host, port=port, threaded=True, debug=False)


if __name__ == "__main__":
    main()
