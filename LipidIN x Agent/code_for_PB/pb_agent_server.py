# -*- coding: utf-8 -*-
"""
PB 双键定位工作流的 Agent。

结构照搬本仓库主管线 `code/agent_server.py` 的形态：Flask + SSE 事件流 + 单后台
worker 串行跑阶段 + 一个会话对象持有配置/缓存，界面也与其一致（16:9 三栏：左侧运行
状态、中间对话+进度条、右侧引用格式）。区别在于这里包的是 PB 的十一个阶段。

阶段（可单独跑也可 `全流程`）：
  1  annotation  常规注释 sample/common        -> master_table_quant.csv
  2  search      PB 衍生化检索 sample/PB       -> *_PB_diagnostic.csv
  3  localize    双键位置重建                   -> *_PB_isomers.csv / *_PB_dbpos.csv
  4  quantify    链×双键级定量 + 总量归一        -> PB_dbposition_quant.csv / PB_species_quant.csv
  5  design      用哪些文件、跟谁比（自然语言→设计，LLM）-> PB_design.json
  6  qc          质控 + 差异分析 + 位置偏好谱     -> PB_qc.csv / PB_*_diff.csv / PB_position_*.csv
  7  plots       出图                          -> figures/Fig1..Fig2 (svg/pdf/png)
  8  interpret   大模型解读质控与差异
  9  hypothesis  科学假设 + 文献检索 + LIPID MAPS 联合分析
       ↑ 假设会回给用户确认，可用自然语言反复修正（`确认` 才进入下一步）
  10 verify      每条假设一个子 agent，基于数据切片再验证（被要求优先反驳）
  11 report      最终 md：验证结果 + 注释可信度 + 强度可信度 + 是否够新 -> agent_report/

大模型：与主界面同一套配置（回复"配置大模型"走向导）。没配也能跑——设计解析回退到样本名
前缀，假设生成回退到基于效应量的规则引擎，只是自然语言那部分用不了。
文献检索（Europe PMC）与 LIPID MAPS 在线接口都不需要 key，且离线时静默降级不影响流程。

用法：
    python pb_agent_server.py            # 起服务并自动开浏览器
    python pb_agent_server.py --no-open --port 8790

注意：阶段 1/2 要读 .raw 并跑全库 C++ 搜索，耗时以分钟计；其余是秒级（8~11 取决于大模型）。
阶段之间有依赖，缺前置产物时会直接报出缺哪个文件，不会假装成功。
"""

import os
import re
import sys
import json
import time
import glob
import queue
import threading
import traceback
import webbrowser
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from flask import Flask, request, Response, jsonify

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))
PROJECT_ROOT = CODE_DIR.parent

import run as pb_run   # 复用 run.py 里已经调好的参数与阶段函数，避免参数在两处漂移
import pb_qc_diff
import pb_evidence as ev
import pb_settings as pbs
import pb_llm_analysis as pla
import pb_verify
import pb_report
from llm_client import LLMClient   # pla 已把主仓库 code/ 追加进 sys.path


def _load_logo_data_uri():
    """Read the LipidIN icon and return a base64 data URI (empty string if unavailable).

    读取 LipidIN 图标并返回 base64 data URI（读不到则返回空串）。PB_demo 下没有图标副本，
    回退到主仓库 code/icon.png，与主界面共用同一枚 logo。
    """
    import base64
    candidates = [CODE_DIR / "icon.png", PROJECT_ROOT.parent / "code" / "icon.png"]
    for icon_path in candidates:
        try:
            data = icon_path.read_bytes()
            return "data:image/png;base64," + base64.b64encode(data).decode("ascii")
        except Exception:
            continue
    return ""


LOGO_DATA_URI = _load_logo_data_uri()


STAGE_ORDER = ["annotation", "search", "localize", "quantify", "design", "qc", "plots",
               "interpret", "hypothesis", "verify", "report"]
STAGE_LABEL = {
    "annotation": "1 常规注释 (sample/common)",
    "search": "2 PB 衍生化检索 (sample/PB)",
    "localize": "3 双键位置重建",
    "quantify": "4 链×双键级定量 + 总量归一",
    "design": "5 实验设计（用哪些文件、跟谁比）",
    "qc": "6 质控 + 差异分析 + 位置偏好",
    "plots": "7 出图",
    "interpret": "8 大模型解读质控与差异",
    "hypothesis": "9 科学假设 + 文献 + LIPID MAPS",
    "verify": "10 子 agent 数据驱动验证",
    "report": "11 最终报告",
}
STAGE_LABEL_EN = {
    "annotation": "1 Annotation (sample/common)",
    "search": "2 PB derivatization search (sample/PB)",
    "localize": "3 Double-bond localization",
    "quantify": "4 Chain x C=C quantification + normalization",
    "design": "5 Design (which files, compared to what)",
    "qc": "6 QC + differential + position preference",
    "plots": "7 Figures",
    "interpret": "8 LLM review of QC and differential",
    "hypothesis": "9 Hypotheses + literature + LIPID MAPS",
    "verify": "10 Sub-agent data-driven verification",
    "report": "11 Final report",
}
# 下面这些路径全部**按当前配置动态算**，不能在 import 时定死：用户随时可能把数据目录
# 指到别处，定死了就会去老地方找文件、或者把报告写回老目录。
def report_dir():
    """报告落在 PB 产物旁边（figures/ 的同级），便于整包拷走。"""
    return os.path.join(str(pb_run.PB_OUT_DIR), "agent_report")


def report_md():
    return os.path.join(report_dir(), "PB_agent_report.md")


def stage_requires(action):
    """每个阶段的前置产物：缺了就直接说缺哪个，不要跑到一半炸在里面。"""
    C, P = str(pb_run.COMMON_OUT_DIR), str(pb_run.PB_OUT_DIR)
    table = {
        "search": [],
        "localize": [(pb_run.common_master_quant, "阶段1 的定量主表"),
                     (os.path.join(P, "*_PB_diagnostic.csv"), "阶段2 的诊断表")],
        "quantify": [(os.path.join(P, "*_PB_dbpos.csv"), "阶段3 的链×双键表")],
        "design": [(os.path.join(P, "PB_dbposition_quant.csv"), "阶段4 的定量表（要从中读样品名）")],
        "qc": [(os.path.join(P, "PB_dbposition_quant.csv"), "阶段4 的定量表"),
               (os.path.join(P, "PB_species_quant.csv"), "阶段4 的物质级定量表")],
        "plots": [(os.path.join(P, "PB_position_profile.csv"), "阶段6 的位置构成谱"),
                  (os.path.join(P, "PB_qc.csv"), "阶段6 的质控表")],
        "interpret": [(os.path.join(P, "PB_dbposition_diff.csv"), "阶段6 的差异表")],
        "hypothesis": [(os.path.join(P, "PB_dbposition_diff.csv"), "阶段6 的差异表")],
    }
    return table.get(action, [])


def stage_outputs(action):
    """每个阶段跑完应该落在盘上的东西，左栏用它显示"产物是否已存在"。"""
    C, P = str(pb_run.COMMON_OUT_DIR), str(pb_run.PB_OUT_DIR)
    table = {
        "annotation": [pb_run.common_master_quant],
        "search": [os.path.join(P, "*_PB_diagnostic.csv")],
        "localize": [os.path.join(P, "*_PB_dbpos.csv")],
        "quantify": [os.path.join(P, "PB_dbposition_quant.csv")],
        "design": [os.path.join(P, "PB_design.json")],
        "qc": [os.path.join(P, "PB_qc.csv")],
        "plots": [os.path.join(P, "figures", "*.svg")],
        "report": [report_md()],
    }
    return table.get(action, [])


class _Tee:
    """把阶段函数里的 print 同时送到控制台和 SSE 事件流。

    管线各模块都是直接 print 日志的，改造成回调注入成本高；这里在阶段执行期间临时
    接管 stdout，是最小侵入的做法。
    """

    def __init__(self, session, original):
        self.session = session
        self.original = original
        self._buf = ""

    def write(self, text):
        self.original.write(text)
        self._buf += text
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            if line.strip():
                self.session.log(line.rstrip())

    def flush(self):
        self.original.flush()


class PBAgentSession:
    """持有配置、事件订阅者和串行 worker。"""

    def __init__(self):
        self.lock = threading.Lock()
        self.subscribers = []
        self.busy = False
        self.current = None
        self.pending = None     # 多轮向导状态：llm_* / design_spec / hypothesis_confirm
        # 配置全部以 pb_settings.SETTINGS 为准，初值直接从 run.py 读——避免两处写同一个
        # 默认值然后慢慢漂移
        self.config = {
            "language": "zh",
            "llm": {"provider": "none", "api_key": None, "model": None, "base_url": None},
            "sample_context": None,
            "max_hypotheses": pla.MAX_HYPOTHESES,
            "lipidmaps_online": True,
            "alpha": ev.ALPHA,
            "min_n": pb_qc_diff.MIN_N,
        }
        for key, spec in pbs.SETTINGS.items():
            if spec["run"]:
                v = getattr(pb_run, spec["run"], None)
                self.config[key] = str(v) if isinstance(v, Path) else v
        # 跨阶段传递的东西：设计、表、通路上下文、解读、假设、验证
        self.cache = {}
        self.done = {k: False for k in STAGE_ORDER}
        self.jobs = queue.Queue()
        threading.Thread(target=self._worker_loop, daemon=True).start()
        # 起服务时先把上次存的设计读回来，接着上次干
        try:
            saved = pb_qc_diff.load_design(str(pb_run.PB_OUT_DIR))
            if saved:
                self.cache["design"] = saved
                self.done["design"] = True
        except Exception:
            pass

    # -------- language / i18n helpers -------- #
    @property
    def lang(self):
        return self.config.get("language", "zh")

    def t(self, zh, en=None):
        return zh if self.lang == "zh" or en is None else en

    def stage_label(self, key):
        return (STAGE_LABEL if self.lang == "zh" else STAGE_LABEL_EN).get(key, key)

    def llm(self):
        return LLMClient.from_dict(self.config.get("llm"))

    def tables(self, force=False):
        """PB 产物表。跑完 qc 会置脏，下次自动重读。"""
        if force or "tables" not in self.cache:
            self.cache["tables"] = ev.load_tables(str(pb_run.PB_OUT_DIR))
        return self.cache["tables"]

    def design(self):
        return self.cache.get("design")

    # ---------- 事件流 ----------
    def subscribe(self):
        q = queue.Queue()
        with self.lock:
            self.subscribers.append(q)
        q.put(self._state_event())
        return q

    def unsubscribe(self, q):
        with self.lock:
            if q in self.subscribers:
                self.subscribers.remove(q)

    def push(self, event):
        with self.lock:
            subs = list(self.subscribers)
        for q in subs:
            q.put(event)

    def assistant(self, zh, en=None):
        self.push({"type": "assistant", "text": self.t(zh, en)})

    def log(self, text):
        self.push({"type": "log", "text": str(text)})

    def progress(self, value, label=None):
        self.push({"type": "progress", "value": max(0, min(100, round(value or 0))),
                   **({"label": label} if label else {})})

    def set_busy(self, busy):
        self.busy = busy
        self.push({"type": "busy", "busy": busy})

    def _state_event(self):
        try:
            panel = _panel_data(self)
        except Exception:
            panel = None
        return {"type": "config", "busy": self.busy, "panel": panel}

    def emit_state(self):
        self.push(self._state_event())

    # ---------- 任务 ----------
    def enqueue(self, action):
        self.jobs.put(action)

    def _worker_loop(self):
        while True:
            action = self.jobs.get()
            actions = STAGE_ORDER if action == "full" else [action]
            for a in actions:
                if a not in STAGES:
                    self.assistant(f"未知阶段：{a}", f"Unknown stage: {a}")
                    continue
                missing = self._missing_inputs(a)
                if missing:
                    self.assistant("缺少前置产物，先跑对应阶段：\n" +
                                   "\n".join(f"  · {d}（{p}）" for p, d in missing),
                                   "Missing prerequisites; run the earlier stage first:\n" +
                                   "\n".join(f"  · {p}" for p, _d in missing))
                    break
                self.current = a
                self.set_busy(True)
                self.emit_state()
                start = time.time()
                # 每次跑之前把配置写回 run.py：用户可能刚在对话里改了路径或容差
                try:
                    pbs.apply_to_run(self.config, pb_run)
                except Exception as e:
                    self.assistant(f"配置写回失败：{type(e).__name__}: {e}",
                                   f"Failed to apply settings: {type(e).__name__}: {e}")
                    self.set_busy(False)
                    self.current = None
                    break
                self.assistant(f"▶ 开始 {STAGE_LABEL[a]}", f"▶ Starting {STAGE_LABEL_EN[a]}")
                original = sys.stdout
                sys.stdout = _Tee(self, original)
                try:
                    STAGES[a](self)
                    self.done[a] = True
                    dt = time.time() - start
                    self.assistant(f"✔ {STAGE_LABEL[a]} 完成，耗时 {dt:.1f}s",
                                   f"✔ {STAGE_LABEL_EN[a]} done in {dt:.1f}s")
                except Exception as e:
                    self.log(traceback.format_exc(limit=8))
                    self.assistant(f"✘ {STAGE_LABEL[a]} 失败：{type(e).__name__}: {e}",
                                   f"✘ {STAGE_LABEL_EN[a]} failed: {type(e).__name__}: {e}")
                    break
                finally:
                    sys.stdout = original
                    self.current = None
                    self.set_busy(False)
                    self.progress(0)
                    self.emit_state()
                # 阶段要用户回话（设计说明 / 假设确认）就停下等着——全流程不能替用户点头
                if self.pending and len(actions) > 1:
                    self.assistant("（全流程在这里停下等你回复，回复后会接着往下跑。）",
                                   "(The full run pauses here for your reply.)")
                    break

    def _missing_inputs(self, action):
        out = []
        for full, desc in stage_requires(action):
            if not glob.glob(full):
                out.append((full, desc))
        return out


def _has_outputs(s, key):
    """该阶段的产物是否已经在盘上——本次会话没跑过、但上次跑出来的也算。"""
    pats = stage_outputs(key)
    for pat in pats:
        if not glob.glob(pat):
            return False
    return bool(pats)


def _panel_data(s):
    """左侧界面面板的结构化数据：阶段状态、实验设计、假设、参数、目录。"""
    cfg = s.config
    stages = []
    for k in STAGE_ORDER:
        if s.current == k:
            state = "run"
        elif s.done[k]:
            state = "done"
        elif _has_outputs(s, k):
            state = "cached"
        else:
            state = "todo"
        stages.append({"key": k, "label": STAGE_LABEL[k], "label_en": STAGE_LABEL_EN[k],
                       "state": state})

    # 设计：把 {样本:组} 反过来变成 {组:[样本]} 给界面画
    design = s.design()
    groups = []
    if design:
        by_group = {}
        for smp, g in (design.get("groups") or {}).items():
            if g:
                by_group.setdefault(g, []).append(smp)
        for g in sorted(by_group):
            role = ("case" if g == design.get("case")
                    else "control" if g == design.get("control") else "other")
            groups.append({"name": g, "n": len(by_group[g]), "role": role,
                           "samples": sorted(by_group[g])})

    # 假设 + 验证结论（有验证就带上裁决）
    hyps = []
    vmap = {v["hypothesis"].get("id"): v for v in (s.cache.get("verifications") or [])}
    for h in (s.cache.get("hypotheses") or []):
        v = vmap.get(h.get("id"))
        nov = h.get("novelty") or {}
        hyps.append({
            "id": h.get("id"),
            "statement": h.get("statement"),
            "n_refs": nov.get("n_references", 0),
            "is_novel": bool(nov.get("is_novel")),
            "verdict": (v["data_verdict"]["verdict"] if v else None),
            "confidence": (v["confidence_score"] if v else None),
            "annotation": (v["annotation_confidence"]["level"] if v else None),
            "intensity": (v["intensity_confidence"]["level"] if v else None),
            "recommendation": (v["recommendation"] if v else None),
        })

    llm = s.config.get("llm") or {}
    return {
        "project_root": cfg.get("project_root"),
        "out_dir": cfg.get("pb_out_dir"),
        "paths": {k: cfg.get(k) for k, sp in pbs.SETTINGS.items() if sp["group"] == "path"},
        "report": report_md() if os.path.exists(report_md()) else None,
        "stages": stages,
        "design": ({"case": design.get("case"), "control": design.get("control"),
                    "groups": groups, "excluded": design.get("excluded") or [],
                    "note": design.get("note") or "", "source": design.get("source") or ""}
                   if design else None),
        "hypotheses": hyps,
        "context": cfg.get("sample_context") or "",
        "llm": {"provider": llm.get("provider", "none"), "model": llm.get("model") or "",
                "enabled": bool(s.llm().is_enabled())},
        "params": {k: cfg[k] for k in ("pb_ppm_err1", "pb_da_err1", "pb_da_err2",
                                       "rt_tol", "l12_min_coverage", "rt_cluster_tol")},
        "pending": s.pending,
        "current": s.current,
        "busy": s.busy,
    }


# ------------------------------------------------------------------ #
# 阶段：直接调 run.py 的函数，参数由 run.py 的模块级常量决定（单一事实来源）
# ------------------------------------------------------------------ #
def stage_annotation(s):
    s.progress(10, "annotation")
    pb_run.stage1_common_annotation()
    s.progress(100)


def stage_search(s):
    s.progress(10, "search")
    pb_run.stage2_pb_search()
    s.progress(100)


def stage_localize(s):
    s.progress(20, "localize")
    pb_run.stage3_pb_localize()
    s.progress(100)


def stage_quantify(s):
    s.progress(20, "quantify")
    pb_run.stage4_pb_quantify()
    s.progress(100)


def stage_design(s):
    """用哪些文件、跟谁比。有 pending 的自然语言就用它，否则问用户。"""
    s.progress(20, "design")
    samples = pb_qc_diff.available_samples(str(pb_run.PB_OUT_DIR))
    spec = s.cache.pop("design_spec", None)
    if not spec and s.cache.get("design") and s.cache.get("design_confirmed"):
        # 全流程里已经定过设计就别再问一遍
        s.assistant("沿用已确认的设计：\n" + pla.design_text(s.cache["design"], s.lang),
                    "Using the confirmed design:\n" + pla.design_text(s.cache["design"], s.lang))
        s.progress(100)
        return
    if not spec:
        # 没给说明就先按前缀给一个，并请用户用自然语言确认/修改
        design = pla.parse_design_spec(None, samples, "", s.lang)
        s.cache["design"] = design
        pb_qc_diff.save_design(str(pb_run.PB_OUT_DIR), design["groups"],
                               design["case"], design["control"], design.get("note", ""))
        s.assistant(
            "当前可用的样品：" + ", ".join(samples) + "\n\n" +
            pla.design_text(design, s.lang) +
            "\n\n这是按样本名前缀自动分的。要改就直接说，例如：\n"
            "  · 「只用 A1 到 A3 和 B1 到 B3，A 是病例 B 是对照」\n"
            "  · 「把 B119 排除掉」\n"
            "  · 「A 组叫 Tumor，B 组叫 Normal，Tumor 跟 Normal 比」",
            "Available samples: " + ", ".join(samples) + "\n\n" +
            pla.design_text(design, s.lang) +
            "\n\nThis is the sample-name-prefix default. Tell me in plain words to change it, e.g.\n"
            "  · 'use only A1-A3 and B1-B3, A is case, B is control'\n"
            "  · 'drop B119'")
        s.pending = "design_spec"
        s.progress(100)
        return
    design = pla.parse_design_spec(s.llm(), samples, spec, s.lang, log=s.log)
    s.cache["design"] = design
    s.cache["design_confirmed"] = True
    pb_qc_diff.save_design(str(pb_run.PB_OUT_DIR), design["groups"],
                           design["case"], design["control"], design.get("note", ""))
    src = ("大模型解析" if design.get("source") == "llm" else "按样本名前缀") if s.lang == "zh" \
        else ("parsed by LLM" if design.get("source") == "llm" else "by name prefix")
    s.assistant(f"（{src}）\n" + pla.design_text(design, s.lang) +
                ("\n\n设计已存入 PB_design.json。接着跑 `质控` 就会按这个设计做差异分析。"
                 if s.lang == "zh" else
                 "\n\nSaved to PB_design.json. Run `qc` to apply it."))
    s.progress(100)


def stage_qc(s):
    from pb_qc_diff import run_pb_qc_diff
    s.progress(30, "qc")
    design = s.design()
    if design:
        s.assistant(f"按设计 {design['case']} vs {design['control']} 做差异分析",
                    f"Differential per design: {design['case']} vs {design['control']}")
        run_pb_qc_diff(str(pb_run.PB_OUT_DIR), design=design.get("groups"),
                       case=design.get("case"), control=design.get("control"))
    else:
        run_pb_qc_diff(str(pb_run.PB_OUT_DIR))
    s.tables(force=True)      # 表变了，缓存作废
    s.progress(100)


def stage_plots(s):
    from pb_plots import run_pb_plots
    s.progress(40, "plots")
    run_pb_plots(str(pb_run.PB_OUT_DIR))
    fig_dir = os.path.join(str(pb_run.PB_OUT_DIR), "figures")
    s.assistant(f"图已出到 {fig_dir}（每张含 svg / pdf / png）",
                f"Figures written to {fig_dir} (svg / pdf / png each)")
    s.progress(100)


def stage_interpret(s):
    """大模型解读质控与差异分析。"""
    s.progress(30, "interpret")
    text = pla.interpret_results(s.llm(), s.tables(force=True), s.design(), s.lang, log=s.log)
    s.cache["interpretation"] = text
    s.cache["interpretation_is_llm"] = bool(s.llm().is_enabled())
    s.progress(90)
    s.assistant(("**质控与差异分析解读**\n\n" if s.lang == "zh"
                 else "**QC & differential review**\n\n") + text)
    if not s.llm().is_enabled():
        s.assistant("（这段是确定性摘要。回复 `配置大模型` 后可得到真正的解读。）",
                    "(Deterministic digest. Reply `configure model` for a real review.)")
    s.progress(100)


def stage_hypothesis(s):
    """科学假设 + 文献检索 + LIPID MAPS 联合分析，然后交给用户确认。"""
    s.progress(15, "hypothesis")
    tables = s.tables(force=True)
    s.assistant("正在做 LIPID MAPS 通路富集并连在线接口……",
                "Running LIPID MAPS pathway enrichment and querying the online API ...")
    ctx = pla.pathway_context(tables, report_dir(), s.lang,
                              online=s.config.get("lipidmaps_online", True), log=s.log)
    s.cache["pathway_ctx"] = ctx
    n_pb = len((ctx.get("pb") or {}).get("enrichment") or [])
    n_online = len(ctx.get("online") or {})
    s.assistant(f"通路富集：PB 层 {n_pb} 条通路；LIPID MAPS 在线注释 {n_online} 个类别。",
                f"Pathways: {n_pb} enriched on the PB layer; {n_online} subclasses annotated online.")
    s.progress(45)
    s.assistant("正在提出科学假设并逐条检索文献……",
                "Proposing hypotheses and searching the literature ...")
    hyps = pla.generate_hypotheses(s.llm(), tables, s.design(), ctx, s.lang,
                                   s.config.get("max_hypotheses", pla.MAX_HYPOTHESES), log=s.log)
    s.cache["hypotheses"] = hyps
    s.cache["verifications"] = []          # 假设变了，旧的验证作废
    s.done["verify"] = s.done["report"] = False
    s.progress(90)
    if not hyps:
        s.assistant("没能从这批数据里提出可支撑的假设。", "No supportable hypothesis from this data.")
        return
    s.assistant(("**科学假设（请确认）**\n\n" if s.lang == "zh" else "**Hypotheses (please confirm)**\n\n")
                + pla.hypotheses_text(hyps, s.lang) +
                ("\n\n回复 `确认` 开始逐条验证；或直接用自然语言说要怎么改，例如"
                 "「删掉第 2 条」「第 1 条聚焦到 PC 类」「再提一条关于 n-3 的」；"
                 "也可以 `重新生成`。" if s.lang == "zh" else
                 "\n\nReply `confirm` to verify, or just say how to change them "
                 "(e.g. 'drop 2', 'focus 1 on PC'), or `regenerate`."))
    s.pending = "hypothesis_confirm"
    s.progress(100)


def stage_verify(s):
    """每条假设一个子 agent，基于数据切片再验证。"""
    hyps = s.cache.get("hypotheses") or []
    if not hyps:
        s.assistant("还没有假设可验证，先跑 `假设`。", "No hypotheses yet; run `hypothesis` first.")
        return
    s.progress(20, "verify")
    s.assistant(f"启动 {len(hyps)} 个验证子 agent（每条假设一个，各自只看该假设名下的真实数据行，"
                f"并被要求优先反驳）……",
                f"Launching {len(hyps)} verification sub-agents (one per hypothesis; each sees only "
                f"its own data rows and is asked to refute) ...")
    vers = pb_verify.verify_hypotheses(hyps, s.tables(), s.llm(), s.design(), s.lang, log=s.log)
    s.cache["verifications"] = vers
    s.progress(90)
    s.assistant(("**验证结果**\n\n" if s.lang == "zh" else "**Verification**\n\n")
                + pb_verify.verification_text(vers, s.lang) +
                ("\n\n继续 `报告` 生成最终 md。" if s.lang == "zh" else "\n\nNext: `report`."))
    s.progress(100)


def stage_report(s):
    """最终 md：验证结果 + 三种可信度，并把关键结论凝练到对话窗口。"""
    vers = s.cache.get("verifications") or []
    if not vers:
        s.assistant("还没有验证结果，先跑 `验证`。", "No verifications yet; run `verify` first.")
        return
    s.progress(30, "report")
    design = s.design()
    params = {
        "pb_ppm_err1": s.config["pb_ppm_err1"], "pb_da_err1": s.config["pb_da_err1"],
        "pb_da_err2": s.config["pb_da_err2"], "rt_tol": s.config["rt_tol"],
        "l12_min_coverage": s.config["l12_min_coverage"],
        "rt_cluster_tol": s.config["rt_cluster_tol"],
        "case": (design or {}).get("case"), "control": (design or {}).get("control"),
        "llm_provider": s.config["llm"].get("provider"),
        "llm_model": s.config["llm"].get("model"),
        "alpha": ev.ALPHA, "min_n": pb_qc_diff.MIN_N,
        "out_dir": s.config.get("pb_out_dir"),
        "pb_raw_dir": s.config.get("pb_raw_dir"),
        "common_raw_dir": s.config.get("common_raw_dir"),
    }
    md = pb_report.build_pb_report(
        report_md(), s.tables(), design, s.cache.get("interpretation") or "",
        s.cache.get("pathway_ctx"), vers, s.lang, params,
        figures_dir=os.path.join(str(pb_run.PB_OUT_DIR), "figures"),
        sample_context=s.config.get("sample_context"),
        interpretation_is_llm=bool(s.cache.get("interpretation_is_llm")))
    s.progress(85)
    s.assistant(("**报告已生成**：" if s.lang == "zh" else "**Report written**: ") + md)
    s.assistant(pb_report.summarize_for_chat(vers, design, s.lang, s.tables()))
    s.progress(100)


STAGES = {
    "annotation": stage_annotation,
    "search": stage_search,
    "localize": stage_localize,
    "quantify": stage_quantify,
    "design": stage_design,
    "qc": stage_qc,
    "plots": stage_plots,
    "interpret": stage_interpret,
    "hypothesis": stage_hypothesis,
    "verify": stage_verify,
    "report": stage_report,
}


# ------------------------------------------------------------------ #
# 消息解析：确定性命令，不接 LLM
# ------------------------------------------------------------------ #
ALIASES = {
    "注释": "annotation", "annotation": "annotation", "annotate": "annotation", "1": "annotation",
    "搜库": "search", "检索": "search", "search": "search", "2": "search",
    "定位": "localize", "localize": "localize", "3": "localize",
    "定量": "quantify", "quantify": "quantify", "4": "quantify",
    "设计": "design", "分组": "design", "design": "design", "5": "design",
    "质控": "qc", "差异": "qc", "差异分析": "qc", "qc": "qc", "6": "qc",
    "出图": "plots", "画图": "plots", "plots": "plots", "figures": "plots", "7": "plots",
    "解读": "interpret", "interpret": "interpret", "8": "interpret",
    "假设": "hypothesis", "科学假设": "hypothesis", "hypothesis": "hypothesis",
    "hypotheses": "hypothesis", "9": "hypothesis",
    "验证": "verify", "verify": "verify", "10": "verify",
    "报告": "report", "report": "report", "11": "report",
    "全流程": "full", "运行全流程": "full", "full": "full", "all": "full",
    "full workflow": "full", "run all": "full",
}
CONFIRM_WORDS = ("确认", "确定", "接受", "confirm", "yes", "y", "ok")
REGEN_WORDS = ("重新生成", "再生成", "regenerate", "redo")


def interpret(s, text):
    raw = (text or "").strip()
    t = raw.lower()
    if not t:
        return

    # 多轮向导优先：大模型配置 / 设计说明 / 假设确认
    if s.pending:
        if s.pending.startswith("llm_"):
            _llm_wizard_step(s, raw)
        elif s.pending == "design_spec":
            _design_spec_step(s, raw)
        elif s.pending == "hypothesis_confirm":
            _hypothesis_confirm_step(s, raw)
        return

    if t in ("lang en", "english", "英文"):
        s.config["language"] = "en"
        s.assistant("语言已切换为英文。", "Switched to English.")
        s.emit_state()
        return
    if t in ("lang zh", "chinese", "中文"):
        s.config["language"] = "zh"
        s.assistant("已切换为中文。", "Switched to Chinese.")
        s.emit_state()
        return
    if re.search(r"配置大模型|配置模型|设置大模型|设置模型|接入大模型|配置\s*api|"
                 r"configure\s+(the\s+)?(llm|model)|llm\s*config", t):
        _start_llm_wizard(s)
        return
    # 研究背景：让假设知道这是什么样本，不然它只能泛泛而谈
    m = re.match(r"^(背景|context|研究背景)\s*[:：]\s*(.+)$", raw, re.I | re.S)
    if m:
        s.config["sample_context"] = m.group(2).strip()
        s.emit_state()
        s.assistant(f"研究背景已记下：{s.config['sample_context']}",
                    f"Study context noted: {s.config['sample_context']}")
        return
    if t in ("状态", "status", "config", "配置"):
        done = [s.stage_label(k) for k in STAGE_ORDER if s.done[k]]
        pend = [s.stage_label(k) for k in STAGE_ORDER if not s.done[k]]
        s.assistant("已完成：" + (", ".join(done) or "无") + "\n待运行：" + (", ".join(pend) or "无"),
                    "Done: " + (", ".join(done) or "none") + "\nPending: " + (", ".join(pend) or "none"))
        return
    if t in ("帮助", "help", "?"):
        _reply_help(s)
        return

    # 「设计 用A1到A3……」这种把说明写在同一句里的
    m = re.match(r"^(设计|分组|design)\s*[:：]?\s+(.+)$", raw, re.I | re.S)
    if m and len(m.group(2).strip()) > 2:
        s.cache["design_spec"] = m.group(2).strip()
        _enqueue_checked(s, "design")
        return

    if t in ("参数", "配置", "settings", "show config", "看参数"):
        s.assistant(pbs.describe(s.config, s.lang))
        return

    action = ALIASES.get(t)
    if action is not None:
        _enqueue_checked(s, action)
        return

    # 到这里说明不是固定命令 —— 交给自由对话解析（改参数 / 改路径 / 跑阶段）
    _free_form(s, raw)


# ------------------------------------------------------------------ #
# 自由对话：任何一句话都尽量落成"改哪个参数 / 跑哪个阶段"
# ------------------------------------------------------------------ #
def _apply_config(s, patch, note_prefix=""):
    """校验 -> 写回 run.py -> 提示哪些阶段的产物过期了。返回是否有改动。"""
    good, errors = pbs.validate_patch(patch)
    if errors:
        s.assistant((note_prefix + "这些没能改：\n" + "\n".join(f"  · {e}" for e in errors))
                    if s.lang == "zh" else
                    (note_prefix + "Could not apply:\n" + "\n".join(f"  · {e}" for e in errors)))
    if not good:
        return False
    s.config.update(good)
    try:
        changed = pbs.apply_to_run(s.config, pb_run)
    except Exception as e:
        s.assistant(f"写回 run.py 失败：{type(e).__name__}: {e}",
                    f"Failed to apply to run.py: {type(e).__name__}: {e}")
        return False
    # 分析层的参数不在 run.py 里，得各自送到用它们的模块
    if "alpha" in good:
        ev.ALPHA = pb_qc_diff.ALPHA = float(good["alpha"])
    if "min_n" in good:
        pb_qc_diff.MIN_N = int(good["min_n"])
    lines = [f"  · {k} = {v}" for k, v in good.items()]
    s.assistant(("已更新：\n" if s.lang == "zh" else "Updated:\n") + "\n".join(lines))
    if changed:
        s.log(f"[config] run.py 已更新：{changed}")
    stale = [x for x in pbs.stale_stages(good.keys()) if s.done.get(x) or _has_outputs(s, x)]
    if stale:
        names = "、".join(s.stage_label(x) for x in stale)
        s.assistant((f"注意：这些阶段的产物是用旧参数跑的，要重跑才作数 —— {names}")
                    if s.lang == "zh" else
                    (f"Note: these stages were run with the old settings and need a re-run — {names}"))
        for x in stale:
            s.done[x] = False
    s.emit_state()
    return True


def _free_form(s, raw):
    """没匹配上固定命令的话，走这里。

    先确定性地抠路径（模型转述路径会出错），再让大模型判断意图；没有大模型时给出
    能直接照抄的命令格式，而不是一句"没听懂"。
    """
    path = pla.grab_path(raw)
    intent = pla.parse_intent(s.llm(), raw, pbs.settings_catalog(s.lang),
                              list(STAGE_ORDER) + ["full", "status", "help", "show_config", "none"],
                              s.lang, log=s.log)

    if intent is None:                       # 没配大模型
        _free_form_offline(s, raw, path)
        return

    cfg = dict(intent.get("config") or {})
    # 路径以正则抠到的为准：模型很爱把 D:\x\y 规范化成 D:/x/y 或补个不存在的子目录
    if path:
        for k, v in list(cfg.items()):
            key = pbs.canonical(k)
            if key and pbs.SETTINGS.get(key, {}).get("kind") in ("dir", "dir_out", "file"):
                cfg[k] = path
    if intent.get("reply"):
        s.assistant(intent["reply"])
    applied = _apply_config(s, cfg) if cfg else False

    if intent.get("ask") and not applied:
        # 模型认出用户想改什么、但没给值 —— 问清楚，并给个能照抄的例子
        s.assistant(intent["ask"] + ("\n例如：`PB raw 目录 D:\\data\\PB`" if s.lang == "zh"
                                     else "\ne.g. `PB raw folder D:\\data\\PB`"))
        return
    act = intent.get("action") or "none"
    if act == "show_config":
        s.assistant(pbs.describe(s.config, s.lang))
        return
    if act == "help":
        _reply_help(s)
        return
    if act == "status":
        _reply_status(s)
        return
    if act in STAGES or act == "full":
        _enqueue_checked(s, act)
        return
    if not applied and not intent.get("reply"):
        s.assistant("我没把这句话对应到任何参数或阶段。输入 `帮助` 看我能做什么，"
                    "或 `参数` 看现在的配置。",
                    "I couldn't map that to a setting or a stage. Try `help` or `settings`.")


def _free_form_offline(s, raw, path):
    """没有大模型时的确定性兜底：认路径、认「键 = 值」，其余给出照抄格式。"""
    low = raw.lower()
    patch = {}
    # 「pb raw 目录 D:\...」这类：先看提到哪个路径项，再把抠到的路径给它
    if path:
        hint = None
        if re.search(r"\bpb\b.*(raw|原始|目录|地址|路径)|(raw|原始).*\bpb\b", low):
            hint = "pb_raw_dir"
        elif re.search(r"common|常规", low):
            hint = "common_raw_dir"
        elif re.search(r"(所有|全部|都).{0,6}(raw|原始)|(raw|原始).{0,6}(都|全部)|"
                       r"raw\s*(根|root)|原始(数据|文件)?(根目录|总目录|目录|文件夹|在)", low):
            hint = "raw_root"
        elif re.search(r"库|library", low):
            hint = "library_dir"
        elif re.search(r"输出|output|out", low):
            hint = "pb_out_dir"
        elif re.search(r"项目|根目录|数据(目录|地址|路径|文件夹)|data\s*root|project", low):
            hint = "project_root"
        if hint:
            patch[hint] = path
        else:
            s.assistant(f"我看到路径 {path}，但不确定要改哪一个。请照这样说：\n"
                        f"  · `raw 根目录 {path}` —— 其下要有 common/ 和 PB/ 两个子目录，"
                        f"我会自动分派并各自转 pkl（最常用）\n"
                        f"  · `项目根目录 {path}`\n  · `PB raw 目录 {path}`\n"
                        f"  · `common raw 目录 {path}`\n  · `谱库目录 {path}`\n"
                        f"  · `输出目录 {path}`\n（配置大模型后就能直接用自然语言说了。）",
                        f"I see the path {path} but not which setting it is. Say e.g. "
                        f"`raw root {path}` (with common/ and PB/ inside) / `project root {path}` / "
                        f"`PB raw folder {path}` / `library folder {path}`.")
            return
    # 「键 = 值」/「键 值」：键必须是表里认识的
    for key in pbs.SETTINGS:
        if pbs.SETTINGS[key]["kind"] in ("dir", "dir_out", "file"):
            continue
        m = re.search(re.escape(key) + r"\s*(?:=|＝|:|：|\s)\s*([-\d.]+|true|false|是|否)", low)
        if m:
            patch[key] = m.group(1)
    if patch:
        _apply_config(s, patch)
        return
    s.assistant(
        "我没配大模型，所以只认固定说法。你可以：\n"
        "  · `配置大模型` —— 配好之后这句话我就能直接听懂了\n"
        "  · 告诉我原始文件在哪：`raw 根目录 D:\\你的数据`（其下 common/ 与 PB/ 两个子目录）\n"
        "  · 改路径：`项目根目录 D:\\你的数据` / `PB raw 目录 D:\\...` / `谱库目录 D:\\...`\n"
        "  · 改参数：`pb_ppm_err1 = 40`、`alpha = 0.1`、`max_hypotheses = 3`\n"
        "  · `参数` 看当前配置，`帮助` 看全部命令",
        "No LLM configured, so I only understand fixed phrasings. Try:\n"
        "  · `configure model` — then I can understand this directly\n"
        "  · paths: `project root D:\\your\\data` / `PB raw folder D:\\...`\n"
        "  · params: `pb_ppm_err1 = 40`, `alpha = 0.1`\n"
        "  · `settings` for the current config, `help` for commands")


def _reply_status(s):
    done = [s.stage_label(k) for k in STAGE_ORDER if s.done[k]]
    pend = [s.stage_label(k) for k in STAGE_ORDER if not s.done[k]]
    s.assistant("已完成：" + (", ".join(done) or "无") + "\n待运行：" + (", ".join(pend) or "无"),
                "Done: " + (", ".join(done) or "none") + "\nPending: " + (", ".join(pend) or "none"))


def _enqueue_checked(s, action):
    if s.busy:
        s.assistant(f"正在跑 {s.stage_label(s.current)}，请等它跑完。",
                    f"{s.stage_label(s.current)} is running; please wait.")
        return
    s.enqueue(action)


def _reply_help(s):
    s.assistant(
        "我可以把 PB 双键定位从原始数据一路跑到带验证的报告。\n"
        "• 先告诉我原始文件在哪：`raw 根目录 D:\\...`（其下要有 common/ 与 PB/ 两个子目录），"
        "raw→pkl 由阶段 1/2 自动完成\n"
        "• 分步：" + " / ".join(f"`{k}`" for k in STAGE_ORDER) + "（也可输入序号 1~11）\n"
        "• `全流程` 依次跑完（跑到 `假设` 会停下来等你确认）\n"
        "• `配置大模型` —— 配好才能用自然语言说设计、才有真正的解读与假设\n"
        "• `设计` —— 说明用哪些文件、跟谁比，例如「只用 A1 到 A3 和 B1 到 B3，A 是病例」\n"
        "• `背景：肝癌 vs 癌旁` —— 告诉我这是什么样本，假设会更贴题\n"
        "• 假设列出来后：`确认` 开始验证，或直接说怎么改（「删掉第2条」「聚焦 PC 类」）\n"
        "• `状态` 看进度",
        "I run PB double-bond localization from raw data to a verified report.\n"
        "• Stages: " + " / ".join(f"`{k}`" for k in STAGE_ORDER) + " (or numbers 1-11)\n"
        "• `full` runs everything (pauses at `hypothesis` for your confirmation)\n"
        "• `configure model` — needed for natural-language design, real review and hypotheses\n"
        "• `design` — say which files to use and what to compare\n"
        "• `context: HCC vs adjacent` — tell me what the samples are\n"
        "• After hypotheses: `confirm` to verify, or say how to change them\n"
        "• `status` for progress")


# ---------- 设计说明 ----------
def _design_spec_step(s, raw):
    """用户对"用哪些文件、跟谁比"的自然语言回复。"""
    s.pending = None
    if raw.strip().lower() in CONFIRM_WORDS + ("默认", "default", "跳过", "skip"):
        s.assistant("好，沿用当前设计。接着可以跑 `质控`。",
                    "Keeping the current design. Next: `qc`.")
        s.emit_state()
        return
    if not s.llm().is_enabled():
        s.assistant("要用自然语言描述设计需要先配置大模型。回复 `配置大模型` 走一遍向导；"
                    "在那之前我只能按样本名前缀分组。",
                    "Natural-language design needs an LLM. Reply `configure model`; "
                    "until then I can only group by sample-name prefix.")
        s.emit_state()
        return
    s.cache["design_spec"] = raw
    _enqueue_checked(s, "design")


# ---------- 假设确认与修正 ----------
def _hypothesis_confirm_step(s, raw):
    """假设列出来之后的多轮：确认 / 删除 / 重新生成 / 自然语言修正。"""
    t = raw.strip().lower()
    if t in CONFIRM_WORDS:
        s.pending = None
        s.assistant("好，开始逐条验证。", "Confirmed. Starting verification.")
        s.enqueue("verify")
        s.enqueue("report")     # 验证完直接出报告，省一次来回
        return
    if t in REGEN_WORDS:
        s.pending = None
        s.enqueue("hypothesis")
        return
    # 「删除 2,4」是确定性的，不必惊动大模型
    m = re.search(r"(删除|删掉|去掉|去除|drop|remove)\s*([\d,，\s和and]+)", t)
    if m:
        ids = {int(x) for x in re.split(r"[^\d]+", m.group(2)) if x.strip().isdigit()}
        keep = [h for h in (s.cache.get("hypotheses") or []) if h.get("id") not in ids]
        if not keep:
            s.assistant("那样就一条都不剩了，没有执行。", "That would drop everything; not applied.")
            return
        for i, h in enumerate(keep, start=1):
            h["id"] = i
        s.cache["hypotheses"] = keep
        s.emit_state()
        s.assistant(f"已删除 {sorted(ids)}，剩 {len(keep)} 条：\n\n"
                    + pla.hypotheses_text(keep, s.lang) +
                    "\n\n回复 `确认` 开始验证，或继续说要怎么改。",
                    f"Dropped {sorted(ids)}, {len(keep)} left:\n\n"
                    + pla.hypotheses_text(keep, s.lang) +
                    "\n\nReply `confirm` to verify, or keep refining.")
        return
    # 其余交给大模型修正
    if not s.llm().is_enabled():
        s.assistant("要用自然语言修正假设需要先配置大模型。可以先用 `删除 2` 这种命令，"
                    "或 `确认` 直接验证。",
                    "Refining in natural language needs an LLM. Use `drop 2`, or `confirm`.")
        return
    s.assistant("正在按你的意思修正假设……", "Revising the hypotheses ...")
    new, err = pla.refine_hypotheses(s.llm(), s.cache.get("hypotheses") or [], raw,
                                     s.tables(), s.design(), s.cache.get("pathway_ctx"),
                                     s.lang, log=s.log)
    if err:
        s.assistant(err)
        return
    s.cache["hypotheses"] = new
    s.cache["verifications"] = []
    s.emit_state()
    s.assistant("**修正后的假设**\n\n" + pla.hypotheses_text(new, s.lang) +
                "\n\n回复 `确认` 开始验证，或继续说要怎么改。",
                "**Revised hypotheses**\n\n" + pla.hypotheses_text(new, s.lang) +
                "\n\nReply `confirm` to verify, or keep refining.")


# ---------- 大模型配置向导（与主界面同一套） ----------
def _start_llm_wizard(s):
    s.pending = "llm_provider"
    s.emit_state()
    s.assistant("好的，先配置大模型。请选择提供方：anthropic / openai / deepseek / custom"
                "（或回复 none 关闭）。",
                "Let's configure the LLM. Provider: anthropic / openai / deepseek / custom "
                "(or 'none' to disable).")


def _is_skip(text):
    return text.strip().lower() in ("默认", "default", "跳过", "skip", "-", "空")


def _llm_wizard_step(s, raw):
    step = s.pending
    if step == "llm_provider":
        provider = raw.strip().lower()
        if provider in ("none", "关闭", "取消", "cancel"):
            s.config["llm"] = {"provider": "none", "api_key": None, "model": None, "base_url": None}
            s.pending = None
            s.emit_state()
            s.assistant("已关闭大模型，将使用离线规则引擎。", "LLM disabled; offline rules only.")
            return
        if provider not in ("anthropic", "openai", "codex", "deepseek", "custom"):
            s.assistant("请从 anthropic / openai / deepseek / custom / none 中选。",
                        "Choose one of anthropic / openai / deepseek / custom / none.")
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
                    "Enter the model name (reply 'default' for the provider default).")
        return
    if step == "llm_model":
        s.config["llm"]["model"] = None if _is_skip(raw) else raw.strip()
        if s.config["llm"]["provider"] == "custom":
            s.pending = "llm_baseurl"
            s.assistant("请输入 Base URL（需兼容 OpenAI /chat/completions）。",
                        "Enter the Base URL (OpenAI /chat/completions compatible).")
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
        s.assistant(f"大模型已就绪（{llm.provider} · {llm.model}）。现在可以用自然语言说实验设计、"
                    f"也能得到真正的解读与科学假设。",
                    f"LLM ready ({llm.provider} · {llm.model}). Natural-language design, review "
                    f"and hypotheses are now available.")
    else:
        s.assistant("大模型未启用：缺少密钥或模型。可回复 `配置大模型` 重来，或继续用命令。",
                    "LLM not enabled: missing key or model. Reply `configure model` to retry.")


SESSION = PBAgentSession()
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
    text = (request.get_json(silent=True) or {}).get("text", "")
    if text.strip():
        threading.Thread(target=interpret, args=(SESSION, text), daemon=True).start()
    return jsonify({"ok": True})


@app.route("/api/state")
def api_state():
    return jsonify(SESSION._state_event())


@app.route("/api/stream")
def api_stream():
    def gen():
        q = SESSION.subscribe()
        try:
            yield "retry: 3000\n\n"
            while True:
                try:
                    ev = q.get(timeout=15)
                    yield f"data: {json.dumps(ev, ensure_ascii=False)}\n\n"
                except queue.Empty:
                    yield ": keepalive\n\n"
        finally:
            SESSION.unsubscribe(q)
    return Response(gen(), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no",
                             "Connection": "keep-alive"})


INDEX_HTML = r"""<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<link rel="icon" href="__FAVICON__"/>
<title>LipidIN PB Agent</title>
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
.grp{display:flex;align-items:center;gap:7px;padding:6px 8px;border-radius:9px;font-size:12px}
.grp:nth-child(odd){background:rgba(128,128,140,.06)}
.grp .n{margin-left:auto;color:var(--muted);font-variant-numeric:tabular-nums}
.role{font-size:10px;padding:1px 7px;border-radius:999px;font-weight:650;flex:none}
.role.run{background:rgba(109,94,252,.14);color:var(--accent)}
.role.done{background:rgba(24,194,156,.14);color:var(--ok)}
.role.cached{background:rgba(232,161,58,.14);color:var(--warn)}
.role.todo{background:rgba(128,128,140,.12);color:var(--muted)}
.stagerow{display:flex;align-items:center;gap:7px;padding:6px 8px;border-radius:9px;font-size:12px}
.stagerow:nth-child(odd){background:rgba(128,128,140,.06)}
.stagerow.active{background:rgba(109,94,252,.10);border:1.5px solid rgba(109,94,252,.22)}
.stagerow .nm{word-break:break-word}
.pip{width:8px;height:8px;border-radius:50%;background:var(--line);flex:none}
.pip.done{background:var(--ok)} .pip.cached{background:var(--warn)}
.pip.run{background:var(--accent);animation:pulse 1.2s infinite}
.note{font-size:11px;color:var(--muted);line-height:1.5;margin:6px 2px}
.mini{font-size:11px;color:var(--accent);cursor:pointer;border:1.5px solid var(--line);border-radius:999px;
  padding:3px 9px;background:var(--panel);transition:.15s;margin-top:6px;display:inline-block}
.mini:hover{border-color:var(--accent)}
.comp{display:flex;align-items:center;gap:6px;padding:5px 9px;margin:4px 0;border-radius:9px;font-size:12px;
  background:rgba(109,94,252,.08);border:1.5px solid rgba(109,94,252,.18)}
.comp b{color:var(--accent)}
.role.cc{background:rgba(109,94,252,.14);color:var(--accent)}
.role.group{background:rgba(128,128,140,.12);color:var(--muted)}
.bioq{font-size:12px;line-height:1.5;color:var(--ink);background:rgba(128,128,140,.06);
  border:1.5px dashed var(--line);border-radius:11px;padding:8px 10px;white-space:pre-wrap}
.bioq.empty{color:var(--muted)}
.hyp{border:1.5px solid var(--line);border-radius:11px;padding:7px 9px;margin:6px 0}
.hyp-h{display:flex;align-items:center;gap:6px}
.hyp-t{font-size:11.5px;line-height:1.5;color:var(--ink);margin-top:4px}
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
      <div class="subline"><span class="subtitle" id="sub">PB 双键定位 · 对话式 · 六阶段流程</span><span class="dot" id="dot"></span></div>
    </div>
    <div class="spacer"></div>
    <div class="toggle" id="langbtn">中 / EN</div>
  </div>
  <div class="cols">
    <!-- LEFT: stages / parameters / outputs -->
    <div class="col">
      <div class="col-title"><span id="lt-param">运行状态</span> <span class="en" id="lt-param-en">Live status</span></div>
      <div class="col-body" id="leftbody">
        <div class="note" id="left-empty">尚未连接。<br/>Not connected yet.</div>
      </div>
    </div>
    <!-- MIDDLE: chat -->
    <div class="col mid">
      <div class="progress"><div class="bar idle" id="bar"></div></div>
      <div class="chat" id="chat"></div>
      <div class="suggest" id="suggest"></div>
      <div class="inputbar">
        <div class="wrap"><textarea id="box" rows="1" placeholder="输入指令，如：全流程 / 定位 / 4 / 状态 / 帮助"></textarea></div>
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
const SUGG={zh:['配置大模型','全流程','设计','质控','解读','假设','确认','报告','状态','帮助'],
             en:['configure model','full','design','qc','interpret','hypothesis','confirm','report','status','help']};
const L={zh:{param:'运行状态',cite:'引用格式',st:'流程阶段',rp:'运行参数',io:'目录',
    dsg:'实验设计',hyp:'科学假设',llm:'大模型',bq:'研究背景',
    root:'项目根目录',out:'输出目录',rpt:'报告',
    ppm1:'PB MS1 ppm',da1:'PB MS1 Da',da2:'PB MS2 Da',rt:'RT 容差',cov:'L1/L2 最低覆盖',clu:'RT 聚类容差',
    stateNote:'阶段 1/2 要读 .raw 并全库搜索，耗时以分钟计；其余是秒级（8~11 取决于大模型）。',
    depNote:'阶段之间有依赖；缺前置产物时会直接告诉你缺哪个文件。',
    runFull:'▶ 运行全流程',
    noDsg:'尚未设置。点这里用一句话说明用哪些文件、跟谁比。',
    setDsg:'✎ 设置实验设计',
    noHyp:'尚未生成。跑完 `质控` 后点 `假设`。',
    bqEmpty:'未设置。可输入：背景：肝癌 vs 癌旁',
    editBq:'✎ 设置研究背景',
    llmOff:'未配置（自然语言与假设解读不可用）',
    lgDone:'本次已跑',lgCached:'盘上已有',lgRun:'运行中',lgTodo:'待运行',
    sub:'PB 双键定位 · 对话式 · 假设生成与子 agent 验证',
    ph:'输入指令，如：全流程 / 设计 / 假设 / 确认；或直接说要怎么改',
    copy:'复制',copied:'已复制',openArt:'查看文章'},
  en:{param:'Live status',cite:'How to cite',st:'Stages',rp:'Run parameters',io:'Folders',
    dsg:'Design',hyp:'Hypotheses',llm:'LLM',bq:'Study context',
    root:'Project root',out:'Output folder',rpt:'Report',
    ppm1:'PB MS1 ppm',da1:'PB MS1 Da',da2:'PB MS2 Da',rt:'RT tolerance',cov:'L1/L2 min coverage',clu:'RT cluster tol.',
    stateNote:'Stages 1-2 read .raw and search the full library (minutes); the rest are seconds.',
    depNote:'Stages depend on each other; a missing prerequisite is reported by file name.',
    runFull:'▶ Run full workflow',
    noDsg:'Not set. Click to say which files to use and what to compare.',
    setDsg:'✎ Set the design',
    noHyp:'None yet. Run `qc`, then `hypothesis`.',
    bqEmpty:'Not set. Type: context: HCC vs adjacent',
    editBq:'✎ Set study context',
    llmOff:'Not configured (no natural language / real review)',
    lgDone:'done',lgCached:'on disk',lgRun:'running',lgTodo:'pending',
    sub:'PB double-bond localization · conversational · hypotheses + sub-agent verification',
    ph:'Type a command: full / design / hypothesis / confirm — or just say what to change',
    copy:'Copy',copied:'Copied',openArt:'Open article'}};
const VERDICT={zh:{supported:'数据支持',partly_supported:'部分支持',not_supported:'数据不支持',
    insufficient_evidence:'证据不足',error:'验证失败'},
  en:{supported:'supported',partly_supported:'partly',not_supported:'not supported',
    insufficient_evidence:'insufficient',error:'error'}};

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
  let h='';
  // Stages
  h+='<div class="sec"><h4>'+T.st+'</h4>';
  (p.stages||[]).forEach(s=>{ const nm=lang==='zh'?s.label:s.label_en;
    const tag={done:T.lgDone,cached:T.lgCached,run:T.lgRun,todo:T.lgTodo}[s.state];
    h+=`<div class="stagerow${s.state==='run'?' active':''}" onclick="hint('${esc(s.key)}')">`+
       `<span class="pip ${s.state}"></span><span class="nm">${esc(nm)}</span>`+
       `<span class="role ${s.state}" style="margin-left:auto">${esc(tag)}</span></div>`; });
  h+='<div class="mini" onclick="sendMsg(\''+(lang==='zh'?'全流程':'full')+'\')">'+T.runFull+'</div>';
  h+='<div class="note">'+T.stateNote+'</div>';
  h+='<div class="note">'+T.depNote+'</div>';
  h+='</div>';
  // Design
  h+='<div class="sec"><h4>'+T.dsg+'</h4>';
  if(p.design){
    h+=`<div class="comp">⚖ <b>${esc(p.design.case)}</b> vs ${esc(p.design.control)}</div>`;
    (p.design.groups||[]).forEach(g=>{
      const rc=g.role==='case'||g.role==='control'?'cc':'group';
      h+=`<div class="grp"><span class="role ${rc}">${esc(g.role)}</span><b>${esc(g.name)}</b>`+
         `<span class="n">×${g.n}</span></div>`;
      h+=`<div class="note" style="margin:0 2px 4px">${esc((g.samples||[]).join(', '))}</div>`; });
    if((p.design.excluded||[]).length)
      h+='<div class="note">'+(lang==='zh'?'排除：':'Excluded: ')+esc(p.design.excluded.join(', '))+'</div>';
    if(p.design.note) h+='<div class="note">'+esc(p.design.note)+'</div>';
  } else { h+='<div class="note">'+T.noDsg+'</div>'; }
  h+='<div class="mini" onclick="sendMsg(\''+(lang==='zh'?'设计':'design')+'\')">'+T.setDsg+'</div>';
  h+='</div>';
  // Hypotheses
  h+='<div class="sec"><h4>'+T.hyp+'</h4>';
  if((p.hypotheses||[]).length){
    p.hypotheses.forEach(x=>{
      const vd=x.verdict?(VERDICT[lang][x.verdict]||x.verdict):null;
      const cls=x.verdict==='supported'?'done':(x.verdict==='not_supported'||x.verdict==='error'?'cached':
                 (x.verdict?'run':'todo'));
      h+='<div class="hyp">';
      h+=`<div class="hyp-h"><span class="role ${cls}">${vd?esc(vd):'#'+x.id}</span>`+
         (x.confidence!=null?`<span class="n">${(lang==='zh'?'可信度 ':'conf ')}${x.confidence.toFixed(2)}</span>`:'')+
         '</div>';
      h+=`<div class="hyp-t">${esc(x.statement)}</div>`;
      let meta=[];
      if(x.annotation) meta.push((lang==='zh'?'注释 ':'annot ')+esc(x.annotation));
      if(x.intensity) meta.push((lang==='zh'?'强度 ':'int ')+esc(x.intensity));
      meta.push((lang==='zh'?'文献 ':'lit ')+x.n_refs+(x.is_novel?(lang==='zh'?'（新）':' (novel)'):''));
      h+=`<div class="note" style="margin:3px 2px 0">${meta.join(' ｜ ')}</div>`;
      if(x.recommendation) h+=`<div class="note" style="margin:1px 2px 0;color:var(--accent)">${esc(x.recommendation)}</div>`;
      h+='</div>'; });
  } else { h+='<div class="note">'+T.noHyp+'</div>'; }
  h+='</div>';
  // Study context
  h+='<div class="sec"><h4>'+T.bq+'</h4>';
  h+= p.context?('<div class="bioq">'+esc(p.context)+'</div>')
                :('<div class="bioq empty">'+T.bqEmpty+'</div>');
  h+='<div class="mini" onclick="hint(\''+(lang==='zh'?'背景：':'context: ')+'\')">'+T.editBq+'</div>';
  h+='</div>';
  // LLM
  h+='<div class="sec"><h4>'+T.llm+'</h4>';
  if(p.llm && p.llm.enabled) h+=kv(esc(p.llm.provider), esc(p.llm.model||''));
  else { h+='<div class="note">'+T.llmOff+'</div>';
    h+='<div class="mini" onclick="sendMsg(\''+(lang==='zh'?'配置大模型':'configure model')+'\')">'+
       (lang==='zh'?'✎ 配置大模型':'✎ Configure')+'</div>'; }
  h+='</div>';
  // Run parameters
  const q=p.params||{};
  h+='<div class="sec"><h4>'+T.rp+'</h4>';
  h+=kv(T.ppm1, esc(q.pb_ppm_err1));
  h+=kv(T.da1, esc(q.pb_da_err1));
  h+=kv(T.da2, esc(q.pb_da_err2));
  h+=kv(T.rt, esc(q.rt_tol));
  h+=kv(T.cov, esc(q.l12_min_coverage));
  h+=kv(T.clu, esc(q.rt_cluster_tol));
  h+='</div>';
  // Folders
  h+='<div class="sec"><h4>'+T.io+'</h4>';
  h+=kv(T.root, esc(p.project_root));
  h+=kv(T.out, esc(p.out_dir));
  if(p.report) h+=kv(T.rpt, esc(p.report));
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
  else if(m.type==='config'){ lastPanel=m.panel; renderPanel(m.panel); setBusy(!!m.busy); }
};
function sendMsg(t){ t=(t||box.value).trim(); if(!t)return; add('user',t); box.value=''; auto();
  fetch('/api/message',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({text:t})}); }
window.sendMsg=sendMsg;
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
add('assistant', '你好，我是 PB 双键定位助手。左栏实时显示流程阶段、实验设计、科学假设与验证结论；右栏是多期刊引用格式。\n\n'
  +'完整链路：注释 → 搜库 → 定位 → 定量 → 设计 → 质控 → 出图 → 解读 → 假设(含文献+LIPID MAPS) → 子agent验证 → 报告。\n\n'
  +'建议先“配置大模型”，然后用一句话告诉我用哪些文件、跟谁比（点左栏“设置实验设计”）。\n'
  +'假设生成后我会先给你确认，你可以直接说“删掉第2条”“聚焦PC类”，改到满意再“确认”验证。\n'
  +'没配大模型也能跑，只是设计按样本名前缀、假设按效应量规则生成。\n\n'
  +'Hi, I am the PB double-bond localization assistant. Left: stages, design, hypotheses and verdicts. Right: citation formats.');
</script>
</body>
</html>
"""


def main():
    import argparse
    ap = argparse.ArgumentParser(description="PB 双键定位工作流 Agent")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8790)
    ap.add_argument("--no-open", action="store_true", help="不自动打开浏览器")
    args = ap.parse_args()
    url = f"http://{args.host}:{args.port}/"
    print(f"[INFO] LipidIN PB Agent UI at {url}\n[信息] PB 双键定位智能体界面： {url}")
    if not args.no_open:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host=args.host, port=args.port, threaded=True, debug=False)


if __name__ == "__main__":
    main()
