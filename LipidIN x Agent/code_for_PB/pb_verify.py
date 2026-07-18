# -*- coding: utf-8 -*-
"""
假设验证：每条科学假设派一个子 agent，拿**真实数据切片**再验一遍。

为什么要再验一遍：提假设的那次大模型调用看到的是全局摘要，它有动力把故事讲圆。验证子
agent 只看**这条假设名下的具体数据行**（份额、q 值、象限、n），并且被明确要求"优先反驳"。
两次调用的立场是对立的，所以能筛掉"讲得通但数据不支持"的假设——这正是提假设的模型
自己做不到的事。

每条假设 1 个子 agent，并发跑（IO 密集，线程池够用）。子 agent 之间互不可见，避免相互
传染。子 agent 挂了不影响别人，那条给"验证失败"而不是假装通过。

最终每条假设给四个独立的数：
  data_verdict      —— 子 agent 基于数据的裁决（支持/部分/不支持/证据不足）
  annotation_conf   —— 注释可信度（pb_evidence 算，确定性）
  intensity_conf    —— 强度可信度（pb_evidence 算，确定性）
  novelty           —— 是否够新（文献命中数）
四个不合并成一个分数：一条"数据强、注释好、但早有人报道"的结论和一条"数据弱、但没人报过"
的结论，价值完全不同，合并成一个数就把这个区别抹掉了。只在最后给一个**推荐等级**做排序。
"""

import concurrent.futures as futures

import pb_evidence as ev

VERIFY_SYSTEM = """You are an adversarial reviewer. Your job is to REFUTE a hypothesis drawn
from a Paternò–Büchi (PB) double-bond localization experiment, not to confirm it.

You are shown ONLY the real data rows that belong to this hypothesis, plus the literature that
was retrieved for it. Judge whether those rows actually support the claim.

Return ONLY JSON:
{"verdict": "supported" | "partly_supported" | "not_supported" | "insufficient_evidence",
 "score": <0.0-1.0>,
 "comment": "<2-4 sentences in {lang}, citing the actual numbers you were shown>",
 "contradictions": ["<any way the data fails to support the claim>"],
 "what_would_settle_it": "<the concrete experiment/analysis that would decide this>"}

Hard rules:
- Default to skepticism. If the rows do not clearly support the claim, say not_supported or
  insufficient_evidence. A hypothesis that merely "sounds plausible" is NOT supported.
- If nothing reached FDR<0.05, the strongest verdict you may give is "partly_supported", and
  your comment must say the evidence is effect-size only at this sample size.
- If the data slice shows 0 matching rows, the verdict MUST be "insufficient_evidence" with
  score <= 0.1.
- Check the direction: if the claim says a position is favored in one group but the numbers
  show the opposite, that is a contradiction — report it.
- Small n (n<3 per group -> the row was never tested) is itself a limitation worth stating.
- score meaning: 1.0 = data clearly supports; 0.5 = effect present but underpowered;
  0.0 = data contradicts the claim."""

VERDICT_ZH = {
    "supported": "数据支持",
    "partly_supported": "部分支持（效应量层面）",
    "not_supported": "数据不支持",
    "insufficient_evidence": "证据不足",
    "error": "验证失败",
}


def _log(log, msg):
    (log or print)(msg)


def verify_one(hypothesis, tables, llm, design_json=None, lang="zh", log=None):
    """单条假设的验证子 agent：数据切片 + 文献 -> 裁决 + 三种可信度。"""
    slice_text = ev.evidence_slice(tables, hypothesis, design_json)
    ann = ev.annotation_confidence(tables, hypothesis, lang)
    inten = ev.intensity_confidence(tables, hypothesis, lang)
    nov = hypothesis.get("novelty") or ev.novelty_score(hypothesis.get("references"), lang)

    verdict = _judge(hypothesis, slice_text, ann, llm, lang, log)
    conf = round(0.40 * verdict["score"] + 0.35 * ann["score"] + 0.25 * inten["score"], 3)
    return {
        "hypothesis": hypothesis,
        "data_slice": slice_text,
        "data_verdict": verdict,
        "annotation_confidence": ann,
        "intensity_confidence": inten,
        "novelty": nov,
        "confidence_score": conf,
        "recommendation": _recommend(conf, verdict, nov, lang),
    }


def _judge(hypothesis, slice_text, ann, llm, lang="zh", log=None):
    """子 agent 的裁决。没有大模型时用确定性规则，绝不假装通过。"""
    if ann.get("n", 0) == 0:
        return {"verdict": "insufficient_evidence", "score": 0.0,
                "comment": ("这条假设在数据里匹配不到任何行。" if lang == "zh"
                            else "No data rows match this hypothesis."),
                "contradictions": ["0 matching rows"], "what_would_settle_it": "",
                "source": "rule"}
    if llm is None or not llm.is_enabled():
        return _judge_rule(hypothesis, slice_text, lang)
    refs = hypothesis.get("references") or []
    lit = "\n".join(f"- {r.get('title')} ({r.get('year')}, {r.get('journal')}), "
                    f"cited {r.get('cited_by')}" for r in refs[:5]) or "（未检索到文献）"
    system = VERIFY_SYSTEM.replace("{lang}", "Chinese" if lang == "zh" else "English")
    user = (f"Hypothesis:\n{hypothesis.get('statement')}\n\n"
            f"Its stated rationale:\n{hypothesis.get('rationale')}\n\n"
            f"Proposed mechanism:\n{hypothesis.get('biological_mechanism') or '(none given)'}\n\n"
            f"=== THE ONLY DATA YOU MAY USE ===\n{slice_text}\n\n"
            f"=== Literature retrieved for this hypothesis ===\n{lit}\n\n"
            f"Remember: your job is to try to refute this.")
    try:
        data = llm.chat_json(system, user)
        v = str(data.get("verdict") or "insufficient_evidence")
        if v not in VERDICT_ZH:
            v = "insufficient_evidence"
        try:
            score = float(data.get("score"))
        except (TypeError, ValueError):
            score = 0.0
        return {"verdict": v, "score": max(0.0, min(1.0, score)),
                "comment": str(data.get("comment") or ""),
                "contradictions": [str(x) for x in (data.get("contradictions") or [])],
                "what_would_settle_it": str(data.get("what_would_settle_it") or ""),
                "source": "llm"}
    except Exception as e:
        _log(log, f"[verify] 子 agent 失败：{type(e).__name__}: {e}")
        return {"verdict": "error", "score": 0.0, "comment": f"{type(e).__name__}: {e}",
                "contradictions": [], "what_would_settle_it": "", "source": "error"}


def _judge_rule(hypothesis, slice_text, lang="zh"):
    """无大模型的回退裁决：只认数据里数得出来的东西。"""
    has_sig = "share_q=" in slice_text and any(
        f"share_q={v}" in slice_text for v in ("0.0", "1e-", "2e-", "3e-", "4e-"))
    has_flip = "翻转=True" in slice_text or "preference_flip" in slice_text
    if has_sig:
        return {"verdict": "partly_supported", "score": 0.5,
                "comment": ("匹配到数据行且存在显著项，但未经大模型独立复核。" if lang == "zh"
                            else "Matching rows with significance, but no independent review."),
                "contradictions": [], "what_would_settle_it": "", "source": "rule"}
    if has_flip:
        return {"verdict": "partly_supported", "score": 0.4,
                "comment": ("存在位置偏好翻转（效应量层面），未达统计显著。" if lang == "zh"
                            else "Position-preference flip present; effect size only."),
                "contradictions": [], "what_would_settle_it": "", "source": "rule"}
    return {"verdict": "insufficient_evidence", "score": 0.2,
            "comment": ("匹配到数据行，但没有达到效应量或显著性门槛。" if lang == "zh"
                        else "Rows matched but no effect size or significance threshold met."),
            "contradictions": [], "what_would_settle_it": "", "source": "rule"}


def _recommend(conf, verdict, nov, lang="zh"):
    """推荐等级：可信度 × 新颖度。两个轴都高才值得优先跟进。"""
    novel = bool(nov.get("is_novel"))
    if verdict["verdict"] in ("not_supported", "error") or conf < 0.3:
        return "不建议跟进" if lang == "zh" else "do not pursue"
    if conf >= 0.6 and novel:
        return "优先跟进（可信且少有人报道）" if lang == "zh" else "high priority (solid + novel)"
    if conf >= 0.6:
        return "可跟进（可信但已有报道）" if lang == "zh" else "worth pursuing (solid but reported)"
    if novel:
        return "谨慎跟进（新颖但证据弱，需扩样本）" if lang == "zh" else "cautious (novel but weak)"
    return "低优先级" if lang == "zh" else "low priority"


def verify_hypotheses(hypotheses, tables, llm, design_json=None, lang="zh", log=None,
                      max_workers=None):
    """给每条假设派一个子 agent，并发验证。返回顺序与输入一致。"""
    if not hypotheses:
        return []
    n = len(hypotheses)
    workers = max_workers or min(n, 5)
    _log(log, f"[verify] 启动 {n} 个验证子 agent（每条假设一个，并发 {workers}）……")
    results = [None] * n
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        fut = {pool.submit(verify_one, h, tables, llm, design_json, lang, log): i
               for i, h in enumerate(hypotheses)}
        for f in futures.as_completed(fut):
            i = fut[f]
            try:
                results[i] = f.result()
                v = results[i]["data_verdict"]["verdict"]
                _log(log, f"[verify] 假设 {hypotheses[i].get('id', i+1)} 裁决："
                          f"{VERDICT_ZH.get(v, v)}（可信度 {results[i]['confidence_score']}）")
            except Exception as e:
                _log(log, f"[verify] 假设 {i+1} 子 agent 异常：{type(e).__name__}: {e}")
                results[i] = {
                    "hypothesis": hypotheses[i], "data_slice": "",
                    "data_verdict": {"verdict": "error", "score": 0.0, "comment": str(e),
                                     "contradictions": [], "what_would_settle_it": "",
                                     "source": "error"},
                    "annotation_confidence": {"n": 0, "score": 0.0, "level": "低", "detail": ""},
                    "intensity_confidence": {"n": 0, "score": 0.0, "level": "低", "detail": ""},
                    "novelty": hypotheses[i].get("novelty") or ev.novelty_score([], lang),
                    "confidence_score": 0.0,
                    "recommendation": "不建议跟进" if lang == "zh" else "do not pursue",
                }
    return results


def verification_text(verifications, lang="zh"):
    """验证结果凝练成对话窗口里的一段话——只留能决定下一步的信息。"""
    if not verifications:
        return "（没有验证结果）" if lang == "zh" else "(no verifications)"
    lines = []
    for v in verifications:
        h = v["hypothesis"]
        dv = v["data_verdict"]
        nov = v["novelty"]
        vd = VERDICT_ZH.get(dv["verdict"], dv["verdict"]) if lang == "zh" else dv["verdict"]
        lines.append(f"[{h.get('id')}] {vd}｜可信度 {v['confidence_score']:.2f}"
                     f"（注释 {v['annotation_confidence']['level']}/强度 "
                     f"{v['intensity_confidence']['level']}）｜文献 {nov.get('n_references', 0)} 篇"
                     f"｜{v['recommendation']}")
        lines.append(f"     {h.get('statement')}")
        if dv.get("comment"):
            lines.append(f"     子 agent：{dv['comment']}" if lang == "zh"
                         else f"     Sub-agent: {dv['comment']}")
        for c in (dv.get("contradictions") or [])[:2]:
            lines.append(f"     ⚠ {c}")
    return "\n".join(lines)
