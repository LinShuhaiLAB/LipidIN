# -*- coding: utf-8 -*-
"""
最终报告：把设计、质控、两层差异、通路、假设与验证写成一份 md。

形态参考主仓库的 `report_writer.py`（同样是"参数可复现 + 图 + 逐条结论 + 附录"的结构），
但内容是 PB 特有的：两层差异、双键位置偏好、以及每条结论的四个独立评分。

核心原则：**不把四个分数合并成一个总分**。
  可信度（数据裁决 + 注释 + 强度）回答"这事是不是真的"；
  新颖度回答"值不值得写"。
两者正交，合并会把"扎实但老生常谈"和"新鲜但证据薄"混成同一个数，而这两条的下一步动作
完全不同。所以报告里并排列出，另给一个推荐等级用于排序。
"""

import os
import datetime

import pb_evidence as ev
from pb_verify import VERDICT_ZH


def _img(md_dir, path, caption):
    """图片存在才插入，路径相对 md 自身。"""
    if not path or not os.path.exists(path):
        return []
    try:
        rel = os.path.relpath(path, md_dir).replace("\\", "/")
    except ValueError:
        rel = path.replace("\\", "/")
    return [f"![{caption}]({rel})", f"*{caption}*", ""]


def _score_bar(score):
    """0-1 分数画成 5 格方块，扫一眼就知道强弱。"""
    filled = int(round(max(0.0, min(1.0, score)) * 5))
    return "█" * filled + "░" * (5 - filled)


def _demote_headings(text):
    """把插入文本里的 markdown 标题降级成粗体。

    大模型偶尔会自作主张加 ### 标题，直接插进来会把报告的章节结构冲掉（目录里冒出
    莫名其妙的层级）。降级成粗体即可，既保留强调又不参与章节树。
    """
    out = []
    for line in (text or "").split("\n"):
        m = line.lstrip()
        if m.startswith("#"):
            stripped = m.lstrip("#").strip()
            out.append(f"**{stripped}**" if stripped else "")
        else:
            out.append(line)
    return "\n".join(out)


def build_pb_report(md_path, tables, design_json, interpretation, pathway_ctx,
                    verifications, lang="zh", params=None, figures_dir=None,
                    sample_context=None, interpretation_is_llm=True):
    """写出最终 md，返回路径。

    interpretation_is_llm=False 时说明 2.1 拿到的只是确定性摘要（没配大模型），
    此时不重复贴一遍——完整摘要 2.3 里已经有了。
    """
    os.makedirs(os.path.dirname(md_path), exist_ok=True)
    md_dir = os.path.dirname(md_path)
    case, control = ev.design_label(design_json)
    zh = lang == "zh"
    L = []

    # ---------- 头 ----------
    L.append("# PB 双键定位分析报告" if zh else "# PB double-bond localization report")
    L.append("")
    L.append(f"*{'生成时间' if zh else 'Generated'}: "
             f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}*")
    L.append("")
    if sample_context:
        L.append(f"> **{'研究背景' if zh else 'Study context'}**：{sample_context}")
        L.append("")

    # ---------- 0 关键结论（先给能用的） ----------
    L.append("## 0. " + ("关键结论" if zh else "Key conclusions"))
    L.append("")
    L.append(summarize_for_chat(verifications, design_json, lang, tables))
    L.append("")

    # ---------- 1 实验设计 ----------
    L.append("## 1. " + ("实验设计" if zh else "Design"))
    L.append("")
    by_group = {}
    for s, g in (design_json or {}).get("groups", {}).items():
        if g:
            by_group.setdefault(g, []).append(s)
    L.append(f"| {'组' if zh else 'Group'} | {'角色' if zh else 'Role'} | n | "
             f"{'样本' if zh else 'Samples'} |")
    L.append("|---|---|---|---|")
    for g in sorted(by_group):
        role = "case" if g == case else ("control" if g == control else "—")
        L.append(f"| {g} | {role} | {len(by_group[g])} | {', '.join(sorted(by_group[g]))} |")
    if (design_json or {}).get("excluded"):
        L.append(f"| — | {'排除' if zh else 'excluded'} | "
                 f"{len(design_json['excluded'])} | {', '.join(design_json['excluded'])} |")
    L.append("")
    L.append(f"log2FC > 0 {'表示在' if zh else 'means higher in'} **{case}** {'中更高' if zh else ''}。")
    if (design_json or {}).get("note"):
        L.append("")
        L.append(f"> {design_json['note']}")
    L.append("")

    # ---------- 2 质控与差异（大模型解读） ----------
    L.append("## 2. " + ("质控与差异分析" if zh else "QC and differential analysis"))
    L.append("")
    L.append("### 2.1 " + ("解读" if zh else "Interpretation"))
    L.append("")
    if interpretation and interpretation_is_llm:
        L.append(_demote_headings(interpretation))
    else:
        L.append(("*未配置大模型，本节无解读；下面 2.2 的统计现实与 2.3 的完整数据摘要"
                  "是确定性计算的，照样可用。*") if zh else
                 ("*No LLM configured, so there is no narrative review here; 2.2 and 2.3 below "
                  "are computed deterministically and remain valid.*"))
    L.append("")
    L.append("### 2.2 " + ("统计现实" if zh else "Statistical reality"))
    L.append("")
    L.append("```")
    L.append(ev.power_reality(tables))
    L.append("```")
    L.append("")
    L.append("### 2.3 " + ("数据摘要" if zh else "Digest"))
    L.append("")
    L.append("<details><summary>" + ("展开完整数据摘要" if zh else "Expand full digest") +
             "</summary>")
    L.append("")
    L.append("```")
    L.append(ev.full_digest(tables, design_json))
    L.append("```")
    L.append("")
    L.append("</details>")
    L.append("")
    if figures_dir:
        for name, cap in (("Fig1_position_dimension", "位置维度（PB 独有）" if zh else
                           "Position dimension (PB-only)"),
                          ("Fig2_quality_control", "质控" if zh else "Quality control")):
            for ext in ("png", "svg"):
                p = os.path.join(figures_dir, f"{name}.{ext}")
                if os.path.exists(p):
                    L += _img(md_dir, p, cap)
                    break

    # ---------- 3 通路 / LIPID MAPS ----------
    L.append("## 3. " + ("LIPID MAPS 联合分析" if zh else "LIPID MAPS joint analysis"))
    L.append("")
    ctx = pathway_ctx or {}
    any_pw = False
    for layer, title in (("pb", "PB 层（双键位置份额变化）" if zh else "PB layer (position share)"),
                         ("species", "常规层（物质总量变化）" if zh else "Conventional layer (species)")):
        rows = (ctx.get(layer) or {}).get("enrichment") or []
        if not rows:
            continue
        any_pw = True
        L.append(f"### {title}")
        L.append("")
        L.append(f"| {'通路' if zh else 'Pathway'} | {'命中/特征' if zh else 'hits/features'} | "
                 f"mean log2FC | p | {'涉及类别' if zh else 'subclasses'} |")
        L.append("|---|---|---|---|---|")
        for r in rows[:8]:
            nm = r.get("name_zh") if zh else r.get("name_en")
            L.append(f"| {nm} | {r.get('n_sig')}/{r.get('n_features')} | "
                     f"{ev._fmt(r.get('mean_log2fc'))} | {ev._fmt(r.get('enrich_p'))} | "
                     f"{r.get('subclasses')} |")
        L.append("")
        plot = (ctx.get(layer) or {}).get("plot")
        L += _img(md_dir, plot, f"{title} — pathway enrichment")
    online = ctx.get("online") or {}
    if online:
        any_pw = True
        L.append("### " + ("LIPID MAPS 在线注释" if zh else "LIPID MAPS online annotation"))
        L.append("")
        L.append(f"| {'类别' if zh else 'Subclass'} | LM_ID | {'链接' if zh else 'Link'} |")
        L.append("|---|---|---|")
        for k, v in online.items():
            L.append(f"| {k} | {v.get('lm_id')} | [{v.get('lm_id')}]({v.get('url')}) |")
        L.append("")
        L.append(f"*{'经 lipidmaps.org REST 接口实时获取' if zh else 'Fetched live from the lipidmaps.org REST API'}*")
        L.append("")
    if not any_pw:
        L.append("（未获得通路信息）" if zh else "(no pathway information)")
        L.append("")

    # ---------- 4 假设与验证 ----------
    L.append("## 4. " + ("科学假设与数据驱动验证" if zh else "Hypotheses and data-driven verification"))
    L.append("")
    if not verifications:
        L.append("（没有假设）" if zh else "(none)")
    else:
        L.append(("每条假设由一个独立的验证子 agent 复核。子 agent 只看该假设名下的真实数据行，"
                  "并被要求优先反驳。四个评分相互独立，不合并为总分。") if zh else
                 ("Each hypothesis was reviewed by an independent sub-agent that saw only that "
                  "hypothesis's real data rows and was asked to refute it. The four scores are "
                  "independent and deliberately not merged."))
        L.append("")
        # 排序表
        L.append(f"| # | {'裁决' if zh else 'Verdict'} | {'可信度' if zh else 'Confidence'} | "
                 f"{'注释' if zh else 'Annot.'} | {'强度' if zh else 'Intensity'} | "
                 f"{'新颖度' if zh else 'Novelty'} | {'推荐' if zh else 'Recommendation'} |")
        L.append("|---|---|---|---|---|---|---|")
        for v in verifications:
            h, dv, nov = v["hypothesis"], v["data_verdict"], v["novelty"]
            vd = VERDICT_ZH.get(dv["verdict"], dv["verdict"]) if zh else dv["verdict"]
            L.append(f"| {h.get('id')} | {vd} | {_score_bar(v['confidence_score'])} "
                     f"{v['confidence_score']:.2f} | {v['annotation_confidence']['level']} | "
                     f"{v['intensity_confidence']['level']} | "
                     f"{'新' if nov.get('is_novel') else '已报道'}"
                     f"（{nov.get('n_references', 0)} 篇） | {v['recommendation']} |")
        L.append("")

        for v in verifications:
            h, dv = v["hypothesis"], v["data_verdict"]
            nov, ann, inten = v["novelty"], v["annotation_confidence"], v["intensity_confidence"]
            L.append(f"### {'假设' if zh else 'Hypothesis'} {h.get('id')}")
            L.append("")
            L.append(f"**{h.get('statement')}**")
            L.append("")
            if h.get("biological_mechanism"):
                L.append(f"- {'提出的机制' if zh else 'Proposed mechanism'}："
                         f"{h['biological_mechanism']}")
            L.append(f"- {'依据' if zh else 'Basis'}：{h.get('rationale', '—')}")
            L.append(f"- {'涉及脂质' if zh else 'Lipids'}："
                     f"{', '.join(h.get('involved_lipids') or []) or '—'}"
                     f"｜{'类别' if zh else 'classes'}："
                     f"{', '.join(h.get('involved_classes') or []) or '—'}"
                     f"｜{'位置' if zh else 'positions'}："
                     f"{', '.join(h.get('involved_positions') or []) or '—'}")
            if h.get("pathways"):
                L.append(f"- {'通路' if zh else 'Pathways'}：{', '.join(h['pathways'])}")
            L.append("")
            L.append(f"**{'子 agent 裁决' if zh else 'Sub-agent verdict'}："
                     f"{VERDICT_ZH.get(dv['verdict'], dv['verdict']) if zh else dv['verdict']}"
                     f"**（score {dv['score']:.2f}，{dv.get('source')}）")
            L.append("")
            if dv.get("comment"):
                L.append(f"> {dv['comment']}")
                L.append("")
            if dv.get("contradictions"):
                L.append(f"{'反驳点' if zh else 'Contradictions'}：")
                for c in dv["contradictions"]:
                    L.append(f"- ⚠ {c}")
                L.append("")
            if dv.get("what_would_settle_it"):
                L.append(f"{'如何定论' if zh else 'What would settle it'}："
                         f"{dv['what_would_settle_it']}")
                L.append("")
            L.append(f"| {'维度' if zh else 'Axis'} | {'评分' if zh else 'Score'} | "
                     f"{'依据' if zh else 'Detail'} |")
            L.append("|---|---|---|")
            L.append(f"| {'注释可信度' if zh else 'Annotation'} | {_score_bar(ann['score'])} "
                     f"{ann['score']:.2f}（{ann['level']}） | {ann.get('detail', '')} |")
            L.append(f"| {'强度可信度' if zh else 'Intensity'} | {_score_bar(inten['score'])} "
                     f"{inten['score']:.2f}（{inten['level']}） | {inten.get('detail', '')} |")
            L.append(f"| {'是否够新' if zh else 'Novelty'} | {_score_bar(nov['score'])} "
                     f"{nov['score']:.2f} | {nov.get('detail', '')} |")
            L.append("")
            refs = nov.get("references") or []
            if refs:
                L.append(f"{'相关报道' if zh else 'Related reports'}：")
                for r in refs:
                    doi = f" doi:[{r.get('doi')}](https://doi.org/{r.get('doi')})" if r.get("doi") else ""
                    L.append(f"- {r.get('title')} — *{r.get('journal')}* {r.get('year')}"
                             f"（{'被引' if zh else 'cited'} {r.get('cited_by', 0)}，"
                             f"{r.get('source')}）{doi}")
                L.append("")
            L.append("<details><summary>" +
                     ("子 agent 看到的数据切片" if zh else "Data slice seen by the sub-agent") +
                     "</summary>")
            L.append("")
            L.append("```")
            L.append(v.get("data_slice") or "")
            L.append("```")
            L.append("")
            L.append("</details>")
            L.append("")

    # ---------- 5 附录 ----------
    L.append("## 5. " + ("附录：参数与产物" if zh else "Appendix: parameters and outputs"))
    L.append("")
    L.append("```")
    for k, v in (params or {}).items():
        L.append(f"{k} = {v}")
    L.append("```")
    L.append("")
    L.append(("评分口径：\n"
              "  注释可信度 = 0.40×master匹配率 + 0.35×MS2注释率 + 0.25×MS2样本覆盖\n"
              "  强度可信度 = 0.35×强度分位 + 0.30×定量样本覆盖 + 0.15×离子对占比 + 0.20×两针相关\n"
              "  是否够新   = 文献命中数（0 篇→0.85，1-2→0.6，3-5→0.4，>5→0.2）\n"
              "  可信度     = 0.40×子agent裁决 + 0.35×注释 + 0.25×强度（**不含新颖度**）")
             if zh else
             ("Scoring:\n"
              "  Annotation = 0.40*master-matched + 0.35*MS2-annotated + 0.25*MS2 coverage\n"
              "  Intensity  = 0.35*intensity pct + 0.30*quant coverage + 0.15*pair share + 0.20*replicate r\n"
              "  Novelty    = literature hits (0->0.85, 1-2->0.6, 3-5->0.4, >5->0.2)\n"
              "  Confidence = 0.40*sub-agent + 0.35*annotation + 0.25*intensity (novelty excluded)"))
    L.append("")

    text = "\n".join(L)
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(text)
    return md_path


def summarize_for_chat(verifications, design_json=None, lang="zh", tables=None):
    """凝练成对话窗口里的几行——只说能决定下一步的话。"""
    zh = lang == "zh"
    case, control = ev.design_label(design_json)
    if not verifications:
        return "（没有可凝练的结论）" if zh else "(nothing to summarize)"

    ranked = sorted(verifications, key=lambda v: v["confidence_score"], reverse=True)
    n_sup = sum(1 for v in verifications
                if v["data_verdict"]["verdict"] in ("supported", "partly_supported"))
    n_novel = sum(1 for v in verifications if (v["novelty"] or {}).get("is_novel"))
    lines = []
    if tables is not None:
        db = tables.get("dbpos_diff")
        n_sig = int((db["share_q"] < ev.ALPHA).sum()) if db is not None and len(db) else 0
        if n_sig == 0:
            lines.append(("⚠ 前提：这批数据在 FDR<0.05 下**没有任何显著项**（小样本功效上限），"
                          "所以下面全部是效应量层面的探索性结论，需要扩样本才能定论。")
                         if zh else
                         ("⚠ Nothing reached FDR<0.05 at this sample size; everything below is "
                          "exploratory effect-size evidence."))
            lines.append("")
    lines.append((f"{case} vs {control}：{len(verifications)} 条假设经独立子 agent 验证，"
                  f"{n_sup} 条获数据支持（含部分支持），{n_novel} 条文献未见报道。")
                 if zh else
                 (f"{case} vs {control}: {len(verifications)} hypotheses verified, "
                  f"{n_sup} supported, {n_novel} without prior reports."))
    lines.append("")
    top = [v for v in ranked if v["data_verdict"]["verdict"] != "not_supported"][:3]
    if top:
        lines.append("**" + ("最值得跟进的" if zh else "Most worth pursuing") + "**：")
        for v in top:
            h, nov = v["hypothesis"], v["novelty"]
            lines.append(f"{h.get('id')}. {h.get('statement')}")
            lines.append(f"   → {'可信度' if zh else 'confidence'} {v['confidence_score']:.2f}"
                         f"（{'注释' if zh else 'annot'} {v['annotation_confidence']['level']}/"
                         f"{'强度' if zh else 'int'} {v['intensity_confidence']['level']}）"
                         f"｜{'文献' if zh else 'lit'} {nov.get('n_references', 0)}"
                         f"｜{v['recommendation']}")
    dropped = [v for v in ranked if v["data_verdict"]["verdict"] == "not_supported"]
    if dropped:
        lines.append("")
        lines.append(("**被子 agent 反驳的**：" if zh else "**Refuted by the sub-agent**: ") +
                     ", ".join(f"#{v['hypothesis'].get('id')}" for v in dropped))
    return "\n".join(lines)
