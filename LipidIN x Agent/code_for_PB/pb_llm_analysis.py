# -*- coding: utf-8 -*-
"""
PB 工作流的大模型层：解析实验设计、解读质控/差异、提出科学假设、接文献与 LIPID MAPS。

复用主仓库 `../../code` 里已经写好的三块（不重复造轮子，也避免两处漂移）：
  literature_search.find_similar_reports —— Europe PMC，失败回退 CrossRef，**不需要 key**，
      网络不通时返回 [] 而不是抛异常，所以离线也能跑完。
  lipid_pathways —— 通路富集用**本地策展表**（离线可用）；online_enrich 才会真的去打
      LIPID MAPS 的 REST 接口（https://www.lipidmaps.org/rest/compound/abbrev/...），
      这就是"接上 lipidmaps 在线版本"的那一步。
  llm_client.LLMClient —— 多提供方（anthropic/openai/deepseek/custom），与主界面同一套配置。

没配大模型也能跑：设计解析回退到样本名前缀，假设生成回退到基于效应量的规则引擎。
只是自然语言那部分（"用A1到A3当对照"）就用不了。
"""

import os
import re
import sys
import json
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
# 主仓库的 code 目录追加到**末尾**：PB 自己的 run.py/quantify.py 等同名模块优先，
# 只有 PB 这边没有的（literature_search / lipid_pathways / llm_client）才落到主仓库。
MAIN_CODE = CODE_DIR.parent.parent / "code"
if str(MAIN_CODE) not in sys.path:
    sys.path.append(str(MAIN_CODE))

import pb_evidence as ev

try:
    import literature_search
except Exception:
    literature_search = None
try:
    import lipid_pathways
except Exception:
    lipid_pathways = None


MAX_HYPOTHESES = 5


def _log(log, msg):
    (log or print)(msg)


# ------------------------------------------------------------------ #
# 0. 自由对话意图解析：用户随便说一句话，落成"改哪个参数 / 跑哪个阶段"
# ------------------------------------------------------------------ #
INTENT_SYSTEM = """You are the command interpreter for a Paternò–Büchi (PB) lipidomics agent.
The user speaks freely (usually Chinese). Turn their message into a concrete intent.

Return ONLY JSON:
{"config": {"<setting key>": <value>, ...},
 "action": "<stage key | full | status | help | show_config | none>",
 "ask": "<if you need one specific missing thing (e.g. a path), the question to ask; else empty>",
 "reply": "<one short sentence in the user's language confirming what you understood>"}

The ONLY valid setting keys are these (type, range and meaning given):
{catalog}

The ONLY valid action values are:
{actions}

Rules:
- Put a key in "config" ONLY if the user actually supplied a value. Never guess a path,
  never invent a number. If they clearly want to change something but gave no value
  (e.g. "我要改文件地址"), leave "config" empty and put the concrete question in "ask".
- Paths: copy them EXACTLY as written, including drive letters and backslashes. Do not
  normalize, quote, or complete them.
- "action" is what to RUN. Merely mentioning a stage is not a request to run it — only set
  "action" when the user is clearly asking to run/redo something. Otherwise use "none".
- Changing a parameter is NOT a request to re-run: if the user says "把 ppm 改成 10",
  return config only, action "none".
- If the user asks what the current settings are, action "show_config".
- Keep "reply" factual: restate the change. Do not add advice."""


def parse_intent(llm, text, catalog, actions, lang="zh", log=None):
    """把一句自由文本变成 {config, action, ask, reply}。没有大模型就返回 None。"""
    if llm is None or not llm.is_enabled() or not (text or "").strip():
        return None
    system = (INTENT_SYSTEM.replace("{catalog}", catalog)
              .replace("{actions}", ", ".join(actions)))
    try:
        data = llm.chat_json(system, f"User message ({'Chinese' if lang == 'zh' else 'English'}):\n{text}")
    except Exception as e:
        _log(log, f"[intent] 大模型解析失败：{type(e).__name__}: {e}")
        return None
    if not isinstance(data, dict):
        return None
    action = str(data.get("action") or "none")
    if action not in actions:
        action = "none"
    cfg = data.get("config")
    return {"config": cfg if isinstance(cfg, dict) else {},
            "action": action,
            "ask": str(data.get("ask") or ""),
            "reply": str(data.get("reply") or "")}


def grab_path(text):
    """从一句话里抠出第一个绝对路径（Windows 盘符 / UNC / POSIX），没有就 None。

    路径**绝不交给大模型转述**——模型会顺手规范化、补全、甚至把反斜杠吃掉，而路径错一个
    字符就是 FileNotFoundError。所以先用正则确定性地抠，抠到就以正则的为准。
    """
    if not text:
        return None
    m = re.search(r'["\']([A-Za-z]:[\\/][^"\']*|\\\\[^"\']+|/[^"\']+)["\']', text)
    if m:
        return m.group(1).strip()
    m = re.search(r'([A-Za-z]:[\\/][^\s，,；;]*|\\\\[^\s，,；;]+|/(?:[^\s/]+/)+[^\s/]*)', text)
    if m:
        return m.group(1).strip().strip('"').strip("'").rstrip("。.,，；;、)）」』】")
    return None


# ------------------------------------------------------------------ #
# 1. 实验设计：用哪些文件、跟哪些比
# ------------------------------------------------------------------ #
DESIGN_SYSTEM = """You map lipidomics sample names to an experimental design.

You are given the list of available sample names from the quantification table and a
free-form instruction (Chinese or English) describing which samples to use and what to
compare against what.

Return ONLY JSON:
{
  "groups": {"<sample name>": "<group label or null>"},
  "case": "<group label>",
  "control": "<group label>",
  "excluded": ["<sample name>", ...],
  "note": "<one short sentence, in the user's language, restating the design>"
}

Rules:
- Every key in "groups" MUST be an exact string from the provided sample list. Never invent names.
- Samples the user excludes, or that are QC/pool/MSMS-only runs, map to null and appear in "excluded".
- "case" and "control" MUST be group labels that appear as values in "groups".
- log2FC is computed as case over control, so "case" is the group the user treats as the
  condition/treatment and "control" is the reference/baseline.
- If the instruction is ambiguous, choose the most literal reading and say so in "note"."""


def parse_design_spec(llm, samples, text, lang="zh", log=None):
    """把'用哪些文件、跟谁比'的自然语言解析成设计。

    返回 {"groups": {...}, "case", "control", "excluded", "note", "source"}。
    没有大模型时回退到样本名前缀（A1->A），与 pb_qc_diff 的历史约定一致。
    """
    if llm is not None and llm.is_enabled() and (text or "").strip():
        user = (f"Available sample names (exact strings):\n{json.dumps(samples, ensure_ascii=False)}\n\n"
                f"Instruction from the user:\n{text}\n\n"
                f"Reply language for 'note': {'Chinese' if lang == 'zh' else 'English'}")
        try:
            data = llm.chat_json(DESIGN_SYSTEM, user)
            design = _sanitize_design(data, samples)
            if design:
                design["source"] = "llm"
                return design
            _log(log, "[design] 大模型返回的设计不合法，回退到样本名前缀")
        except Exception as e:
            _log(log, f"[design] 大模型解析失败（{type(e).__name__}: {e}），回退到样本名前缀")
    return _design_from_prefix(samples, note_lang=lang)


def _sanitize_design(data, samples):
    """只信任真实存在的样本名和自洽的 case/control——模型编的一律丢掉。"""
    if not isinstance(data, dict):
        return None
    raw = data.get("groups") or {}
    if not isinstance(raw, dict):
        return None
    valid = set(samples)
    groups = {}
    for k, v in raw.items():
        if k in valid:
            groups[k] = None if v in (None, "", "null", "None") else str(v)
    if not groups:
        return None
    labels = {v for v in groups.values() if v}
    case, control = data.get("case"), data.get("control")
    if case not in labels or control not in labels or case == control:
        # 模型给的 case/control 不在组里就自己挑：出现的前两个标签，字母序稳定
        ordered = sorted(labels)
        if len(ordered) < 2:
            return None
        case, control = ordered[0], ordered[1]
    excluded = [s for s in samples if not groups.get(s)]
    return {"groups": groups, "case": str(case), "control": str(control),
            "excluded": excluded, "note": str(data.get("note") or "")}


def _design_from_prefix(samples, note_lang="zh"):
    """回退：样本名前缀即组名（A1/A2->A，B119->B），认不出的排除。"""
    import re
    groups, labels = {}, []
    for s in samples:
        m = re.match(r"^([A-Za-z]+)(\d+)$", s)
        g = m.group(1) if m else None
        groups[s] = g
        if g and g not in labels:
            labels.append(g)
    labels = sorted(labels)
    case = labels[0] if labels else "A"
    control = labels[1] if len(labels) > 1 else "B"
    note = (f"按样本名前缀自动分组：{case} vs {control}" if note_lang == "zh"
            else f"Grouped by sample-name prefix: {case} vs {control}")
    return {"groups": groups, "case": case, "control": control,
            "excluded": [s for s in samples if not groups.get(s)],
            "note": note, "source": "prefix"}


def design_text(design, lang="zh"):
    """把设计渲染成给用户看的确认文本。"""
    by_group = {}
    for s, g in (design.get("groups") or {}).items():
        if g:
            by_group.setdefault(g, []).append(s)
    case, control = design.get("case"), design.get("control")
    if lang == "zh":
        lines = [f"实验设计：**{case}** vs **{control}**（log2FC>0 表示在 {case} 中更高）"]
        for g in sorted(by_group):
            role = "（case）" if g == case else ("（control）" if g == control else "（不参与比较）")
            lines.append(f"  · {g}{role}：{len(by_group[g])} 个 —— {', '.join(sorted(by_group[g]))}")
        if design.get("excluded"):
            lines.append(f"  · 排除：{', '.join(design['excluded'])}")
        if design.get("note"):
            lines.append(f"  说明：{design['note']}")
    else:
        lines = [f"Design: **{case}** vs **{control}** (log2FC>0 means higher in {case})"]
        for g in sorted(by_group):
            role = " (case)" if g == case else (" (control)" if g == control else " (not compared)")
            lines.append(f"  · {g}{role}: {len(by_group[g])} — {', '.join(sorted(by_group[g]))}")
        if design.get("excluded"):
            lines.append(f"  · Excluded: {', '.join(design['excluded'])}")
        if design.get("note"):
            lines.append(f"  Note: {design['note']}")
    return "\n".join(lines)


# ------------------------------------------------------------------ #
# 2. 质控 + 差异分析的大模型解读
# ------------------------------------------------------------------ #
INTERPRET_SYSTEM = """You are a lipidomics QC and statistics reviewer. You are given the
real numbers from a Paternò–Büchi (PB) double-bond localization run.

Write a short, honest review in {lang}. Cover, in this order:
1. QC verdict — is this data usable? Point at the actual numbers (missing rate, total-area
   spread, injection replicate correlation). Say plainly if something looks bad.
2. What the differential analysis can and cannot support at this sample size. If the digest
   says nothing reached FDR<0.05, say so in the FIRST sentence of this part and do not
   describe anything as "significant".
3. Where the real signal is (effect sizes, position-preference flips) and which specific
   lipids carry it.

Hard rules:
- Never claim statistical significance the digest does not show.
- Quote concrete numbers from the digest. No generic lipidomics boilerplate.
- Be concise: at most ~250 words. No headings, no bullet spam — short paragraphs.
- If the data looks weak, saying so IS the useful answer."""


def interpret_results(llm, tables, design_json=None, lang="zh", log=None):
    """让大模型对质控与差异分析给一段人话解读。没有大模型就给确定性摘要。"""
    digest = ev.full_digest(tables, design_json)
    if llm is None or not llm.is_enabled():
        head = ("（未配置大模型，下面是确定性摘要；配置后可得到解读）\n\n" if lang == "zh"
                else "(No LLM configured; deterministic digest below.)\n\n")
        return head + digest
    system = INTERPRET_SYSTEM.replace("{lang}", "Chinese" if lang == "zh" else "English")
    try:
        return llm.chat(system, f"Data digest:\n\n{digest}")
    except Exception as e:
        _log(log, f"[interpret] 大模型失败：{type(e).__name__}: {e}")
        return digest


# ------------------------------------------------------------------ #
# 3. LIPID MAPS 联合分析
# ------------------------------------------------------------------ #
def pathway_context(tables, output_dir, lang="zh", online=True, log=None):
    """通路富集（本地策展表）+ LIPID MAPS 在线注释（REST）。

    两层都做：PB 层（位置份额变了的）和常规层（物质总量变了的），因为"同一通路上常规
    看不到、PB 才看到"本身就是结论的一部分。
    """
    if lipid_pathways is None:
        return {"error": "lipid_pathways 不可用", "pb": {}, "species": {}, "online": {}}
    out = {"pb": {}, "species": {}, "online": {}, "subclasses": []}
    os.makedirs(output_dir, exist_ok=True)
    for layer in ("pb", "species"):
        frame = ev.pathway_frame(tables, layer=layer)
        if not len(frame):
            continue
        try:
            out[layer] = lipid_pathways.analyze_pathways(
                frame, output_dir, tag=f"PB_{layer}", language=("zh" if lang == "zh" else "en"),
                log=(log or print))
        except Exception as e:
            _log(log, f"[pathway] {layer} 层富集失败：{type(e).__name__}: {e}")
            out[layer] = {}
    # 在线：把出现过的 subclass 拿去 LIPID MAPS 换 LM_ID / 正式名 / 链接
    db = tables.get("dbpos_diff")
    subs = []
    if db is not None and len(db):
        for s in db["subclass"].dropna().unique():
            norm = lipid_pathways.normalize_subclass(s)
            if norm and norm not in subs:
                subs.append(norm)
    out["subclasses"] = subs
    if online and subs:
        try:
            _log(log, f"[pathway] 连 LIPID MAPS 在线接口，查询 {len(subs)} 个 subclass ……")
            out["online"] = lipid_pathways.online_enrich(subs, log=(log or print)) or {}
            _log(log, f"[pathway] LIPID MAPS 返回 {len(out['online'])} 条")
        except Exception as e:
            _log(log, f"[pathway] LIPID MAPS 在线查询失败（不影响后续）：{type(e).__name__}: {e}")
            out["online"] = {}
    return out


def pathway_digest_text(ctx, lang="zh"):
    """通路上下文渲染成给大模型/报告用的文本。"""
    if not ctx:
        return "（无通路信息）"
    lines = []
    for layer, title in (("pb", "PB 层（双键位置份额变化）"), ("species", "常规层（物质总量变化）")):
        res = ctx.get(layer) or {}
        rows = res.get("enrichment") or []
        if not rows:
            continue
        lines.append(f"{title} 富集到的通路：")
        for r in rows[:6]:
            name = r.get("name_zh") if lang == "zh" else r.get("name_en")
            lines.append(f"  · {name}：{r.get('n_sig')}/{r.get('n_features')} 个特征命中，"
                         f"涉及 {r.get('subclasses')}，mean log2FC={r.get('mean_log2fc')}")
    online = ctx.get("online") or {}
    if online:
        lines.append("LIPID MAPS 在线注释：" +
                     "；".join(f"{k} -> {v.get('lm_id')} {v.get('name')}"
                               for k, v in list(online.items())[:10]))
    return "\n".join(lines) if lines else "（无通路信息）"


# ------------------------------------------------------------------ #
# 4. 科学假设
# ------------------------------------------------------------------ #
HYPOTHESIS_SYSTEM = """You are a lipid biologist reading the results of a Paternò–Büchi (PB)
double-bond localization experiment. PB reveals WHERE the C=C sits inside a fatty acyl chain —
information ordinary lipidomics cannot see. The share of each double-bond position is
normalized within a species, so it is independent of how much of that species is present.

Propose at most {max_n} scientific hypotheses that this specific dataset actually supports.

Return ONLY JSON:
{"hypotheses": [
  {"statement": "<one falsifiable sentence in {lang}>",
   "rationale": "<why this data supports it; cite the actual numbers>",
   "involved_lipids": ["<total_chain_name exactly as in the data, e.g. 'PC 40:5'>"],
   "involved_classes": ["<subclass, e.g. 'PC'>"],
   "involved_positions": ["<e.g. 'n-9'>"],
   "evidence_tier": "<position_flip|position_only|both_change|omega|species_only>",
   "biological_mechanism": "<the enzyme/pathway that would explain it, e.g. SCD1 vs FADS2>",
   "search_query": "<English literature query, 3-8 words, no lipid-name punctuation>"}]}

Hard rules:
- Ground every hypothesis in lipids that APPEAR IN THE DIGEST. Never invent a lipid name.
- Obey the "统计现实" section. If nothing is FDR-significant, your statements must be phrased
  as exploratory ("suggests", "is consistent with") and must never say "significant".
- Prefer hypotheses that only PB could raise (position flips, share changes) over ones ordinary
  lipidomics already answers (total-amount changes).
- The mechanism matters: a double-bond position shift usually points at a specific desaturase
  (SCD1 -> n-9 from 18:0; FADS2 -> n-6/n-3 series; elongase-driven n-x shifts). Name the enzyme.
- Fewer, better hypotheses beat more. If the data only supports 2, return 2.
- 'search_query' must be about the BIOLOGY (e.g. "SCD1 n-9 desaturation cancer"), not about
  the method."""


def generate_hypotheses(llm, tables, design_json=None, pathway_ctx=None, lang="zh",
                        max_n=MAX_HYPOTHESES, log=None):
    """提出科学假设，并给每条检索文献 + 判断新颖度。"""
    digest = ev.full_digest(tables, design_json)
    pathways = pathway_digest_text(pathway_ctx, lang) if pathway_ctx else ""
    hyps = []
    if llm is not None and llm.is_enabled():
        system = (HYPOTHESIS_SYSTEM.replace("{max_n}", str(max_n))
                  .replace("{lang}", "Chinese" if lang == "zh" else "English"))
        user = f"Data digest:\n\n{digest}\n\nPathway context:\n{pathways}"
        try:
            data = llm.chat_json(system, user)
            raw = data.get("hypotheses") if isinstance(data, dict) else data
            hyps = _sanitize_hypotheses(raw, tables, max_n)
        except Exception as e:
            _log(log, f"[hypothesis] 大模型失败（{type(e).__name__}: {e}），回退到规则引擎")
    if not hyps:
        hyps = _hypotheses_from_rules(tables, design_json, lang, max_n)
    for i, h in enumerate(hyps, start=1):
        h["id"] = i
        refs = find_literature(h.get("search_query") or h.get("statement", ""), log=log)
        h["references"] = refs
        h["novelty"] = ev.novelty_score(refs, lang)
        h["pathways"] = _pathways_for(h, pathway_ctx, lang)
    return hyps


def _sanitize_hypotheses(raw, tables, max_n):
    """丢掉模型编的脂质名——数据里不存在的名字一律剔除，剔空了就丢掉整条假设。"""
    if not isinstance(raw, list):
        return []
    db = tables.get("dbpos_diff")
    known_lipids = set(db["total_chain_name"].dropna().astype(str)) if db is not None else set()
    known_classes = set(db["subclass"].dropna().astype(str)) if db is not None else set()
    out = []
    for item in raw:
        if not isinstance(item, dict) or not str(item.get("statement") or "").strip():
            continue
        lipids = [x for x in (item.get("involved_lipids") or []) if str(x) in known_lipids]
        classes = [x for x in (item.get("involved_classes") or []) if str(x) in known_classes]
        if not lipids and not classes:
            continue   # 一个真实脂质都对不上 -> 这条是编的，不要
        out.append({
            "statement": str(item.get("statement")).strip(),
            "rationale": str(item.get("rationale") or "").strip(),
            "involved_lipids": lipids,
            "involved_classes": classes,
            "involved_positions": [str(x) for x in (item.get("involved_positions") or [])],
            "evidence_tier": str(item.get("evidence_tier") or "position_flip"),
            "biological_mechanism": str(item.get("biological_mechanism") or "").strip(),
            "search_query": str(item.get("search_query") or "").strip(),
            "source": "llm",
        })
        if len(out) >= max_n:
            break
    return out


def _hypotheses_from_rules(tables, design_json, lang="zh", max_n=MAX_HYPOTHESES):
    """没有大模型时的回退：直接把效应量最大的位置翻转写成假设。"""
    case, control = ev.design_label(design_json)
    pref = tables.get("position_preference")
    out = []
    if pref is None or not len(pref):
        return out
    flips = pref[pref["preference_flip"] == True].head(max_n)  # noqa: E712
    for _, r in flips.iterrows():
        sub = str(r["total_chain_name"]).split()[0]
        if lang == "zh":
            st = (f"{r['total_chain_name']} 在 {case} 组中双键更偏向 {r['dominant_A']}，"
                  f"而在 {control} 组中偏向 {r['dominant_B']}，提示两组间该物质的去饱和/延长"
                  f"路径不同（探索性，最大份额差 {r['max_abs_delta']:.0%}）。")
            ra = (f"位置偏好表显示该 (物质,f={r['f']}) 层优势位置在两组间翻转；"
                  f"该层不做假设检验，仅效应量。")
        else:
            st = (f"In {case}, {r['total_chain_name']} favors {r['dominant_A']} while in "
                  f"{control} it favors {r['dominant_B']}, suggesting different desaturation/"
                  f"elongation routing (exploratory; max share difference {r['max_abs_delta']:.0%}).")
            ra = f"Dominant position flips within the (species, f={r['f']}) layer; effect size only."
        out.append({
            "statement": st, "rationale": ra,
            "involved_lipids": [str(r["total_chain_name"])],
            "involved_classes": [sub],
            "involved_positions": [str(r["dominant_A"]), str(r["dominant_B"])],
            "evidence_tier": "position_flip",
            "biological_mechanism": "",
            "search_query": f"{sub} double bond position {r['dominant_A']} desaturase",
            "source": "rule",
        })
    return out


def _pathways_for(hypothesis, ctx, lang="zh"):
    """这条假设涉及的脂质类别命中了哪些 LIPID MAPS 通路。"""
    if not ctx or lipid_pathways is None:
        return []
    names = []
    for cls in hypothesis.get("involved_classes") or []:
        info = lipid_pathways.classify_subclass(cls)
        if not info:
            continue
        for pw in info.get("pathways") or []:
            for row in (ctx.get("pb") or {}).get("enrichment") or []:
                if row.get("key") == pw:
                    nm = row.get("name_zh") if lang == "zh" else row.get("name_en")
                    if nm and nm not in names:
                        names.append(nm)
    return names


# ------------------------------------------------------------------ #
# 5. 文献
# ------------------------------------------------------------------ #
_STOP = {"the", "and", "for", "with", "from", "that", "this", "are", "was", "were", "its",
         "into", "via", "role", "study", "analysis", "effect", "effects", "level", "levels",
         "human", "mice", "mouse", "cell", "cells", "acid", "acids", "lipid", "lipids"}


def _query_terms(query):
    """查询里真正有区分度的词。太短的（PC/PE）和万金油词（lipid/study）不算。"""
    words = re.findall(r"[a-z0-9][a-z0-9\-]{2,}", (query or "").lower())
    return [w for w in dict.fromkeys(words) if w not in _STOP]


def _is_relevant(ref, terms, min_hits):
    text = f"{ref.get('title', '')} {ref.get('abstract', '')}".lower()
    return sum(1 for t in terms if t in text) >= min_hits


def find_literature(query, limit=5, log=None):
    """检索相关报道，并按词命中过滤掉明显不相关的。

    为什么要自己再过滤一道：literature_search 用的是 Europe PMC 的 `sort="CITED desc"`，
    查询词一宽，返回的就是"任意沾边、但被引极高"的综述——实测用"PC double bond position
    n-15 desaturase"能搜出《Heart Disease and Stroke Statistics》这种完全无关的文献。
    这些垃圾命中会直接污染新颖度打分（明明没人报道过，却因为 5 篇无关文献被判成"已有报道"），
    所以宁可少给几篇，也不能让"是否够新"这个轴变成噪声。

    那个 sort 在主仓库的共享模块里，主界面的智能体也在用，故不动它，只在 PB 这侧过滤。
    """
    if literature_search is None or not (query or "").strip():
        return []
    try:
        raw = literature_search.find_similar_reports(query, limit=max(limit * 3, 10)) or []
    except Exception as e:
        _log(log, f"[literature] 「{query}」检索失败：{type(e).__name__}: {e}")
        return []
    terms = _query_terms(query)
    if not terms:
        return raw[:limit]
    # 命中一半以上的区分词才算相关，且至少两个——只中一个词多半是巧合
    min_hits = max(2, (len(terms) + 1) // 2)
    keep = [r for r in raw if _is_relevant(r, terms, min_hits)]
    _log(log, f"[literature] 「{query}」返回 {len(raw)} 篇，"
              f"按关键词命中(≥{min_hits}/{len(terms)})过滤后剩 {len(keep)} 篇")
    return keep[:limit]


# ------------------------------------------------------------------ #
# 6. 假设修正（多轮对话）
# ------------------------------------------------------------------ #
REFINE_SYSTEM = """You are revising scientific hypotheses from a PB double-bond localization
experiment, based on what the user just told you.

You get: the current hypotheses (JSON), the data digest, and the user's instruction.
Apply the user's instruction faithfully — they may want a hypothesis dropped, merged, narrowed,
re-worded, re-focused on a lipid class, or a new one added.

Return ONLY JSON in the same shape as the input:
{"hypotheses": [{"statement", "rationale", "involved_lipids", "involved_classes",
                 "involved_positions", "evidence_tier", "biological_mechanism", "search_query"}]}

Hard rules:
- Keep hypotheses the user did not ask you to change EXACTLY as they were.
- Still never invent lipid names, and still obey the 统计现实 section of the digest.
- If the user asks for something the data cannot support, return the hypotheses unchanged and
  explain the refusal in the "rationale" of the affected item."""


def refine_hypotheses(llm, hypotheses, user_text, tables, design_json=None,
                      pathway_ctx=None, lang="zh", log=None):
    """按用户的话修正假设。改完重新检索文献与新颖度（因为假设内容变了）。"""
    if llm is None or not llm.is_enabled():
        return None, ("需要先配置大模型才能用自然语言修正假设。" if lang == "zh"
                      else "Configure an LLM first to refine hypotheses in natural language.")
    digest = ev.full_digest(tables, design_json)
    slim = [{k: h.get(k) for k in ("statement", "rationale", "involved_lipids",
                                   "involved_classes", "involved_positions",
                                   "evidence_tier", "biological_mechanism", "search_query")}
            for h in hypotheses]
    user = (f"Current hypotheses:\n{json.dumps(slim, ensure_ascii=False, indent=1)}\n\n"
            f"Data digest:\n{digest}\n\nUser instruction:\n{user_text}")
    try:
        data = llm.chat_json(REFINE_SYSTEM, user)
        raw = data.get("hypotheses") if isinstance(data, dict) else data
        new = _sanitize_hypotheses(raw, tables, MAX_HYPOTHESES)
    except Exception as e:
        return None, (f"修正失败：{type(e).__name__}: {e}" if lang == "zh"
                      else f"Refine failed: {type(e).__name__}: {e}")
    if not new:
        return None, ("模型没有给出有效的修正结果（可能是要求超出了数据能支持的范围）。"
                      if lang == "zh" else "No valid revision returned.")
    old_q = {h.get("statement"): h for h in hypotheses}
    for i, h in enumerate(new, start=1):
        h["id"] = i
        prev = old_q.get(h["statement"])
        if prev and prev.get("references") is not None:
            h["references"], h["novelty"] = prev["references"], prev["novelty"]   # 没变就不重查
        else:
            refs = find_literature(h.get("search_query") or h.get("statement", ""), log=log)
            h["references"], h["novelty"] = refs, ev.novelty_score(refs, lang)
        h["pathways"] = _pathways_for(h, pathway_ctx, lang)
    return new, None


def hypotheses_text(hypotheses, lang="zh"):
    """假设列表渲染给用户确认。"""
    if not hypotheses:
        return "（没有假设）" if lang == "zh" else "(no hypotheses)"
    lines = []
    for h in hypotheses:
        nov = h.get("novelty") or {}
        tag = ("新" if nov.get("is_novel") else "已有报道") if lang == "zh" else \
              ("novel" if nov.get("is_novel") else "reported")
        lines.append(f"[{h['id']}] {h['statement']}")
        if h.get("biological_mechanism"):
            lines.append(f"     机制：{h['biological_mechanism']}" if lang == "zh"
                         else f"     Mechanism: {h['biological_mechanism']}")
        lines.append(f"     依据：{h.get('rationale', '')}" if lang == "zh"
                     else f"     Basis: {h.get('rationale', '')}")
        lines.append(f"     脂质：{', '.join(h.get('involved_lipids') or []) or '—'}"
                     f"｜文献 {nov.get('n_references', 0)} 篇（{tag}）"
                     + (f"｜通路：{', '.join(h['pathways'])}" if h.get("pathways") else ""))
    return "\n".join(lines)
