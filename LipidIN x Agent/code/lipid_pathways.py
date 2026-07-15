"""Map lipid subclasses onto the LIPID MAPS classification and metabolic pathways, then run a
pathway-level enrichment on the differential result.

把脂质亚类映射到 LIPID MAPS 分类体系与代谢通路，并对差异分析结果做通路层面的富集。

The offline core is a curated table linking each lipid subclass (PC, TG, Cer, ...) to its LIPID MAPS
category (Fatty Acyls / Glycerolipids / Glycerophospholipids / Sphingolipids / Sterols ...) and to the
classic biosynthetic / catabolic pathways it takes part in (de novo FA synthesis, the Kennedy pathway,
the Lands remodeling cycle, cardiolipin synthesis, sphingolipid de novo, sterol ester metabolism, ...).
`pathway_enrichment` aggregates the differential features per pathway and scores over-representation of
significant lipids with a hypergeometric test. `online_enrich` optionally queries the LIPID MAPS REST
API to attach canonical class records/URLs when the machine is online; it never blocks the offline path.

离线内核是一张策展表，把每个脂质亚类（PC、TG、Cer……）关联到其 LIPID MAPS 大类（脂肪酰/甘油脂/甘油
磷脂/鞘脂/固醇……）以及其参与的经典合成/分解代谢通路（脂肪酸从头合成、Kennedy 通路、Lands 重塑循环、
心磷脂合成、鞘脂从头合成、固醇酯代谢……）。`pathway_enrichment` 按通路汇总差异特征，并用超几何检验对
显著脂质的富集进行打分。`online_enrich` 可在联网时调用 LIPID MAPS REST API 附上规范的分类记录/链接，
但绝不阻塞离线流程。
"""

import os
import json
import urllib.request

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import hypergeom
except Exception:  # pragma: no cover - scipy is a hard dependency elsewhere
    hypergeom = None

from qc_analysis import _configure_cjk_font

_configure_cjk_font()


# --------------------------------------------------------------------------- #
# LIPID MAPS categories (the 8 top-level classes of the LIPID MAPS ontology)
# LIPID MAPS 大类（本体的 8 个顶层类别）
# --------------------------------------------------------------------------- #
CATEGORIES = {
    "FA": ("Fatty Acyls", "脂肪酰"),
    "GL": ("Glycerolipids", "甘油脂"),
    "GP": ("Glycerophospholipids", "甘油磷脂"),
    "SP": ("Sphingolipids", "鞘脂"),
    "ST": ("Sterol Lipids", "固醇脂"),
    "PR": ("Prenol Lipids", "异戊烯醇脂"),
    "SL": ("Saccharolipids", "糖脂"),
    "PK": ("Polyketides", "聚酮"),
}

# Curated subclass -> LIPID MAPS mapping.
# lm_main = LIPID MAPS main-class code, used to build a stable classification URL.
# 策展的亚类 -> LIPID MAPS 映射；lm_main 为主类代码，用于拼接稳定的分类 URL。
SUBCLASS_INFO = {
    # ---- Fatty acyls ----
    "FA":     {"category": "FA", "lm_main": "FA01", "name_en": "Fatty acids", "name_zh": "脂肪酸"},
    "FAHFA":  {"category": "FA", "lm_main": "FA0110", "name_en": "Branched fatty-acid esters", "name_zh": "支链脂肪酸酯"},
    "CAR":    {"category": "FA", "lm_main": "FA07", "name_en": "Acyl carnitines", "name_zh": "酰基肉碱"},
    "AcCa":   {"category": "FA", "lm_main": "FA07", "name_en": "Acyl carnitines", "name_zh": "酰基肉碱"},
    # ---- Glycerolipids ----
    "MG":     {"category": "GL", "lm_main": "GL01", "name_en": "Monoacylglycerols", "name_zh": "单酰甘油"},
    "DG":     {"category": "GL", "lm_main": "GL02", "name_en": "Diacylglycerols", "name_zh": "二酰甘油"},
    "TG":     {"category": "GL", "lm_main": "GL03", "name_en": "Triacylglycerols", "name_zh": "三酰甘油"},
    "MGDG":   {"category": "GL", "lm_main": "GL05", "name_en": "Monogalactosyldiacylglycerol", "name_zh": "单半乳糖二酰甘油"},
    "DGDG":   {"category": "GL", "lm_main": "GL05", "name_en": "Digalactosyldiacylglycerol", "name_zh": "双半乳糖二酰甘油"},
    # ---- Glycerophospholipids ----
    "PC":     {"category": "GP", "lm_main": "GP01", "name_en": "Phosphatidylcholines", "name_zh": "磷脂酰胆碱"},
    "PE":     {"category": "GP", "lm_main": "GP02", "name_en": "Phosphatidylethanolamines", "name_zh": "磷脂酰乙醇胺"},
    "PS":     {"category": "GP", "lm_main": "GP03", "name_en": "Phosphatidylserines", "name_zh": "磷脂酰丝氨酸"},
    "PG":     {"category": "GP", "lm_main": "GP04", "name_en": "Phosphatidylglycerols", "name_zh": "磷脂酰甘油"},
    "PGP":    {"category": "GP", "lm_main": "GP05", "name_en": "Phosphatidylglycerophosphates", "name_zh": "磷脂酰甘油磷酸"},
    "PI":     {"category": "GP", "lm_main": "GP06", "name_en": "Phosphatidylinositols", "name_zh": "磷脂酰肌醇"},
    "PIP":    {"category": "GP", "lm_main": "GP07", "name_en": "Phosphatidylinositol phosphates", "name_zh": "磷脂酰肌醇磷酸"},
    "PA":     {"category": "GP", "lm_main": "GP10", "name_en": "Phosphatidic acids", "name_zh": "磷脂酸"},
    "CL":     {"category": "GP", "lm_main": "GP12", "name_en": "Cardiolipins", "name_zh": "心磷脂"},
    "BMP":    {"category": "GP", "lm_main": "GP04", "name_en": "Bis(monoacylglycero)phosphate", "name_zh": "双单酰甘油磷酸"},
    "LPC":    {"category": "GP", "lm_main": "GP01", "name_en": "Lyso-phosphatidylcholines", "name_zh": "溶血磷脂酰胆碱"},
    "LPE":    {"category": "GP", "lm_main": "GP02", "name_en": "Lyso-phosphatidylethanolamines", "name_zh": "溶血磷脂酰乙醇胺"},
    "LPS":    {"category": "GP", "lm_main": "GP03", "name_en": "Lyso-phosphatidylserines", "name_zh": "溶血磷脂酰丝氨酸"},
    "LPG":    {"category": "GP", "lm_main": "GP04", "name_en": "Lyso-phosphatidylglycerols", "name_zh": "溶血磷脂酰甘油"},
    "LPI":    {"category": "GP", "lm_main": "GP06", "name_en": "Lyso-phosphatidylinositols", "name_zh": "溶血磷脂酰肌醇"},
    "LPA":    {"category": "GP", "lm_main": "GP10", "name_en": "Lyso-phosphatidic acids", "name_zh": "溶血磷脂酸"},
    # ---- Sphingolipids ----
    "SPB":    {"category": "SP", "lm_main": "SP01", "name_en": "Sphingoid bases", "name_zh": "鞘氨醇碱基"},
    "Cer":    {"category": "SP", "lm_main": "SP02", "name_en": "Ceramides", "name_zh": "神经酰胺"},
    "CerP":   {"category": "SP", "lm_main": "SP02", "name_en": "Ceramide-1-phosphates", "name_zh": "神经酰胺磷酸"},
    "SM":     {"category": "SP", "lm_main": "SP03", "name_en": "Sphingomyelins", "name_zh": "鞘磷脂"},
    "HexCer": {"category": "SP", "lm_main": "SP05", "name_en": "Hexosylceramides", "name_zh": "己糖基神经酰胺"},
    "LacCer": {"category": "SP", "lm_main": "SP05", "name_en": "Lactosylceramides", "name_zh": "乳糖基神经酰胺"},
    "GM3":    {"category": "SP", "lm_main": "SP06", "name_en": "Gangliosides", "name_zh": "神经节苷脂"},
    "SL":     {"category": "SP", "lm_main": "SP01", "name_en": "Sphingolipids", "name_zh": "鞘脂"},
    # ---- Sterol lipids ----
    "ST":     {"category": "ST", "lm_main": "ST01", "name_en": "Sterols", "name_zh": "固醇"},
    "CE":     {"category": "ST", "lm_main": "ST0102", "name_en": "Cholesteryl esters", "name_zh": "胆固醇酯"},
    "Chol":   {"category": "ST", "lm_main": "ST0101", "name_en": "Cholesterol", "name_zh": "胆固醇"},
    "BA":     {"category": "ST", "lm_main": "ST04", "name_en": "Bile acids", "name_zh": "胆汁酸"},
    # ---- Prenols ----
    "CoQ":    {"category": "PR", "lm_main": "PR02", "name_en": "Ubiquinones", "name_zh": "泛醌"},
    "VitE":   {"category": "PR", "lm_main": "PR02", "name_en": "Tocopherols", "name_zh": "生育酚"},
}

# Common alternative spellings normalized to the canonical key above.
# 常见的别名，统一到上面的规范键。
_SYNONYMS = {
    "TAG": "TG", "DAG": "DG", "MAG": "MG",
    "PIP2": "PIP", "PIP3": "PIP", "PI4P": "PIP",
    "LYSOPC": "LPC", "LYSOPE": "LPE", "LYSOPI": "LPI", "LYSOPS": "LPS", "LYSOPA": "LPA", "LYSOPG": "LPG",
    "CERAMIDE": "Cer", "SPHINGOMYELIN": "SM", "CHOLESTERYL ESTER": "CE",
    "SO": "SPB", "SA": "SPB", "SPH": "SPB", "DHCER": "Cer", "DHSM": "SM",
    "GLCCER": "HexCer", "GALCER": "HexCer", "HEXCER": "HexCer",
    "COH": "Chol", "CHOLESTEROL": "Chol", "FFA": "FA",
    "ACAR": "CAR", "ACYLCARNITINE": "CAR",
}


# --------------------------------------------------------------------------- #
# Canonical lipid metabolic pathways (curated from LIPID MAPS pathway maps and
# standard lipid biochemistry). Each names the subclasses that flow through it.
# 经典脂质代谢通路（依据 LIPID MAPS 通路图与标准脂质生化策展）。每条列出流经的亚类。
# --------------------------------------------------------------------------- #
PATHWAYS = [
    {
        "key": "fa_synthesis",
        "name_en": "Fatty-acid de novo synthesis & elongation",
        "name_zh": "脂肪酸从头合成与延长",
        "members": ["FA", "FAHFA"],
        "desc_en": "Acetyl-CoA -> malonyl-CoA -> palmitate (FASN), then ELOVL elongation and SCD/FADS desaturation.",
        "desc_zh": "乙酰辅酶A→丙二酰辅酶A→棕榈酸（FASN），再经 ELOVL 延长与 SCD/FADS 去饱和。",
        "lm": "LIPID MAPS Fatty Acyls (FA)",
    },
    {
        "key": "beta_oxidation",
        "name_en": "Fatty-acid beta-oxidation",
        "name_zh": "脂肪酸 β-氧化",
        "members": ["FA", "CAR", "AcCa"],
        "desc_en": "Carnitine shuttle (CPT1/CPT2) imports acyl chains for mitochondrial beta-oxidation.",
        "desc_zh": "肉碱穿梭（CPT1/CPT2）转运酰基链进入线粒体 β-氧化。",
        "lm": "LIPID MAPS Fatty Acyls (FA)",
    },
    {
        "key": "glycerolipid_tag",
        "name_en": "Glycerolipid / triacylglycerol synthesis",
        "name_zh": "甘油脂 / 三酰甘油合成",
        "members": ["PA", "DG", "TG", "MG", "LPA"],
        "desc_en": "G3P -> LPA -> PA -> DAG -> TAG (GPAT/AGPAT/PAP/DGAT); DAG is the branch point to phospholipids.",
        "desc_zh": "G3P→LPA→PA→DAG→TAG（GPAT/AGPAT/PAP/DGAT）；DAG 是通向磷脂的分叉点。",
        "lm": "LIPID MAPS Glycerolipids (GL)",
    },
    {
        "key": "kennedy",
        "name_en": "Kennedy pathway (CDP-choline / CDP-ethanolamine)",
        "name_zh": "Kennedy 通路（CDP-胆碱 / CDP-乙醇胺）",
        "members": ["PC", "PE", "DG"],
        "desc_en": "DAG + CDP-choline/CDP-ethanolamine -> PC / PE (CEPT1/CHPT1/EPT1).",
        "desc_zh": "DAG + CDP-胆碱/CDP-乙醇胺 → PC / PE（CEPT1/CHPT1/EPT1）。",
        "lm": "LIPID MAPS Glycerophospholipids (GP)",
    },
    {
        "key": "lands_cycle",
        "name_en": "Glycerophospholipid remodeling (Lands cycle)",
        "name_zh": "甘油磷脂重塑（Lands 循环）",
        "members": ["LPC", "LPE", "LPS", "LPI", "LPG", "LPA", "PC", "PE"],
        "desc_en": "Deacylation (PLA2) to lyso-phospholipids and re-acylation (LPCAT/MBOAT) reshape acyl composition.",
        "desc_zh": "PLA2 脱酰生成溶血磷脂，再经 LPCAT/MBOAT 再酰化，重塑酰基组成。",
        "lm": "LIPID MAPS Glycerophospholipids (GP)",
    },
    {
        "key": "ps_pe_pc",
        "name_en": "PS synthesis & PS -> PE -> PC interconversion",
        "name_zh": "PS 合成与 PS→PE→PC 互变",
        "members": ["PS", "PE", "PC"],
        "desc_en": "PSS1/2 make PS; PISD decarboxylates PS -> PE; PEMT methylates PE -> PC.",
        "desc_zh": "PSS1/2 合成 PS；PISD 脱羧 PS→PE；PEMT 甲基化 PE→PC。",
        "lm": "LIPID MAPS Glycerophospholipids (GP)",
    },
    {
        "key": "pi_signaling",
        "name_en": "Phosphatidylinositol signaling",
        "name_zh": "磷脂酰肌醇信号",
        "members": ["PI", "PIP", "LPI"],
        "desc_en": "PI is phosphorylated to PIP/PIP2/PIP3 (PI kinases) driving membrane signaling.",
        "desc_zh": "PI 被磷酸化为 PIP/PIP2/PIP3（PI 激酶），驱动膜信号转导。",
        "lm": "LIPID MAPS Glycerophospholipids (GP)",
    },
    {
        "key": "cardiolipin",
        "name_en": "Cardiolipin synthesis (mitochondrial)",
        "name_zh": "心磷脂合成（线粒体）",
        "members": ["PG", "PGP", "CL", "PA"],
        "desc_en": "PA -> CDP-DAG -> PGP -> PG -> CL (CRLS1); CL supports mitochondrial cristae.",
        "desc_zh": "PA→CDP-DAG→PGP→PG→CL（CRLS1）；心磷脂维持线粒体嵴结构。",
        "lm": "LIPID MAPS Glycerophospholipids (GP)",
    },
    {
        "key": "sphingolipid",
        "name_en": "Sphingolipid de novo synthesis & metabolism",
        "name_zh": "鞘脂从头合成与代谢",
        "members": ["SPB", "Cer", "CerP", "SM", "HexCer", "LacCer", "GM3", "SL"],
        "desc_en": "Serine + palmitoyl-CoA -> sphinganine -> ceramide (CERS) -> SM (SGMS) / glycosphingolipids.",
        "desc_zh": "丝氨酸 + 棕榈酰辅酶A→鞘氨醇→神经酰胺（CERS）→鞘磷脂（SGMS）/ 糖鞘脂。",
        "lm": "LIPID MAPS Sphingolipids (SP)",
    },
    {
        "key": "sterol_ester",
        "name_en": "Sterol & cholesteryl-ester metabolism",
        "name_zh": "固醇与胆固醇酯代谢",
        "members": ["ST", "Chol", "CE", "BA"],
        "desc_en": "Cholesterol esterification (ACAT/SOAT) to CE for storage; hydrolysis regenerates free sterol.",
        "desc_zh": "胆固醇经 ACAT/SOAT 酯化为 CE 储存；水解再生游离固醇。",
        "lm": "LIPID MAPS Sterol Lipids (ST)",
    },
]

_PATHWAY_BY_KEY = {p["key"]: p for p in PATHWAYS}


def lm_class_url(lm_main):
    """Build a stable LIPID MAPS structure-database URL for a main-class code.

    为主类代码拼接稳定的 LIPID MAPS 结构数据库 URL。
    """
    if not lm_main:
        return ""
    return f"https://www.lipidmaps.org/databases/lmsd/{lm_main}"


def _strip_modifiers(head):
    """Strip common lipidomics modifier decorations to reveal the base subclass.

    去掉常见的脂质修饰装饰，还原基础亚类。Handles ether/plasmalogen (O-/P-, -O), oxidation
    (Ox/O2x/+O/-OH) and methylation (DM) markers, e.g. 'PG-O' -> 'PG', 'O2xFA' -> 'FA', 'DMPE' -> 'PE'.
    """
    h = head.strip()
    up = h.upper()
    # Oxidized fatty acids: OxFA / O2xFA / O3xFA -> FA (oxylipins are still Fatty Acyls).
    # 氧化脂肪酸：OxFA/O2xFA/O3xFA -> FA（氧化脂仍属脂肪酰）。
    if up.endswith("XFA") or up in ("OXFA", "OFA"):
        return "FA"
    # Dimethyl-PE (PEMT intermediate) -> PE.
    if up == "DMPE":
        return "PE"
    base = h
    # Leading ether markers: 'O-PC', 'P-PE'.
    for pre in ("O-", "P-"):
        if base.upper().startswith(pre):
            base = base[len(pre):]
    # Trailing modifiers: ether '-O', oxidation '+O'/'+2O'/'-OH'/'-O'.
    for suf in ("-O2", "+2O", "-OH", "-O", "+O", "-P", "+2OH", "-2O"):
        if base.upper().endswith(suf):
            base = base[: -len(suf)]
            break
    return base.strip()


def normalize_subclass(value):
    """Map a raw subclass string to a canonical SUBCLASS_INFO key (or None).

    把原始亚类字符串映射到规范的 SUBCLASS_INFO 键（映射不到则返回 None）。
    """
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return None
    # Drop any acyl detail after a space/parenthesis, e.g. "PC 34:1" -> "PC".
    # 去掉空格/括号后的酰基细节，如 "PC 34:1" -> "PC"。
    head = s.replace("(", " ").split()[0].strip()

    def _lookup(token):
        for cand in (token, token.upper()):
            if cand in SUBCLASS_INFO:
                return cand
        up_t = token.upper()
        if up_t in _SYNONYMS:
            return _SYNONYMS[up_t]
        for key in SUBCLASS_INFO:
            if key.upper() == up_t:
                return key
        return None

    hit = _lookup(head)
    if hit:
        return hit
    # Retry after stripping ether/oxidation/methylation modifiers.
    # 去掉醚键/氧化/甲基化修饰后重试。
    stripped = _strip_modifiers(head)
    if stripped and stripped != head:
        return _lookup(stripped)
    return None


def classify_subclass(value):
    """Return the LIPID MAPS classification dict for a subclass, or None.

    返回某亚类的 LIPID MAPS 分类信息字典，映射不到则 None。
    """
    key = normalize_subclass(value)
    if key is None:
        return None
    info = dict(SUBCLASS_INFO[key])
    info["subclass"] = key
    cat_en, cat_zh = CATEGORIES.get(info["category"], (info["category"], info["category"]))
    info["category_name_en"] = cat_en
    info["category_name_zh"] = cat_zh
    info["lm_url"] = lm_class_url(info["lm_main"])
    info["pathways"] = [p["key"] for p in PATHWAYS if key in p["members"]]
    return info


def _hypergeom_p(n_hits, n_pathway, n_sig_total, n_total):
    """Over-representation p-value: P(X >= n_hits) for significant lipids in this pathway.

    过表达 p 值：该通路中显著脂质数 >= n_hits 的概率 P(X >= n_hits)。
    """
    if hypergeom is None or n_total <= 0 or n_pathway <= 0 or n_sig_total <= 0:
        return float("nan")
    # X ~ Hypergeom(M=n_total, n=n_sig_total, N=n_pathway); P(X>=n_hits) = sf(n_hits-1).
    return float(hypergeom.sf(n_hits - 1, n_total, n_sig_total, n_pathway))


def pathway_enrichment(diff_table):
    """Aggregate differential features per canonical pathway and score over-representation.

    按经典通路汇总差异特征并对富集进行打分。返回一个 DataFrame（每条通路一行）。
    """
    if diff_table is None or len(diff_table) == 0 or "subclass" not in diff_table.columns:
        return pd.DataFrame()

    df = diff_table.copy()
    df["_norm_subclass"] = df["subclass"].map(normalize_subclass)
    has_sig = "significant" in df.columns
    n_total = int(len(df))
    n_sig_total = int(df["significant"].sum()) if has_sig else 0

    rows = []
    for pathway in PATHWAYS:
        members = set(pathway["members"])
        sub = df[df["_norm_subclass"].isin(members)]
        if sub.empty:
            continue
        n = int(len(sub))
        sig = sub[sub["significant"]] if has_sig else sub.iloc[0:0]
        n_sig = int(len(sig))
        log2fc = pd.to_numeric(sub.get("log2fc"), errors="coerce").dropna()
        mean_log2fc = float(log2fc.mean()) if not log2fc.empty else float("nan")
        n_up = int((sig["direction"] == "up").sum()) if ("direction" in sig.columns and n_sig) else 0
        n_down = int((sig["direction"] == "down").sum()) if ("direction" in sig.columns and n_sig) else 0
        p_enrich = _hypergeom_p(n_sig, n, n_sig_total, n_total) if n_sig else float("nan")
        present = sorted({normalize_subclass(x) for x in sub["_norm_subclass"] if x})
        rows.append({
            "key": pathway["key"],
            "name_en": pathway["name_en"],
            "name_zh": pathway["name_zh"],
            "n_features": n,
            "n_sig": n_sig,
            "n_up": n_up,
            "n_down": n_down,
            "mean_log2fc": round(mean_log2fc, 3) if mean_log2fc == mean_log2fc else float("nan"),
            "enrich_p": p_enrich,
            "neg_log10_p": (-np.log10(p_enrich) if (p_enrich == p_enrich and p_enrich > 0) else float("nan")),
            "subclasses": ", ".join(present),
            "desc_en": pathway["desc_en"],
            "desc_zh": pathway["desc_zh"],
            "lm": pathway["lm"],
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values(
            by=["n_sig", "neg_log10_p", "n_features"], ascending=[False, False, False]
        ).reset_index(drop=True)
    return result


def pathway_plot(enrichment, output_dir, tag, language="en"):
    """Horizontal bar of pathway mean log2FC colored by enrichment significance. Returns PNG path.

    通路平均 log2FC 水平条形图，按富集显著性着色。返回 PNG 路径。
    """
    if enrichment is None or enrichment.empty:
        return None
    data = enrichment.copy()
    data = data[data["n_features"] >= 1]
    if data.empty:
        return None
    data = data.sort_values("mean_log2fc")
    labels = data["name_zh" if language == "zh" else "name_en"].tolist()
    values = data["mean_log2fc"].fillna(0.0).to_numpy()
    sig_color = data["neg_log10_p"].fillna(0.0).to_numpy()

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.5 * len(data) + 1.2)))
    norm = plt.Normalize(vmin=0, vmax=max(1.0, float(np.nanmax(sig_color)) if len(sig_color) else 1.0))
    colors = plt.cm.magma_r(norm(sig_color))
    bars = ax.barh(range(len(data)), values, color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Mean log2 fold change")
    title = "脂质代谢通路富集 (LIPID MAPS)" if language == "zh" else "Lipid metabolic pathway enrichment (LIPID MAPS)"
    ax.set_title(f"{title} ({tag})")
    for i, (bar, n_sig, n_feat) in enumerate(zip(bars, data["n_sig"], data["n_features"])):
        ax.text(bar.get_width(), i, f"  {int(n_sig)}/{int(n_feat)}", va="center",
                fontsize=8, color="#333")
    sm = plt.cm.ScalarMappable(cmap="magma_r", norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, label="-log10 enrichment p")
    fig.tight_layout()
    path = os.path.join(output_dir, f"pathway_enrichment_{tag}.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def category_summary(diff_table):
    """Aggregate differential features by LIPID MAPS top-level category. Returns a list of dicts.

    按 LIPID MAPS 顶层大类汇总差异特征。返回字典列表。
    """
    if diff_table is None or len(diff_table) == 0 or "subclass" not in diff_table.columns:
        return []
    df = diff_table.copy()
    df["_cat"] = df["subclass"].map(lambda v: (classify_subclass(v) or {}).get("category"))
    df = df[df["_cat"].notna()]
    out = []
    for cat, sub in df.groupby("_cat"):
        log2fc = pd.to_numeric(sub.get("log2fc"), errors="coerce").dropna()
        cat_en, cat_zh = CATEGORIES.get(cat, (cat, cat))
        out.append({
            "category": cat,
            "name_en": cat_en,
            "name_zh": cat_zh,
            "n_features": int(len(sub)),
            "n_sig": int(sub["significant"].sum()) if "significant" in sub.columns else 0,
            "mean_log2fc": round(float(log2fc.mean()), 3) if not log2fc.empty else float("nan"),
        })
    return sorted(out, key=lambda r: r["n_sig"], reverse=True)


def analyze_pathways(diff_table, output_dir, tag, language="en", log=print):
    """Full offline pathway analysis: enrichment table + category summary + plot.

    完整的离线通路分析：富集表 + 大类汇总 + 图。返回结构化字典供报告/结论使用。
    """
    os.makedirs(output_dir, exist_ok=True)
    enrichment = pathway_enrichment(diff_table)
    categories = category_summary(diff_table)
    plot = pathway_plot(enrichment, output_dir, tag, language=language)
    csv_path = None
    if not enrichment.empty:
        csv_path = os.path.join(output_dir, f"pathway_enrichment_{tag}.csv")
        enrichment.to_csv(csv_path, index=False, encoding="utf-8-sig")
    n_hit = 0 if enrichment.empty else int((enrichment["n_sig"] > 0).sum())
    log(f"Pathway enrichment: {len(enrichment)} pathways touched, {n_hit} with significant lipids "
        f"/ 通路富集：涉及 {len(enrichment)} 条通路，其中 {n_hit} 条含显著脂质")
    return {
        "tag": tag,
        "enrichment": enrichment.to_dict(orient="records") if not enrichment.empty else [],
        "categories": categories,
        "plot": plot,
        "csv": csv_path,
    }


def pathway_digest(pathway_result, top_n=6):
    """Compact text digest of the pathway result for the conclusion/LLM prompt.

    通路结果的紧凑文本摘要，供结论/大模型提示使用。
    """
    if not pathway_result:
        return ""
    lines = []
    cats = pathway_result.get("categories") or []
    if cats:
        lines.append("LIPID MAPS categories (name: n_sig/n, mean log2FC): " + "; ".join(
            f"{c['name_en']}: {c['n_sig']}/{c['n_features']}, {c['mean_log2fc']}" for c in cats[:6]))
    enrich = pathway_result.get("enrichment") or []
    hits = [e for e in enrich if e.get("n_sig", 0) > 0][:top_n]
    if hits:
        lines.append("Enriched metabolic pathways (pathway: n_sig/n, mean log2FC, p): " + "; ".join(
            f"{e['name_en']}: {e['n_sig']}/{e['n_features']}, {e['mean_log2fc']}, "
            f"{('%.1e' % e['enrich_p']) if e.get('enrich_p') == e.get('enrich_p') else 'NA'}"
            for e in hits))
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Optional online enrichment via the LIPID MAPS REST API (never blocks offline).
# 可选的在线增强：通过 LIPID MAPS REST API（绝不阻塞离线流程）。
# --------------------------------------------------------------------------- #
def online_enrich(subclasses, timeout=6, log=print):
    """Query the LIPID MAPS REST API for canonical records of the given subclasses.

    调用 LIPID MAPS REST API 获取给定亚类的规范记录。返回 {subclass: record} 字典；失败则空字典。
    Returns {subclass_key: {"lm_id":..., "name":..., "url":...}} or {} on any failure.
    """
    out = {}
    for raw in subclasses:
        key = normalize_subclass(raw)
        if not key or key in out:
            continue
        info = SUBCLASS_INFO.get(key)
        if not info:
            continue
        # LIPID MAPS REST: look up by abbreviation/name; endpoint is public and keyless.
        # LIPID MAPS REST：按缩写/名称查询；接口公开且无需密钥。
        url = (f"https://www.lipidmaps.org/rest/compound/abbrev/{key}/all/json")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LipidIN/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = json.loads(resp.read().decode("utf-8", "replace"))
            record = None
            if isinstance(payload, dict) and payload:
                first = next(iter(payload.values())) if not payload.get("lm_id") else payload
                record = first if isinstance(first, dict) else None
            if record:
                out[key] = {
                    "lm_id": record.get("lm_id") or record.get("regno") or "",
                    "name": record.get("name") or info["name_en"],
                    "url": lm_class_url(info["lm_main"]),
                }
        except Exception as e:
            log(f"LIPID MAPS online lookup skipped for {key}: {type(e).__name__} "
                f"/ LIPID MAPS 在线查询跳过 {key}")
            # One network failure means offline; stop hammering the API.
            # 一次联网失败即视为离线，停止继续请求。
            break
    return out


if __name__ == "__main__":
    import sys
    csv = sys.argv[1] if len(sys.argv) > 1 else None
    if csv and os.path.exists(csv):
        table = pd.read_csv(csv, encoding="utf-8-sig")
        result = analyze_pathways(table, os.path.dirname(csv) or ".", "cli", language="en")
        print(json.dumps({"categories": result["categories"],
                          "enrichment": result["enrichment"][:8]}, ensure_ascii=False, indent=2))
    else:
        # Self-test the mapping on a few subclasses.
        for s in ["PC 34:1", "TAG", "Cer", "CE", "LPC", "PI", "CL", "SM", "unknownX"]:
            print(s, "->", classify_subclass(s))
