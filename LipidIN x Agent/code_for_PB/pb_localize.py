# -*- coding: utf-8 -*-
"""
第三阶段：PB 双键位置重建（在 *_PB_diagnostic.csv 之上做后处理）。

两步（对每个 *_PB_diagnostic.csv）：
  步骤1 | 与 master_table_quant.csv 对同一物质（master 只做“细化”，不做“闸门”）：
     - PB 侧的双键定位由 PB 诊断离子自证，**不依赖 master 是否收录该物质**：对不上
       master 的峰一律保留，只在 master_match_status 里记录对不上的原因。
     - 按“总组成(total_chain_name) 相同 + RT 接近”把 PB 峰对到 master 行；对上了才
       用 master 的信息去细化链划分。
     - 若 master 对该物质有明确链划分(sub_chain_name 里带 '_')，则把该 PB 峰的候选
       收窄到与 master 一致的链划分；若 master 只有总组成，则只在总组成层面认定
       （链划分不收窄，但仍可给出两条链各自的双键数/位置）。
     - 收窄若把候选清空（master 的链划分与 PB 诊断行给出的划分无交集），退回不收窄，
       而不是丢弃该峰——PB 的 omega 位置证据与链如何切分无关。

  步骤2 | 双键位置重建（“拼接”诊断行）：
     每条诊断行给出某条链上“某个双键”的 (n=离甲基端位置, f=其甲基侧还有几个双键,
     another_chain=另一条链的双键数)。一条有 D 个双键的链，需要 f=0..D-1 一组、n 严格递增
     的检出来完整定位；允许部分定位（能定几个定几个）。
     置信度 = 该定位所用诊断行 diag_pair_min_intensity 的最小值（最弱一环）。

  输出粒度（取决于 common 能否定链组成）：
     - **master 给了链划分**(sub_chain_name 带 '_'，来自 common 的 MS2 FA 碎片)：链组成
       已被 common 定死，链划分不再是冗余枚举 → 报**链级** 'PC 18:1(n-9)_18:2(n-6,n-9)'
       (collapse_to_chain)。`chain_resolved=True`。
     - **master 只有总组成**：PB 诊断离子定位的是双键“相对甲基端(omega)的位置 n-x”，与链
       如何切分无关，故 PE 15:3(n-3,n-6,n-9)_25:3(...) 与 16:3_24:3、17:3_23:3… 给出完全
       相同的 omega 位置与置信度，纯属链划分冗余。此时才把这些等价划分**合并成总链长:
       总不饱和度**一个标签(如 'PE 40:6')，双键位置用**去重后的 distinct omega 位置集合**
       表达(collapse_to_total)。`chain_resolved=False`。
     两条链在同一 n-x 都有双键时诊断离子质量相同、无法区分来自哪条链，故 distinct 位置
     集只诚实报“n-x 存在”(6 个双键可能只落在 5 个 distinct 位置上)。

  唯一输出 *_PB_isomers.csv（逐异构体一行，与 A1.xlsx 基准同构）：
     每个 PB 完整定位 = 一个竞争结构异构体行；common 层异构体 × PB 层异构体做笛卡尔积。
     - pb_isomer：该异构体标签。common 能定链组成时是链级 'PC 18:1(n-9)_18:2(n-6,n-9)'，
       否则退到总组成级 'PE 40:6(n-3,n-6,n-9,n-12,n-15)'；
     - pb_isomer_total_level：恒为总组成级标签，便于跨行汇总；
     - pb_double_bond_positions：该异构体的 omega 位置串 'n-3,n-6,n-9,...'；
     - all_double_bond_positions：该 (峰,物质) **所有**竞争异构体的位置并集；
     - chain_resolved：该行是否链级(True=common 给了链组成)；
     - pb_rank：该 (峰,物质) 内按 (位置数, 置信度) 的排名，1 = 最优定位；
     - n_isomers_this_peak：该 (峰,物质) 的竞争异构体总数(截断前)；
     - common_isomer / common_pct：common 侧异构体及占比(本 demo 无 MS2 链信息→单一总组成 100%)；
     - pb_pct：该定位的诊断离子总强度在竞争定位间的占比；
     - combined_pct = common_pct × pb_pct / 100：该完整结构在总峰里的占比；
     - isomer_confidence：该定位所用诊断行 diag_pair_min_intensity 的最小值(最弱一环)；
     - position_intensity_ratio：该异构体内部逐位置的特征离子强度占比(%)；
     - l1_coverage / l2_coverage：LEVEL1(母离子对)/LEVEL2(头基对)的碎片覆盖率(见 l12_coverage)；
     - matched_master_rt / master_evidence_type / master_match_status：与 common master 的对应
       情况，用于判断“定物质”这一步本身的可信度；未对上 master 时前两者为空。
"""

import os
import re
import glob
from itertools import product
from collections import defaultdict

import numpy as np
import pandas as pd


def log(msg):
    print(f"[PB-loc] {msg}")


# master_match_status 取值：master 只用于细化链划分，任何一种“对不上”都不丢弃该峰。
MASTER_MATCHED = "master匹配"              # 同物质 + RT 在容差内
MASTER_RT_OFF = "master有该物质但RT超容差"   # 同物质，但最近的 master RT 超出 rt_tol
MASTER_NO_SPECIES = "master无该物质"        # master 侧根本没有这个总组成
MASTER_PART_FALLBACK = "master匹配(链划分无交集,已退回不收窄)"


_CHAIN_TOKEN = re.compile(r"\d+:\d+(?:\([^)]*\))?")
# f 与 another_chain 都可选：
#   二酰基 diacyl：'15:3(n-3,f-0,another_chain-3)'（含 another_chain）
#   溶血 lyso    ：'20:4(n-15,f-3)'（单链无 another_chain）、'17:1(n-9)'（单不饱和连 f 都无）
_ANN_RE = re.compile(r"(\d+):(\d+)\(n-(\d+)(?:,f-(\d+))?(?:,another_chain-(\d+))?\)")
_CU_RE = re.compile(r"(\d+):(\d+)")
_NPOS_RE = re.compile(r"n-(\d+)")


def _split_class_body(name):
    m = re.match(r"^(\S+)\s+(.*)$", str(name).strip())
    if not m:
        return None, None
    return m.group(1), m.group(2)


def parse_pb_subchain(sub_chain_name):
    """解析 PB_diagnostic 的 sub_chain_name（一条链带 (n-,f-,another_chain-)）。

    返回 dict：cls, partition(排序后的 'c:u' 元组), chain('c:u'), n, f, another_chain。
    不含 (n-...) 的（L1/L2 总组成行）返回 None。
    """
    cls, body = _split_class_body(sub_chain_name)
    if cls is None:
        return None
    chains = _CHAIN_TOKEN.findall(body)
    annotated = None
    bare = []
    for tok in chains:
        m = _ANN_RE.search(tok)
        if m:
            annotated = {
                "chain": f"{m.group(1)}:{m.group(2)}",
                "n": int(m.group(3)),
                "f": int(m.group(4)) if m.group(4) is not None else 0,
                "another_chain": int(m.group(5)) if m.group(5) is not None else 0,
            }
        else:
            mm = _CU_RE.search(tok)
            if mm:
                bare.append(f"{mm.group(1)}:{mm.group(2)}")
    if annotated is None:
        return None
    partition = tuple(sorted([annotated["chain"]] + bare))
    return {
        "cls": cls,
        "partition": partition,
        "chain": annotated["chain"],
        "n": annotated["n"],
        "f": annotated["f"],
        "another_chain": annotated["another_chain"],
    }


def _chain_unsat(chain):
    m = _CU_RE.match(chain)
    return int(m.group(2)) if m else 0


def load_master_partitions(master_quant_csv):
    """master_table_quant.csv -> {total_chain_name: [(rt, partition_or_None, evidence), ...]}。

    partition：sub_chain_name 里有 '_' 时的排序 'c:u' 元组，否则 None（只有总组成）。
    evidence：该 master 注释的证据等级(evidence_type，如 'MS2注释'/'MS1推断(MS1_only)')。
    """
    df = pd.read_csv(master_quant_csv, encoding="utf-8-sig")
    out = defaultdict(list)
    for _, r in df.iterrows():
        total = str(r.get("total_chain_name", "")).strip()
        if not total or total.lower() == "nan":
            continue
        sub = str(r.get("sub_chain_name", "")).strip()
        rt = pd.to_numeric(r.get("sample_rt"), errors="coerce")
        evidence = str(r.get("evidence_type", "")).strip()
        partition = None
        if "_" in sub or "/" in sub:
            _, body = _split_class_body(sub)
            if body:
                toks = [f"{m.group(1)}:{m.group(2)}"
                        for m in _CU_RE.finditer(body.replace("/", "_"))]
                if toks:
                    partition = tuple(sorted(toks))
        out[total].append((float(rt) if pd.notna(rt) else None, partition, evidence))
    return dict(out)


def load_l12_ions(pb_library_pkl):
    """PB 库 -> {total_chain_name: {'1': [L1 离子 m/z], '2': [L2 离子 m/z]}}。

    LEVEL1 = 母离子对 (pmz, pmz−H2O)；LEVEL2 = 头基对：胆碱类(PC/PC-O/LPC/LPC-O)
    是游离头基 184.07/104.11，其余类(PE/PE-O/PI/PG/LPE…)沿用 (pmz, pmz2)——**故非胆碱
    类的 L2 与 L1 是同一对离子**，"L1 与 L2 都够"对它们等价于单一条件。
    """
    import pickle
    with open(pb_library_pkl, "rb") as f:
        lib = pickle.load(f)
    lab = lib["lab_mz"].astype(str)
    level = lab.str.extract(r"^LOL_(\d+)_DTL_")[0]
    keep = level.isin(["1", "2"])
    after = lab[keep].str.split("_DTL_", n=1).str[1]
    total = after.str.split(r"_CAH_|_ADD_|_FML_", n=1, regex=True).str[0]
    frame = pd.DataFrame({
        "total": total,
        "level": level[keep],
        "mz": pd.to_numeric(lib.loc[keep, "mz_value"], errors="coerce"),
    }).dropna(subset=["mz"])
    out = {}
    for (t, lv), g in frame.groupby(["total", "level"]):
        out.setdefault(t, {})[lv] = sorted(set(g["mz"].round(4)))
    del lib
    return out


def _parse_joined(text):
    try:
        return [float(x) for x in str(text).split("_") if x not in ("", "nan")]
    except ValueError:
        return []


def l12_coverage(total, sample_ms2_mz, l12_map, ppm_err2=20.0, da_err2=0.02):
    """该峰对 total 的 LEVEL1/LEVEL2 碎片覆盖率 (命中根数/该 level 应有根数)。

    每个 level 各有 2 根碎片，命中 1 根 = 0.5。返回 (l1_cov, l2_cov)；
    库里查不到该物质时返回 (nan, nan)。
    """
    ions = l12_map.get(total) if l12_map else None
    if not ions:
        return float("nan"), float("nan")
    mz = np.asarray(sample_ms2_mz, dtype=float)

    def cov(level_key):
        targets = ions.get(level_key)
        if not targets:
            return float("nan")
        hit = sum(1 for m in targets
                  if np.any(np.abs(mz - m) <= max(m * ppm_err2 * 1e-6, da_err2)))
        return hit / len(targets)

    return cov("1"), cov("2")


def _enumerate_chain_ladders(dets):
    """某条链的检出 -> 该链所有自洽双键定位。

    dets: [(f, n, conf), ...]。按 f 分组，每个被检出的 f 取一个 n、n 随 f 严格递增
    （允许部分：只用被检出的 f 层）。返回 [(positions_tuple_sorted, conf_min), ...]。
    """
    byf = defaultdict(dict)
    for f, n, c in dets:
        if n not in byf[f] or c > byf[f][n]:
            byf[f][n] = c
    fs = sorted(byf)
    results = []

    def rec(i, prev_n, chosen, cmin):
        if i == len(fs):
            if chosen:
                results.append((tuple(sorted(chosen)), cmin))
            return
        for n, c in sorted(byf[fs[i]].items()):
            if n > prev_n:
                rec(i + 1, n, chosen + [n], min(cmin, c))

    rec(0, -1, [], float("inf"))
    return results


def _fmt_chain(chain, positions):
    """('18:2', (6,9)) -> '18:2(n-6,n-9)'；无定位 -> '18:2'。"""
    if not positions:
        return chain
    return f"{chain}(" + ",".join(f"n-{p}" for p in positions) + ")"


def assemble_peak(peak_rows, allowed_partitions=None):
    """把一个 PB 峰的诊断行拼接成分子级双键定位列表。

    peak_rows: list of dict(parse_pb_subchain 结果 + 'conf')。
    allowed_partitions: 若非 None，只处理这些链划分（master 有明确链划分时）。
    返回按 (定位双键数, 置信度) 降序的 [(sub_chain_str, n_db, conf), ...]。
    """
    by_part = defaultdict(list)
    for row in peak_rows:
        part = row["partition"]
        if allowed_partitions is not None and part not in allowed_partitions:
            continue
        # another_chain 一致性：另一条链的双键数须等于 another_chain
        others = [c for c in part if c != row["chain"]]
        # 同名双链时 others 可能为空，跳过该一致性校验
        if others and _chain_unsat(others[0]) != row["another_chain"]:
            continue
        by_part[part].append(row)

    localizations = []
    for part, rows in by_part.items():
        cls = rows[0]["cls"]
        # 每条“链槽”收集检出；同名双链共享同一检出池（近似处理）
        chain_dets = defaultdict(list)
        for r in rows:
            chain_dets[r["chain"]].append((r["f"], r["n"], r["conf"]))

        # 对 partition 里每个不同链值枚举梯子；饱和/无检出链 -> 单一“无定位”
        slot_options = []
        for chain in part:
            dets = chain_dets.get(chain)
            if dets and _chain_unsat(chain) > 0:
                ladders = _enumerate_chain_ladders(dets)
                if not ladders:
                    ladders = [((), float("inf"))]
            else:
                ladders = [((), float("inf"))]
            slot_options.append([(chain, pos, conf) for pos, conf in ladders])

        for combo in product(*slot_options):
            n_db = sum(len(pos) for _, pos, _ in combo)
            if n_db == 0:
                continue
            confs = [conf for _, _, conf in combo if conf != float("inf")]
            conf = min(confs) if confs else 0.0
            chain_strs = [_fmt_chain(ch, pos) for ch, pos, _ in combo]
            sub = f"{cls} " + "_".join(chain_strs)
            localizations.append((sub, n_db, conf))

    localizations.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return localizations


def collapse_to_total(localizations):
    """把 assemble_peak 的“按链划分”定位合并到总组成层面。

    localizations: [(sub_chain_str, n_db, conf), ...]（sub_chain_str 里含各链 n-x 位置）。
    等价链划分(15:3_25:3 / 16:3_24:3 …)给出相同 omega 位置，合并成同一条。
    双键位置取**跨两条链去重后的 distinct omega 集合**（同一 n-x 出现在两条链上
    无法从诊断离子质量区分，只报一次）。同一位置集合取最高置信度。
    返回 [(positions_tuple_sorted, conf), ...]，按 (定位到的 distinct 位置数, conf) 降序。
    """
    best = {}
    for sub, _n_db, conf in localizations:
        positions = tuple(sorted({int(x) for x in _NPOS_RE.findall(sub)}))
        if not positions:
            continue
        if positions not in best or conf > best[positions]:
            best[positions] = conf
    collapsed = list(best.items())
    collapsed.sort(key=lambda pc: (len(pc[0]), pc[1]), reverse=True)
    return collapsed


def collapse_to_chain(localizations):
    """**保留链划分**的定位：仅在 master(common) 已给出该物质链组成时使用。

    与 collapse_to_total 相反——不把 18:1_18:2 压成 36:3。common 的 MS2 已经把链组成
    定死，链划分不再是冗余枚举，此时报 'PC 18:1(n-9)_18:2(n-6,n-9)' 才是应有粒度。
    返回 [(sub_chain_str, distinct_positions_tuple, conf), ...]，按 (位置数, conf) 降序。
    """
    best = {}
    for sub, _n_db, conf in localizations:
        if sub not in best or conf > best[sub]:
            best[sub] = conf
    out = []
    for sub, conf in best.items():
        positions = tuple(sorted({int(x) for x in _NPOS_RE.findall(sub)}))
        if positions:
            out.append((sub, positions, conf))
    out.sort(key=lambda x: (len(x[1]), x[2]), reverse=True)
    return out


def common_isomers_for_peak(total_chain_name):
    """common 层的异构体及各自占比(%)，返回 [(isomer_label, pct), ...]。

    当前实现：本批 PB 物质在 common 里只有单一总组成（无 MS2、无链异构体，见数据现实），
    故恒返回 [(总组成, 100.0)]。

    扩展点（保持接口不变即可自动生效）：若将来某峰在 common 有 MS2 解析出的多条链
    组成异构体及其 FA 碎片相对强度，则在此返回 [(链组成A, %A), (链组成B, %B), ...]
    （各 % 由“该异构体 FA 链特征碎片强度 / 全谱”归一而来）。下游 combined_pct 会自动
    做 common × PB 双层展开：combined_pct = common_pct × pb_pct / 100。
    """
    return [(total_chain_name, 100.0)]


def localize_one_file(pb_diag_csv, master_map, rt_tol=0.5, isomer_top_n=20,
                      l12_map=None, l12_min_coverage=0.5, ppm_err2=20.0, da_err2=0.02):
    try:
        df = pd.read_csv(pb_diag_csv, encoding="utf-8-sig")
    except pd.errors.EmptyDataError:
        # PB 搜库无命中时会写出空的 diagnostic CSV（无列）——正常情况，跳过。
        return pd.DataFrame(), pd.DataFrame()
    if df.empty or "has_double_bond_position" not in df.columns:
        return pd.DataFrame(), pd.DataFrame()

    d = df[df["has_double_bond_position"] == True].copy()
    if d.empty:
        return pd.DataFrame(), pd.DataFrame()

    d["_parsed"] = d["sub_chain_name"].map(parse_pb_subchain)
    d = d[d["_parsed"].notna()].copy()
    d["_conf"] = pd.to_numeric(d.get("diag_pair_min_intensity"), errors="coerce").fillna(0.0)
    d["_tot"] = pd.to_numeric(d.get("diag_pair_total_intensity"), errors="coerce").fillna(0.0)

    iso_records = []
    dbpos_records = []
    n_l12_rejected = 0
    # 必须按 (peak_id, total_chain_name) 分组：一个峰常同时命中多个物质（实测 PB1_A1
    # 419 个带双键证据的峰里 183 个命中 ≥2 个物质）。只按 peak_id 分组会(1)每峰只留下
    # iloc[0] 那一个物质、静默丢弃其余，(2)把别的物质的诊断行混进这个物质的拼接里。
    for (peak_id, total_raw), g in d.groupby(["peak_id", "total_chain_name"]):
        total = str(total_raw).strip()
        peak_rt = float(pd.to_numeric(g["sample_rt"].iloc[0], errors="coerce"))

        # 步骤1：总组成 + RT 接近，对到 master。master 只用于细化链划分，对不上不丢弃。
        master_rows = master_map.get(total, [])
        near = [(rt, part, ev) for (rt, part, ev) in master_rows
                if rt is not None and abs(rt - peak_rt) <= rt_tol]
        if near:
            rt_matched, _p, ev_matched = min(near, key=lambda x: abs(x[0] - peak_rt))
            matched_rt = rt_matched
            matched_evidence = ev_matched
            explicit = {part for _, part, _ in near if part is not None}
            allowed = explicit if explicit else None  # None=只在总组成层面，不收窄链划分
            master_status = MASTER_MATCHED
        else:
            # PB 的 omega 位置由诊断离子自证，master 缺该物质/RT 对不上都不影响该证据成立。
            matched_rt = None
            matched_evidence = ""
            explicit = set()
            allowed = None
            master_status = MASTER_NO_SPECIES if not master_rows else MASTER_RT_OFF

        # common 里没有该物质时，物质本身的正确性只能由 PB 自证：要求 LEVEL1(母离子对)
        # 与 LEVEL2(头基对)各自至少命中 1/2 根碎片。对上了 master 的峰有 common 背书，
        # 不走这道门。链划分仍然给不出——只报总组成 + omega 位置，这是本管线的既定粒度。
        l1_cov, l2_cov = float("nan"), float("nan")
        if l12_map is not None:
            l1_cov, l2_cov = l12_coverage(
                total, _parse_joined(g["sample_ms2_mz"].iloc[0]),
                l12_map, ppm_err2=ppm_err2, da_err2=da_err2,
            )
            if master_status in (MASTER_NO_SPECIES, MASTER_RT_OFF):
                if not (l1_cov >= l12_min_coverage and l2_cov >= l12_min_coverage):
                    n_l12_rejected += 1
                    continue

        # 每个 omega 位置 n-x 的“特征离子总强度”：取该峰所有 n=x 诊断行的最大总强度
        pos_intensity = defaultdict(float)
        for _, r in g.iterrows():
            n_pos = r["_parsed"]["n"]
            pos_intensity[n_pos] = max(pos_intensity[n_pos], float(r["_tot"]))

        peak_rows = []
        for _, r in g.iterrows():
            pr = dict(r["_parsed"])
            pr["conf"] = float(r["_conf"])
            pr["tot"] = float(r["_tot"])
            peak_rows.append(pr)

        locs = assemble_peak(peak_rows, allowed_partitions=allowed)
        partition_constrained = allowed is not None
        if not locs and allowed is not None:
            # master 的链划分与本峰诊断行给出的划分无交集。omega 位置与链如何切分无关，
            # 故退回“不收窄”重拼，而不是把整峰的双键证据丢掉。
            locs = assemble_peak(peak_rows, allowed_partitions=None)
            partition_constrained = False
            if locs:
                master_status = MASTER_PART_FALLBACK
        if not locs:
            continue

        def _fmt_total(positions):
            return f"{total}(" + ",".join(f"n-{p}" for p in positions) + ")"

        # 粒度取决于 common 能否定链组成：
        #   master 给了链划分 -> 报链级 'PC 18:1(n-9)_18:2(n-6,n-9)'；
        #   master 只有总组成 -> 链划分是纯冗余枚举，压成 'PC 36:3(n-6,n-9)'。
        if partition_constrained:
            collapsed = collapse_to_chain(locs)
        else:
            collapsed = [(_fmt_total(pos), pos, conf)
                         for pos, conf in collapse_to_total(locs)]
        if not collapsed:
            continue

        # 该 (峰,物质) 所有竞争异构体的位置并集——单条异构体行看不到全貌，这列给出。
        union_pos = sorted({p for _lb, pos, _c in collapsed for p in pos})

        # ---- 链×单双键位置表（最细粒度，定量用）----
        # 每条诊断行 = 某条链上某一个双键的一对特征离子(L3+L4)，其 diag_pair_total_intensity
        # 就是“确定这个双键的那对特征离子”的强度。按 (链, n) 归并(同 n 多个 f 取最强的那对)，
        # 再在该 (峰,物质) 内归一得到占比 -> 下游 pb_quantify 用 物质总定量 × 占比 定量到
        # “某条链的某个双键位置”。
        # 链归属只在 master 给了链划分(partition_constrained)时才成立；否则链划分是冗余
        # 枚举，chain 记为空，粒度退到 (物质, n)。
        # 键 = (f, n)，**不含链**：诊断离子质量 = pmz − C{x}H{2x−2f}O，只由 (n,f) 与母离子
        # 决定，**与链如何切分无关**。库里同一对离子会被复制到每一种链划分上，按链记录
        # 等于把同一个测量重复计数(实测 PC 42:8 的 18:2/20:4/22:4/24:6 四条链给出完全
        # 相同的强度)。链归属在离子层面从来不是实测量，故这里只到 (物质, f, n)。
        # f 不能并掉——它决定了行与行的关系：
        #   同 f 不同 n = 该链“第 f+1 个双键”落在不同位置 = **真正竞争的位置异构体**
        #                (如 18:1 的 n-9 油酸 vs n-7 异油酸)；
        #   不同 f      = 同一个分子上的不同双键(如 20:4 的 n-6/n-9/n-12/n-15)。
        pair_int = {}   # (f, n) -> 该对特征离子的 total 强度
        for pr in peak_rows:
            if partition_constrained and pr["partition"] not in allowed:
                continue
            key = (pr["f"], pr["n"])
            if key not in pair_int or pr["tot"] > pair_int[key]:
                pair_int[key] = pr["tot"]
        # 同一 f 内归一：该 f 层每个分子有且只有一个双键，故各 n 的份额是对这批分子的
        # 真实划分，**天然加总为 1**。跨 f 归一没有这个含义(会把同一批分子数多次)。
        f_sum = defaultdict(float)
        for (f_val, _n), tot_int in pair_int.items():
            f_sum[f_val] += tot_int
        for (f_val, n_pos), tot_int in sorted(pair_int.items(), key=lambda kv: (kv[0][0], -kv[1])):
            dbpos_records.append({
                "peak_id": int(peak_id),
                "sample_pmz": g["sample_pmz"].iloc[0],
                "sample_rt": peak_rt,
                "subclass": total.split()[0],
                "total_chain_name": total,
                "double_bond_position": f"n-{n_pos}",
                "f": f_val,
                "pair_total_intensity": tot_int,
                # 该 f 层内部的位置份额：同 f 各行之和 = 1，即该层双键的位置分布。
                "share_within_f": (tot_int / f_sum[f_val]) if f_sum[f_val] > 0 else float("nan"),
                # 该对离子在本峰所有特征离子对里的占比（跨 f），用于把物质总量拆到每个双键。
                "pair_share": (tot_int / sum(pair_int.values())) if pair_int else float("nan"),
                "master_match_status": master_status,
                "master_evidence_type": matched_evidence,
            })

        # ---- 异构体反卷积表：common 层 × PB 完整定位层，逐异构体一行 ----
        # 每个 PB 完整定位(collapsed 里一条 omega 位置组合)= 一个竞争结构异构体；
        # pb_pct = 该定位的诊断离子总强度(各位置特征离子总强度之和)在竞争定位间归一。
        loc_intensity = [sum(float(pos_intensity.get(p, 0.0)) for p in pos)
                         for _lb, pos, _c in collapsed]
        total_loc_int = sum(loc_intensity)
        common_isos = common_isomers_for_peak(total)
        for rank, ((iso_label, pos, conf), li) in enumerate(
                zip(collapsed[:isomer_top_n], loc_intensity[:isomer_top_n]), 1):
            pb_pct = (100.0 * li / total_loc_int) if total_loc_int > 0 else float("nan")
            # 该异构体内部逐位置占比
            iso_pos_int = {p: float(pos_intensity.get(p, 0.0)) for p in pos}
            s = sum(iso_pos_int.values())
            iso_pos_ratio = "; ".join(
                f"n-{p}:{(100.0 * iso_pos_int[p] / s) if s > 0 else float('nan'):.1f}%"
                for p in pos
            )
            for common_iso, common_pct in common_isos:
                iso_records.append({
                    "peak_id": int(peak_id),
                    "sample_pmz": g["sample_pmz"].iloc[0],
                    "sample_rt": peak_rt,
                    "subclass": total.split()[0],
                    "total_chain_name": total,
                    "common_isomer": common_iso,
                    "common_pct": round(common_pct, 4),
                    "pb_isomer": iso_label,
                    "pb_isomer_total_level": _fmt_total(pos),
                    "pb_double_bond_positions": ",".join(f"n-{p}" for p in pos),
                    "chain_resolved": partition_constrained,
                    "pb_pct": round(pb_pct, 4),
                    "combined_pct": round(common_pct * pb_pct / 100.0, 4),
                    "pb_rank": rank,                    # 1 = 该(峰,物质)置信度最高的定位
                    "n_isomers_this_peak": len(collapsed),
                    "all_double_bond_positions": ",".join(f"n-{p}" for p in union_pos),
                    "n_positions": len(pos),
                    "isomer_confidence": conf,
                    "isomer_diag_total_intensity": li,
                    "position_intensity_ratio": iso_pos_ratio,
                    "l1_coverage": l1_cov,
                    "l2_coverage": l2_cov,
                    "master_matched_rt": matched_rt,
                    "master_evidence_type": matched_evidence,
                    "master_match_status": master_status,
                })

    if n_l12_rejected:
        log(f"  {os.path.basename(pb_diag_csv)}: common 无背书且 L1/L2 碎片覆盖 <"
            f"{l12_min_coverage} 的峰剔除 {n_l12_rejected} 个")
    return pd.DataFrame(iso_records), pd.DataFrame(dbpos_records)


def localize_all_pb_diagnostic(pb_output_dir, master_quant_csv,
                               rt_tol=0.5, isomer_top_n=20,
                               pb_library_pkl=None, l12_min_coverage=0.5,
                               ppm_err2=20.0, da_err2=0.02,
                               isomer_suffix="_PB_isomers.csv",
                               dbpos_suffix="_PB_dbpos.csv"):
    master_map = load_master_partitions(master_quant_csv)
    log(f"master 物质数 {len(master_map)}；RT 容差 ±{rt_tol} min")

    l12_map = None
    if pb_library_pkl:
        l12_map = load_l12_ions(pb_library_pkl)
        log(f"已载入 PB 库 LEVEL1/2 碎片：{len(l12_map)} 个物质；common 无背书的峰需"
            f"L1/L2 覆盖率各 ≥{l12_min_coverage}")
    else:
        log("未提供 pb_library_pkl：跳过 L1/L2 碎片覆盖判据，common 无背书的峰全部保留")

    files = sorted(glob.glob(os.path.join(pb_output_dir, "*_PB_diagnostic.csv")))
    summary = []
    for path in files:
        iso, dbpos = localize_one_file(path, master_map, rt_tol=rt_tol,
                                       isomer_top_n=isomer_top_n,
                                       l12_map=l12_map, l12_min_coverage=l12_min_coverage,
                                       ppm_err2=ppm_err2, da_err2=da_err2)
        iso.to_csv(path.replace("_PB_diagnostic.csv", isomer_suffix),
                   index=False, encoding="utf-8-sig")
        dbpos.to_csv(path.replace("_PB_diagnostic.csv", dbpos_suffix),
                     index=False, encoding="utf-8-sig")
        n_peaks = iso.groupby(["peak_id", "total_chain_name"]).ngroups if len(iso) else 0
        summary.append({"file": os.path.basename(path), "peak_species": n_peaks,
                        "isomers": len(iso), "dbpos": len(dbpos)})
        log(f"{os.path.basename(path)} -> (峰,物质) {n_peaks}，异构体行 {len(iso)}，"
            f"链×双键行 {len(dbpos)}")
    total_ps = sum(s["peak_species"] for s in summary)
    total_iso = sum(s["isomers"] for s in summary)
    total_db = sum(s["dbpos"] for s in summary)
    log(f"完成：{len(files)} 个文件，共 {total_ps} 个(峰,物质)给出双键定位，"
        f"{total_iso} 条异构体行(*{isomer_suffix})，{total_db} 条链×双键行(*{dbpos_suffix})。")
    return summary


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass
    localize_all_pb_diagnostic(
        pb_output_dir=r"D:\bio_inf\LipidIN_github\PB_demo\sample\PB\pkl",
        master_quant_csv=r"D:\bio_inf\LipidIN_github\PB_demo\sample\common\pkl\master_table_quant.csv",
        pb_library_pkl=r"D:\bio_inf\LipidIN_github\PB_demo\PB library\PB_PX_H.pkl",
        rt_tol=0.5,
    )
