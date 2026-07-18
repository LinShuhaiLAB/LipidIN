# -*- coding: utf-8 -*-
"""
第四阶段：PB 异构体级定量 + 总量归一。

复用 quantify.py 的 EIC 机制（load_ms1_scans / build_eic_matrix / integrate_peak），
把定量靶从 common 主表换成 PB 侧的物质峰。

定量单元 = **(物质, f, 单个双键位置)**：
  1) 物质级：PB 衍生化母离子的 MS1 EIC 峰面积——该物质的总量，由 MS1 测得，
     与 MS2 是否采到无关，故每个样本都有值。
  2) 链×双键级：用 *_PB_dbpos.csv 里的 pair_share 把物质总量拆开——
       area(物质, f, n-x) = species_area × pair_share(f, n-x)
     pair_share = “确定该双键的那对特征离子(L3+L4)”的强度在该峰所有特征离子对里的占比，
     各行之和 = 1，故各行面积之和 = 该物质的 EIC 总面积。

  **不做链归属**：诊断离子质量 = pmz − C{x}H{2x−2f}O，只由 (n,f) 与母离子决定，与链
  如何切分无关，库里同一对离子会被复制到每种链划分上。链级信息在鉴定层(*_PB_isomers.csv
  的 pb_isomer)保留，但定量层只到 (物质, f, n-x)。

  share_within_f：同一 f 层内部各位置的份额，**同 f 各行之和 = 1**。该 f 层每个分子有且
  只有一个双键，故这是对这批分子的真实划分——“A 组偏向 n-9 / B 组偏向 n-7”只在此列成立。

  ★ 关键前提与偏倚（务必看）：pair_share 来自诊断离子强度，**不等于摩尔丰度**。
    PB 裂解效率与双键位置相关——越靠甲基端(f=0，如 n-3)越易裂解、诊断离子越强，实测
    可占 76–89%。未做裂解偏倚校正时，本表会系统性高估靠甲基端的双键位置。
    本模块按用户要求只做总量归一，不做裂解偏倚校正；该偏倚原样传递到输出。

跨样本对齐：
  各样本的 *_PB_dbpos.csv 是各自独立的峰表。本模块按 (total_chain_name, chain,
  double_bond_position) 分组，组内对 sample_rt 做贪心聚类(RT_CLUSTER_TOL)，每个簇 =
  一个跨样本一致的行。同一位置在不同 RT 出现(真的色谱异构体)会分成不同簇，不强行合并。

归一（用户选定“总量归一”）：
  sample_total(S) = 该样本**所有定量靶的物质级 EIC 面积之和**（不是异构体面积之和——
  物质级面积每样本都有，比依赖 MS2 检出的异构体面积稳定）。
  normalized = isomer_area / sample_total(S) × 1e6   （即 ppm of total quantified lipid）

输出 PB_dbposition_quant.csv：一行一个 (物质, 链, 双键位置)，样本列成对给出
  <sample>__area（该链该双键的面积）与 <sample>__norm（总量归一后，×1e6）。
  某样本的 MS2 没检出该链该双键时 pair_share 未知 -> 该样本两列为空(NaN)，
  而不是填 0——“没测到”与“测到且为零”不是一回事。
"""

import os
import re
import glob
import time
import concurrent.futures
from multiprocessing import cpu_count

import numpy as np
import pandas as pd

from quantify import log, load_ms1_scans, build_eic_matrix, integrate_peak

# 同一 (物质, 异构体) 的 sample_rt 相差在此之内视为同一色谱峰（跨样本对齐）。
# PB 衍生化后 RT 在样本间会漂移，实测 common/PB 之间可达 0.34 min，样本间略小。
RT_CLUSTER_TOL = 0.3

# 定量靶去重用的分箱精度（与 quantify.quantify_master_table 一致）
MZ_ROUND = 4
RT_ROUND = 2

NORM_SCALE = 1e6  # 归一后乘的系数：ppm of total quantified lipid


# PB1_A1 与 PB2_A1 是**同一个样品打的两针**（technical injections），不是两个批次、更不是
# 两个生物样本。故在定量层必须合并成一个样品 A1：两针的 MS2 证据取并集(每针只被 DDA 采到
# ~320/1197 个单元，合并后覆盖显著提高)，EIC 面积取两针均值。
# 注：进样级的重复性(针与针之间)在 pb_qc_diff 的质控里单独评估，那里仍看未合并的 PB1/PB2。
INJECTION_RE = re.compile(r"^PB(\d+)_(.+)$")


def split_injection(run_name):
    """PB1_A1 -> ('A1', 'PB1')；不符合该命名的(如 PB_MSMS1)原样返回、无进样号。"""
    m = INJECTION_RE.match(run_name)
    if m:
        return m.group(2), "PB" + m.group(1)
    return run_name, ""


def pool_dbpos_tables(pb_output_dir, dbpos_suffix="_PB_dbpos.csv"):
    """把各次进样的 *_PB_dbpos.csv 汇成一张长表，加 run / sample / injection 列。"""
    files = sorted(glob.glob(os.path.join(pb_output_dir, "*" + dbpos_suffix)))
    if not files:
        raise FileNotFoundError(f"{pb_output_dir} 下没有 *{dbpos_suffix}，请先跑第三阶段")
    frames = []
    for path in files:
        df = pd.read_csv(path, encoding="utf-8-sig")
        if df.empty:
            continue
        # PB1_A1_MS1_MSMS_PB_dbpos.csv -> run=PB1_A1（与 .raw 同名），再拆出 sample=A1
        run = os.path.basename(path)[: -len(dbpos_suffix)]
        if run.endswith("_MS1_MSMS"):
            run = run[: -len("_MS1_MSMS")]
        df["run"] = run
        df["sample"], df["injection"] = split_injection(run)
        frames.append(df)
    pooled = pd.concat(frames, ignore_index=True)
    log(f"汇总 {len(files)} 次进样的链×双键表，共 {len(pooled)} 行；"
        f"进样 {pooled['run'].nunique()} 次 -> 样品 {pooled['sample'].nunique()} 个")
    return pooled


def cluster_rt(values, tol=RT_CLUSTER_TOL):
    """对一维 RT 做贪心聚类，返回每个元素的簇号。相邻值差 > tol 就开新簇。"""
    order = np.argsort(values, kind="mergesort")
    labels = np.empty(len(values), dtype=int)
    cid = 0
    anchor = None
    for i in order:
        v = values[i]
        if anchor is None or v - anchor > tol:
            cid += 1
            anchor = v
        labels[i] = cid
    return labels


def build_consensus_dbpos(pooled, rt_tol=RT_CLUSTER_TOL):
    """跨样本一致的 (物质峰, 链, f, 双键位置) 表。

    **必须先把“物质峰”聚出来，再往峰里放位置**：先按 (total_chain_name) 分组、对 RT 做
    贪心聚类得到物质峰(species_cluster)，峰内再分 (链, f, 位置)。
    不能反过来按 (物质,链,f,位置) 各自聚 RT——那样同一个峰上竞争的 n-6 与 n-9 会各自
    算出略不同的 RT 中位数、拿到不同的 species_key，于是永远分不进同一组，下游“同一物质
    内部各位置谁占优”的比较就全成了单元素组(实测 1878 组里 1860 组只剩 1 个位置)。
    """
    pooled = pooled.copy()
    pooled["_rt"] = pd.to_numeric(pooled["sample_rt"], errors="coerce")
    pooled["_mz"] = pd.to_numeric(pooled["sample_pmz"], errors="coerce")
    pooled = pooled.dropna(subset=["_rt", "_mz"])

    # 第一层：物质峰 = (总组成) 内部按 RT 聚类
    sp_cluster = np.zeros(len(pooled), dtype=int)
    for _, idx in pooled.groupby("total_chain_name", sort=False).indices.items():
        sp_cluster[idx] = cluster_rt(pooled["_rt"].to_numpy()[idx], tol=rt_tol)
    pooled["species_cluster"] = sp_cluster
    sp_info = pooled.groupby(["total_chain_name", "species_cluster"], sort=False).agg(
        species_mz=("_mz", "median"), species_rt=("_rt", "median")).reset_index()
    sp_info["species_key"] = (sp_info["total_chain_name"] + "@"
                              + sp_info["species_rt"].map(lambda x: f"{x:.2f}"))
    pooled = pooled.merge(sp_info, on=["total_chain_name", "species_cluster"], how="left")

    # 第二层：峰内按 (链, f, 位置)。f 不能并掉——同链同 f 不同 n 才是竞争的位置异构体。
    group_keys = ["species_key", "total_chain_name", "f", "double_bond_position"]
    agg = pooled.groupby(group_keys, sort=False).agg(
        subclass=("subclass", "first"),
        consensus_mz=("species_mz", "first"),
        consensus_rt=("species_rt", "first"),
        n_samples_ms2=("sample", "nunique"),
        median_pair_share=("pair_share", "median"),
        median_share_within_f=("share_within_f", "median"),
        max_pair_intensity=("pair_total_intensity", "max"),
        master_match_status=("master_match_status", "first"),
        master_evidence_type=("master_evidence_type", "first"),
    ).reset_index()
    log(f"物质峰 {len(sp_info)} 个；峰内 (f,双键位置) 单元 {len(agg)} 个"
        f"（来自 {len(pooled)} 条样本级行）")
    return pooled, agg


def _quantify_one_raw(task_args):
    raw_path, target_mz_sorted, target_rt_sorted = task_args
    sample = os.path.splitext(os.path.basename(raw_path))[0]
    try:
        start = time.time()
        # cache_dir=None：*_MS1_centroid.npz 缓存可能已被清理，直接读 .raw（慢但正确）
        retention_times, scans = load_ms1_scans(raw_path, cache_dir=None)
        eic = build_eic_matrix(target_mz_sorted, scans)
        areas = np.array([
            integrate_peak(retention_times, eic[i], target_rt_sorted[i])
            for i in range(len(target_mz_sorted))
        ])
        log(f"{sample} 定量完成，MS1扫描 {len(scans)}，非零靶 {(areas > 0).sum()}/{len(areas)}，"
            f"耗时 {time.time() - start:.1f}s")
        return sample, areas
    except Exception as e:
        log(f"{sample} 定量失败 {type(e).__name__}: {e}")
        return sample, np.full(len(target_mz_sorted), np.nan)


def quantify_pb_dbpositions(pb_output_dir, pb_raw_dir, output_csv=None,
                            max_workers=None, rt_tol=RT_CLUSTER_TOL,
                            dbpos_suffix="_PB_dbpos.csv"):
    if output_csv is None:
        output_csv = os.path.join(pb_output_dir, "PB_dbposition_quant.csv")

    log("Step 1 | 汇总各样本链×双键表")
    pooled = pool_dbpos_tables(pb_output_dir, dbpos_suffix)

    log("Step 2 | 跨样本对齐 (物质,链,双键位置)（组内按 RT 贪心聚类）")
    pooled, cons = build_consensus_dbpos(pooled, rt_tol=rt_tol)

    log("Step 3 | 构建物质级定量靶（PB 母离子 EIC）")
    # 一个 species_key = 一个物质峰 = 一个 EIC 定量靶
    cons["_mz_key"] = cons["consensus_mz"].round(MZ_ROUND)
    cons["_rt_key"] = cons["consensus_rt"].round(RT_ROUND)
    targets = cons[["species_key", "_mz_key", "_rt_key"]].drop_duplicates(subset=["species_key"])
    targets = targets.sort_values("_mz_key", kind="mergesort").reset_index(drop=True)
    target_mz = targets["_mz_key"].to_numpy(dtype=float)
    target_rt = targets["_rt_key"].to_numpy(dtype=float)
    log(f"唯一定量靶 {len(targets)}（来自 {len(cons)} 个 链×双键 单元）")

    raw_files = sorted(glob.glob(os.path.join(pb_raw_dir, "*.raw")))
    if not raw_files:
        raise FileNotFoundError(f"未在 {pb_raw_dir} 找到 .raw")
    if max_workers is None:
        max_workers = max(1, int(cpu_count() * 0.85))
    max_workers = max(1, min(int(max_workers), len(raw_files)))

    log(f"Step 4 | 对 {len(raw_files)} 个 raw(进样) 提 EIC 并积分，{max_workers} 进程")
    start = time.time()
    run_areas = {}
    args_list = [(p, target_mz, target_rt) for p in raw_files]
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        for run, areas in ex.map(_quantify_one_raw, args_list):
            run_areas[run] = areas
    log(f"EIC 定量总耗时 {time.time() - start:.1f}s")

    # 同一样品的多针取均值 -> 每个样品一列（PB1_A1 与 PB2_A1 是同一个 A1 打了两次）
    by_sample = {}
    for run, areas in run_areas.items():
        by_sample.setdefault(split_injection(run)[0], []).append(areas)
    for smp, arrs in sorted(by_sample.items()):
        targets[smp] = np.nanmean(np.vstack(arrs), axis=0)
    samples = sorted(by_sample)
    n_multi = sum(1 for v in by_sample.values() if len(v) > 1)
    log(f"进样合并：{len(run_areas)} 次进样 -> {len(samples)} 个样品"
        f"（其中 {n_multi} 个样品有多针、EIC 面积取均值）")

    # 物质级总量单独出一张表：这是“常规脂质组学也能拿到”的那一层，下游差异分析要用它
    # 跟位置组成层做对照，才能挑出“总量没变但双键位置分布变了”这种 PB 独有的发现。
    species_csv = os.path.join(os.path.dirname(output_csv), "PB_species_quant.csv")
    species = (cons[["species_key", "subclass", "total_chain_name", "consensus_mz",
                     "consensus_rt"]]
               .drop_duplicates(subset=["species_key"])
               .merge(targets.drop(columns=["_mz_key", "_rt_key"]), on="species_key", how="left"))
    sp_tot = {s: float(np.nansum(targets[s].to_numpy())) for s in samples}
    for s in samples:
        species[f"{s}__norm"] = (species[s] / sp_tot[s] * NORM_SCALE) if sp_tot[s] > 0 else np.nan
    species = species.rename(columns={s: f"{s}__area" for s in samples})
    species.to_csv(species_csv, index=False, encoding="utf-8-sig")
    log(f"物质级定量表已保存 {species_csv}（{len(species)} 个物质峰）")

    log("Step 5 | 物质总量 × 该样本的 pair_share -> 链×双键级面积")
    cons = cons.merge(targets.drop(columns=["_mz_key", "_rt_key"]), on="species_key", how="left")

    # 多针合并后必须**从强度重算 share**，不能拿某一针的 pair_share 了事——每针的
    # pair_share 是在该针内部归一的，直接取其中一针(或取max)既丢掉另一针的证据，也让
    # 各行之和不再为 1。这里先把两针的同一 (样品,物质峰,f,位置) 强度取最强，再归一。
    key_cols = ["sample", "species_key", "f", "double_bond_position"]
    inten = (pooled.groupby(key_cols, sort=False)["pair_total_intensity"]
             .max().reset_index())
    # pair_share：在该 (样品, 物质峰) 的所有特征离子对之间归一 -> 各行之和 = 1，
    # 故各行面积之和 = 该物质的 EIC 总面积。
    inten["pair_share"] = (inten["pair_total_intensity"]
                           / inten.groupby(["sample", "species_key"])["pair_total_intensity"]
                           .transform("sum"))
    share_lookup = {
        (r.sample, r.species_key, r.f, r.double_bond_position): r.pair_share
        for r in inten.itertuples(index=False)
    }

    # 样本总量 = 该样本所有定量靶的物质级面积之和（每样本都有，不依赖 MS2 检出）
    sample_total = {s: float(np.nansum(targets[s].to_numpy())) for s in samples}
    for s in samples:
        if sample_total[s] <= 0:
            log(f"[警告] 样本 {s} 的物质级面积总和为 0，其归一列将全为 NaN")

    area_cols, norm_cols = [], []
    for s in samples:
        species_area = cons[s].to_numpy(dtype=float)
        shares = np.array([
            share_lookup.get(
                (s, r.species_key, r.f, r.double_bond_position), np.nan)
            for r in cons.itertuples(index=False)
        ], dtype=float)
        db_area = species_area * shares     # share 为 NaN(该样本 MS2 没检出) -> NaN
        cons[f"{s}__area"] = db_area
        tot = sample_total[s]
        cons[f"{s}__norm"] = (db_area / tot * NORM_SCALE) if tot > 0 else np.nan
        area_cols.append(f"{s}__area")
        norm_cols.append(f"{s}__norm")

    cons["n_samples_quantified"] = (~cons[area_cols].isna() & (cons[area_cols] > 0)).sum(axis=1)
    cons = cons.drop(columns=["_mz_key", "_rt_key"] + samples)

    front = ["species_key", "subclass", "total_chain_name", "double_bond_position",
             "f", "consensus_mz", "consensus_rt", "n_samples_ms2", "n_samples_quantified",
             "median_pair_share", "median_share_within_f", "max_pair_intensity",
             "master_match_status", "master_evidence_type"]
    cons = cons[front + area_cols + norm_cols]
    cons = cons.sort_values(["subclass", "total_chain_name", "consensus_rt", "f"],
                            kind="mergesort").reset_index(drop=True)
    cons.to_csv(output_csv, index=False, encoding="utf-8-sig")

    log(f"链×双键定量表已保存 {output_csv}")
    log(f"单元 {len(cons)}（物质×f×双键位置），样本 {len(samples)}"
        f"（每样本 __area 面积 / __norm 总量归一×{NORM_SCALE:.0g}）")
    log(f"至少在 1 个样本定量到的 {int((cons['n_samples_quantified'] > 0).sum())}；"
        f"全部 {len(samples)} 个样本都定量到的 "
        f"{int((cons['n_samples_quantified'] == len(samples)).sum())}")
    log("[提醒] pair_share 源自诊断离子强度，靠甲基端(f=0,如 n-3)的双键被系统性高估"
        "(实测占 76-89%)；本表未做裂解偏倚校正。")
    return cons


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass
    quantify_pb_isomers(
        pb_output_dir=r"D:\bio_inf\LipidIN_github\PB_demo\sample\PB\pkl",
        pb_raw_dir=r"D:\bio_inf\LipidIN_github\PB_demo\sample\PB",
    )
