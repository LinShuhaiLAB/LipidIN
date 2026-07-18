# -*- coding: utf-8 -*-
"""
第二阶段：PB 衍生化检索。与 common 完全独立——直接对全 PB 库搜库，不用主表裁剪。

流程：
  1) raw2pkl：sample/PB 下的 .raw -> *_MS1_MSMS.pkl。
  2) pb_diagnostic_search：用全 PB 库跑 EQ3，判据是“PB 诊断离子成对存在”
     (strip_dtl_prefix 后 L1≡L2、L3≡L4，同 label 命中≥2 根离子才保留)，
     输出 *_PB_diagnostic.csv：每 (峰, 库条目) 一行，含 has_double_bond_position
     与 diag_pair_min/total_intensity。该表是第三阶段 pb_localize 的输入。

母离子容差注意：这批 raw 的 DDA 母离子被量化到 0.1 网格，tm_raw2pkl 的精修只能在
隔离窗内取最强 MS1 峰去猜，同一化合物的精修母离子在扫描间可飘 -144~+84 ppm。故
ppm_err1/da_err1 对 PB 要放宽(run.py 用 60ppm/0.05Da)，让母离子只做粗筛，真正的
鉴别交给高精度诊断离子对(ppm_err2/da_err2)。详见 run.py 里 pb_ppm_err1 的注释。
"""
import os
import re
import glob
import pickle
import numpy as np
import pandas as pd


# ---------- 与 master_table.py 对齐的常量 ----------
EVIDENCE_MS2 = "MS2注释"
CHECK_MISSING = "缺失"

# master_table 里链“c:u”的正则(取整数碳:整数不饱和度)
_CU_RE = re.compile(r"(\d+)\s*:\s*(\d+)")


def log(msg):
    print(f"[PB] {msg}")


def _from_excel_text(value):
    """还原被写成 ="..." 的 Excel 文本(与 master_table 输出兼容)。"""
    if pd.isna(value):
        return ""
    text = str(value)
    if text.startswith('="') and text.endswith('"'):
        return text[2:-1]
    return text


def _parse_total_cu(total_chain_name):
    """'PC 34:1' -> ('PC', 34, 1)；解析不出返回 (class, None, None)。"""
    text = str(total_chain_name).strip()
    if not text or text.lower() == "nan":
        return "", None, None
    parts = text.split(maxsplit=1)
    cls = parts[0]
    if len(parts) < 2:
        return cls, None, None
    m = _CU_RE.search(parts[1])
    if not m:
        return cls, None, None
    return cls, int(m.group(1)), int(m.group(2))


def _chain_multiset(fa_composition):
    """'16:0_18:1' 或 'O-16:0/18:1' -> frozenset 计数元组 (((c,u),n),...)；空则 None。"""
    text = _from_excel_text(fa_composition).replace("/", "_")
    if not text.strip():
        return None
    chains = []
    for token in text.split("_"):
        m = _CU_RE.search(token)          # 忽略可能的 'O-' 醚前缀，只取 c:u
        if m:
            chains.append((int(m.group(1)), int(m.group(2))))
    if not chains:
        return None
    # 用排序元组表示多重集，便于做集合成员判断
    return tuple(sorted(chains))


# ============================================================
# 1) 从主表构建“每峰一标签”的过滤表
# ============================================================
# ============================================================
# 背景：常规注释流程在 postprocess 里用强度加权 Jaccard(JCR>0.2) 过滤匹配。
# PB 的双键诊断离子对(LEVEL3/4)在样品二级谱里天然低丰度(谱被头基/母离子主导)，
# 整谱强度相似度结构上过不了 0.2，导致所有 (n-x) 注释被删。
# 这里为 PB 单独实现“诊断离子存在性”判据：
#   - 复用同一 C++ 搜索引擎；
#   - 某个双键假设(某链、某 n- 位)的 L3、L4 两个诊断离子在容差内**都**被检出，
#     即认定该双键位置成立(去 DTL 前缀后 L3/L4 折叠成同一标签，靠 repeated_mask 表达“成对命中”)；
#   - 额外给每个 (n-) 行算“诊断离子对最小强度”作为置信分，并支持一个强度下限过滤。
# 不改动 EQ2_batchsearchv2.py，第一阶段常规注释流程完全不受影响。

_LOL_RE = re.compile(r"^LOL_(\d+)_DTL_")


def _parse_joined_spectrum(mz_text, int_text):
    """把下划线连接的 mz / intensity 串还原成两个 numpy 数组。"""
    def to_arr(text):
        vals = []
        for tok in str(text).split("_"):
            tok = tok.strip()
            if tok and tok.lower() != "nan":
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
        return np.asarray(vals, dtype=float)
    return to_arr(mz_text), to_arr(int_text)


def _build_diag_ion_map(library_df):
    """去 DTL 前缀后的 LEVEL3/4 标签 -> 该双键假设的诊断离子 m/z 数组。"""
    lab = library_df["lab_mz"].astype(str)
    level = lab.str.extract(_LOL_RE)[0]
    stripped = lab.str.split("_DTL_", n=1).str[1]
    is34 = level.isin(["3", "4"]).to_numpy()
    s_arr = stripped.to_numpy()
    mz_arr = library_df["mz_value"].to_numpy(dtype=float)
    diag = {}
    for i in np.nonzero(is34)[0]:
        diag.setdefault(s_arr[i], []).append(mz_arr[i])
    return {k: np.asarray(v, dtype=float) for k, v in diag.items()}


def _diag_pair_intensities(stripped_label, diag_map, sample_mz, sample_int,
                           ppm_err2, da_err2):
    """双键诊断离子对在样品谱里的 (最小强度, 总强度)。

    最小强度=对里较弱的一根(作置信分/存在性下限)；
    总强度=L3(较高)+L4(较低)两根特征离子各自最强匹配之和(即“该双键位特征离子总强度”)。
    任一离子缺失则两者都返回 NaN。
    """
    dmz = diag_map.get(stripped_label)
    if dmz is None or dmz.size == 0 or sample_mz.size == 0:
        return np.nan, np.nan
    matched = []
    for m0 in dmz:
        tol = max(m0 * ppm_err2 * 1e-6, da_err2)
        hit = np.where(np.abs(sample_mz - m0) <= tol)[0]
        if hit.size == 0:
            return np.nan, np.nan
        matched.append(float(sample_int[hit[np.argmax(sample_int[hit])]]))
    if not matched:
        return np.nan, np.nan
    return float(min(matched)), float(sum(matched))


def pb_diagnostic_search(
    filtered_lib, sample_inputs, cpp_file,
    ppm_err1, da_err1, ppm_err2, da_err2,
    diag_intensity_floor=0.0, output_suffix="_PB_diagnostic.csv",
):
    """对 PB 样品用“诊断离子存在性”判据搜库并逐样品写表。返回每样品统计列表。"""
    import EQ2_batchsearchv2 as E

    if not E.compile_cpp_module(cpp_file):
        raise RuntimeError("无法编译/加载 C++ 搜索模块")
    if not E.compile_cpp_similarity_module("eq3_similarity.cpp"):
        raise RuntimeError("无法编译/加载 C++ 相似度模块")

    search_engine, library_df, _sim = E.load_library_once(filtered_lib)
    diag_map = _build_diag_ion_map(library_df)
    del library_df

    summary = []
    for sample_input in sample_inputs:
        with open(sample_input, "rb") as f:
            data = sorted(pickle.load(f), key=lambda x: x[1])

        peak_ids, pmz, rt, ms2, inten = E.prepare_sample_data(data)
        matches = E.run_single_search(
            search_engine, peak_ids, pmz, rt, ms2, inten,
            ppm_err1, da_err1, ppm_err2, da_err2,
        )
        raw = E.convert_matches_to_dataframe(matches, search_engine)
        del matches

        out_path = sample_input.replace(".pkl", "") + output_suffix
        if raw is None or raw.empty:
            pd.DataFrame().to_csv(out_path, index=False, encoding="utf-8-sig")
            summary.append({"sample": os.path.basename(sample_input), "rows": 0, "n_db": 0})
            continue

        # 去 DTL 前缀 -> L1≡L2、L3≡L4 折叠；repeated_mask=“成对命中”即存在性判据
        raw["label"] = raw["label"].map(E.strip_dtl_prefix)
        subclass0 = raw["label"].astype(str).str.split(" ", n=1).str[0]
        repeated = raw.duplicated(subset=["peak_id", "label"], keep=False)
        kept = raw.loc[repeated | subclass0.isin(E.SINGLE_PEAK_SUBCLASSES)].copy()
        kept.drop_duplicates(subset=["peak_id", "label"], keep="first", inplace=True)
        kept.reset_index(drop=True, inplace=True)

        has_db = kept["label"].astype(str).str.contains(r"\(n-", regex=True)

        # 逐 (n-) 行算诊断离子对最小强度 + 总强度
        diag_min = np.full(len(kept), np.nan)
        diag_tot = np.full(len(kept), np.nan)
        for idx in np.nonzero(has_db.to_numpy())[0]:
            smz, sint = _parse_joined_spectrum(
                kept["sample_ms2_mz"].iat[idx], kept["sample_ms2_int"].iat[idx]
            )
            diag_min[idx], diag_tot[idx] = _diag_pair_intensities(
                kept["label"].iat[idx], diag_map, smz, sint, ppm_err2, da_err2
            )
        kept["diag_pair_min_intensity"] = diag_min
        kept["diag_pair_total_intensity"] = diag_tot

        # 可选强度下限：只对 (n-) 行生效
        if diag_intensity_floor > 0:
            drop = has_db.to_numpy() & ~(
                pd.to_numeric(kept["diag_pair_min_intensity"], errors="coerce")
                .fillna(-1.0).to_numpy() >= diag_intensity_floor
            )
            kept = kept.loc[~drop].copy().reset_index(drop=True)
            has_db = kept["label"].astype(str).str.contains(r"\(n-", regex=True)

        out = E.split_label_columns_before_output(kept)      # 拆出 total/sub_chain_name(含 n-)
        out = out.copy()
        out["has_double_bond_position"] = has_db.to_numpy()
        out["diag_pair_min_intensity"] = kept["diag_pair_min_intensity"].to_numpy()
        out["diag_pair_total_intensity"] = kept["diag_pair_total_intensity"].to_numpy()
        out.to_csv(out_path, index=False, encoding="utf-8-sig")

        n_db = int(has_db.sum())
        n_peak_db = int(out.loc[has_db, "peak_id"].nunique()) if n_db else 0
        summary.append({
            "sample": os.path.basename(sample_input),
            "rows": len(out), "n_db": n_db, "peaks_with_db": n_peak_db,
        })
        log(f"{os.path.basename(sample_input)}: 总行 {len(out)}，含双键位置(n-) 行 {n_db}"
            f"（涉及 {n_peak_db} 个峰）-> {os.path.basename(out_path)}")

    return summary


# ============================================================
# 3) 驱动：raw2pkl(PB) + 过滤 + PB 诊断离子搜库
# ============================================================
def run_pb_search(
    pb_raw_dir,
    pb_library_pkl,
    cpp_file,
    ms2_filter=0.0005,
    ppm_err1=10, da_err1=0.01, ppm_err2=20, da_err2=0.05,
    eq3_workers=1, raw_workers=None,
    pb_output_dir=None,
    diag_intensity_floor=0.0,
):
    """PB 检索：直接对全 PB 库搜库，与 common 完全独立。

    PB 结果与 common 的对应关系由第三阶段 pb_localize 独立处理
    (total_chain_name + rt 匹配)，不在搜库阶段耦合——master 只用于细化链划分，
    不做闸门。曾有过 use_master_filter 分支拿主表当库白名单裁剪 PB 库，默认关闭、
    从未被 run.py 启用，只留下过一个把人误导的陈旧 PB_PX_H_filtered.pkl，已删除。"""
    from tm_raw2pkl import ThermoRawToPklAgent
    from EQ2 import find_ms1_msms_files

    if pb_output_dir is None:
        pb_output_dir = os.path.join(pb_raw_dir, "pkl")

    log("Step A | PB raw -> pkl")
    # 需显式 run()，构造函数只保存参数不执行。
    ThermoRawToPklAgent(
        folder_path=pb_raw_dir, MS2_filter=ms2_filter, max_workers=raw_workers,
        recursive=0, output_folder=pb_output_dir,
    ).run()

    log("Step B | 全 PB 库搜库(PB 诊断离子存在性判据，保留 n- 双键位置)")
    sample_inputs = find_ms1_msms_files(pb_output_dir)
    if not sample_inputs:
        raise FileNotFoundError(f"{pb_output_dir} 下没有 *MS1_MSMS.pkl，raw2pkl 可能失败")
    results = pb_diagnostic_search(
        filtered_lib=pb_library_pkl, sample_inputs=sample_inputs, cpp_file=cpp_file,
        ppm_err1=ppm_err1, da_err1=da_err1, ppm_err2=ppm_err2, da_err2=da_err2,
        diag_intensity_floor=diag_intensity_floor,
    )
    total_db = sum(s.get("n_db", 0) for s in results)
    log(f"PB 搜库完成，样本数 {len(sample_inputs)}；共 {total_db} 条带双键位置(n-)注释；结果表暂不做后处理。")
    return results
