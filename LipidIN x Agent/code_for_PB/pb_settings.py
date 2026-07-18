# -*- coding: utf-8 -*-
"""所有可配置项的**单一事实来源**：一张表驱动校验、大模型词表、写回 run.py、以及界面展示。

为什么要这张表：参数散在 run.py、pb_qc_diff.py、agent 会话里三处，用户说"把 MS2 容差
放宽到 0.05"时，得知道这句话指的是哪个变量、合法范围是多少、改完谁受影响。把这些集中
成数据（而不是散落在 if/else 里），加一个参数就只改一行，大模型的可用词表也自动跟着更新。

`run` 字段 = 该项在 run.py 里的全局名（None 表示只活在 agent 会话里，如 alpha/假设条数）。
`stages` 字段 = 改了它之后**哪些阶段的产物过期了**，据此提醒用户需要重跑什么——参数改了
却不重跑，界面上的数字就是骗人的。
"""

import os

# kind: dir(必须已存在) / dir_out(不存在就建) / file(必须已存在) / float / int / bool
SETTINGS = {
    # ---------------- 路径 ----------------
    "project_root": {
        "kind": "dir", "run": "PROJECT_ROOT", "group": "path",
        "zh": "项目根目录（按标准布局重建 sample/common、sample/PB、PB library）",
        "en": "Project root (rebuilds sample/common, sample/PB, PB library)",
        "stages": ["annotation", "search", "localize", "quantify"]},
    "raw_root": {
        "kind": "dir", "run": "RAW_ROOT", "group": "path",
        "zh": "原始文件根目录（其下 common/ 与 PB/ 两个子目录，自动分派并各自转 pkl）",
        "en": "Raw root (with common/ and PB/ subfolders; auto-assigned and converted to pkl)",
        "stages": ["annotation", "search", "localize", "quantify"]},
    "common_raw_dir": {
        "kind": "dir", "run": "COMMON_RAW_DIR", "group": "path",
        "zh": "常规注释的 raw 目录", "en": "Common raw folder",
        "stages": ["annotation"]},
    "pb_raw_dir": {
        "kind": "dir", "run": "PB_RAW_DIR", "group": "path",
        "zh": "PB 衍生化的 raw 目录", "en": "PB raw folder",
        "stages": ["search", "quantify"]},
    "library_dir": {
        "kind": "dir", "run": "LIBRARY_DIR", "group": "path",
        "zh": "谱库目录", "en": "Library folder", "stages": ["annotation", "search"]},
    "common_library": {
        "kind": "file", "run": "COMMON_LIBRARY", "group": "path",
        "zh": "常规谱库 pkl（负模式 CH3COO）", "en": "Common library pkl",
        "stages": ["annotation"]},
    "pb_library": {
        "kind": "file", "run": "PB_LIBRARY", "group": "path",
        "zh": "PB 谱库 pkl（正模式 [M+H]+）", "en": "PB library pkl",
        "stages": ["search", "localize"]},
    "common_out_dir": {
        "kind": "dir_out", "run": "COMMON_OUT_DIR", "group": "path",
        "zh": "常规注释输出目录", "en": "Common output folder", "stages": ["annotation"]},
    "pb_out_dir": {
        "kind": "dir_out", "run": "PB_OUT_DIR", "group": "path",
        "zh": "PB 输出目录（所有 PB_*.csv 与报告都在这）", "en": "PB output folder",
        "stages": ["search", "localize", "quantify", "qc"]},

    # ---------------- 搜库/注释参数 ----------------
    "ms2_filter": {
        "kind": "float", "run": "MS2_filter", "group": "param", "min": 0, "max": 1,
        "zh": "低丰度二级碎片过滤阈值", "en": "MS2 low-abundance filter",
        "stages": ["annotation", "search"]},
    "ppm_err1": {
        "kind": "float", "run": "ppm_err1_input", "group": "param", "min": 0, "max": 200,
        "zh": "常规 MS1 ppm 容差", "en": "Common MS1 ppm", "stages": ["annotation"]},
    "da_err1": {
        "kind": "float", "run": "da_err1_input", "group": "param", "min": 0, "max": 1,
        "zh": "常规 MS1 Da 容差", "en": "Common MS1 Da", "stages": ["annotation"]},
    "ppm_err2": {
        "kind": "float", "run": "ppm_err2_input", "group": "param", "min": 0, "max": 200,
        "zh": "MS2 ppm 容差（common 与 PB 共用）", "en": "MS2 ppm",
        "stages": ["annotation", "search", "localize"]},
    "da_err2": {
        "kind": "float", "run": "da_err2_input", "group": "param", "min": 0, "max": 1,
        "zh": "常规 MS2 Da 容差", "en": "Common MS2 Da", "stages": ["annotation"]},
    "pb_ppm_err1": {
        "kind": "float", "run": "pb_ppm_err1", "group": "param", "min": 0, "max": 300,
        "zh": "PB 母离子 ppm 容差（故意比 common 松：raw 的 DDA 母离子被量化到 0.1 网格）",
        "en": "PB precursor ppm (deliberately loose)", "stages": ["search"]},
    "pb_da_err1": {
        "kind": "float", "run": "pb_da_err1", "group": "param", "min": 0, "max": 1,
        "zh": "PB 母离子 Da 容差", "en": "PB precursor Da", "stages": ["search"]},
    "pb_da_err2": {
        "kind": "float", "run": "pb_da_err2", "group": "param", "min": 0, "max": 1,
        "zh": "PB 诊断碎片 Da 容差（故意比 common 紧：诊断离子近噪声）",
        "en": "PB diagnostic fragment Da (deliberately tight)",
        "stages": ["search", "localize"]},
    "rt_tol": {
        "kind": "float", "run": "PB_LOC_RT_TOL", "group": "param", "min": 0, "max": 5,
        "zh": "PB 峰与 master 匹配的 RT 容差(min)", "en": "PB-to-master RT tolerance (min)",
        "stages": ["localize"]},
    "l12_min_coverage": {
        "kind": "float", "run": "PB_LOC_L12_MIN_COV", "group": "param", "min": 0, "max": 1,
        "zh": "无 common 背书时 LEVEL1/2 的最低碎片覆盖率", "en": "L1/L2 min fragment coverage",
        "stages": ["localize"]},
    "isomer_top_n": {
        "kind": "int", "run": "PB_LOC_ISOMER_TOP_N", "group": "param", "min": 1, "max": 200,
        "zh": "每个(峰,物质)最多写出多少竞争异构体", "en": "Max competing isomers per peak",
        "stages": ["localize"]},
    "rt_cluster_tol": {
        "kind": "float", "run": "PB_QUANT_RT_CLUSTER_TOL", "group": "param", "min": 0, "max": 5,
        "zh": "跨样本对齐同一异构体的 RT 聚类容差(min)", "en": "RT cluster tolerance (min)",
        "stages": ["quantify"]},
    "master_ms2_only": {
        "kind": "bool", "run": "MASTER_MS2_ONLY", "group": "param",
        "zh": "主表只保留有 MS2 证据的注释", "en": "Master table keeps MS2-backed only",
        "stages": ["annotation"]},
    "stage1_quantify": {
        "kind": "bool", "run": "STAGE1_QUANTIFY", "group": "param",
        "zh": "常规注释后是否定量", "en": "Quantify after common annotation",
        "stages": ["annotation"]},
    "eq3_workers": {
        "kind": "int", "run": "EQ3_WORKERS", "group": "param", "min": 1, "max": 64,
        "zh": "搜库进程数（16GB 内存建议 1）", "en": "Library search workers",
        "stages": ["annotation", "search"]},
    "raw_workers": {
        "kind": "int", "run": "raw_workers", "group": "param", "min": 1, "max": 64,
        "zh": "raw 转换并发数（16GB 建议 6，太大要 OOM）", "en": "raw2pkl workers",
        "stages": ["annotation", "search", "quantify"]},

    # ---------------- 分析层（只活在 agent 会话里）----------------
    "alpha": {
        "kind": "float", "run": None, "group": "analysis", "min": 0, "max": 1,
        "zh": "差异分析 FDR 显著性阈值", "en": "FDR significance threshold",
        "stages": ["qc", "interpret", "hypothesis"]},
    "min_n": {
        "kind": "int", "run": None, "group": "analysis", "min": 2, "max": 100,
        "zh": "每组最少几个有值样本才做检验", "en": "Min samples per group to test",
        "stages": ["qc"]},
    "max_hypotheses": {
        "kind": "int", "run": None, "group": "analysis", "min": 1, "max": 12,
        "zh": "最多提几条科学假设", "en": "Max hypotheses", "stages": ["hypothesis"]},
    "lipidmaps_online": {
        "kind": "bool", "run": None, "group": "analysis",
        "zh": "是否连 LIPID MAPS 在线接口", "en": "Query LIPID MAPS online",
        "stages": ["hypothesis"]},
    "sample_context": {
        "kind": "text", "run": None, "group": "analysis",
        "zh": "研究背景（这是什么样本、想比什么）", "en": "Study context",
        "stages": ["hypothesis"]},
    "language": {
        "kind": "text", "run": None, "group": "analysis",
        "zh": "界面与回复语言（zh/en）", "en": "Language (zh/en)", "stages": []},
}

ALIAS = {   # 用户/模型常用的别名 -> 正式键
    "data_root": "project_root", "root": "project_root", "数据目录": "project_root",
    "raw_dir": "raw_root", "raw_root_dir": "raw_root", "原始文件目录": "raw_root",
    "原始数据目录": "raw_root", "raw根目录": "raw_root",
    "pb_dir": "pb_raw_dir", "common_dir": "common_raw_dir",
    "library": "library_dir", "谱库": "library_dir",
    "out_dir": "pb_out_dir", "output_dir": "pb_out_dir",
    "ppm_err1_input": "ppm_err1", "da_err1_input": "da_err1",
    "ppm_err2_input": "ppm_err2", "da_err2_input": "da_err2",
    "MS2_filter": "ms2_filter", "PB_LOC_RT_TOL": "rt_tol",
    "PB_LOC_L12_MIN_COV": "l12_min_coverage",
    "PB_QUANT_RT_CLUSTER_TOL": "rt_cluster_tol",
    "PB_LOC_ISOMER_TOP_N": "isomer_top_n",
    "fdr": "alpha", "q_threshold": "alpha", "context": "sample_context",
}


def canonical(key):
    k = str(key or "").strip()
    if k in SETTINGS:
        return k
    if k in ALIAS:
        return ALIAS[k]
    low = k.lower()
    if low in SETTINGS:
        return low
    return ALIAS.get(low)


def coerce(key, value):
    """把用户/模型给的值转成该项该有的类型并校验范围。

    返回 (ok, 值或错误信息)。路径类会真的去看盘上存不存在——参数写错了当场说，
    比等跑到一半报 FileNotFoundError 强。
    """
    spec = SETTINGS.get(key)
    if not spec:
        return False, f"未知配置项 {key}"
    kind = spec["kind"]
    try:
        if kind in ("dir", "dir_out", "file"):
            p = str(value).strip().strip('"').strip("'")
            if not p:
                return False, "路径为空"
            p = os.path.normpath(os.path.expanduser(p))
            if kind == "dir" and not os.path.isdir(p):
                return False, f"目录不存在：{p}"
            if kind == "file" and not os.path.isfile(p):
                return False, f"文件不存在：{p}"
            return True, p
        if kind == "float":
            v = float(value)
        elif kind == "int":
            v = int(float(value))
        elif kind == "bool":
            s = str(value).strip().lower()
            if s in ("true", "1", "yes", "on", "是", "开", "要"):
                v = True
            elif s in ("false", "0", "no", "off", "否", "关", "不要"):
                v = False
            else:
                return False, f"{value} 不是是/否"
            return True, v
        else:
            return True, str(value)
    except (TypeError, ValueError):
        return False, f"{value} 不是合法的 {kind}"
    lo, hi = spec.get("min"), spec.get("max")
    if lo is not None and v < lo:
        return False, f"{v} 小于下限 {lo}"
    if hi is not None and v > hi:
        return False, f"{v} 超过上限 {hi}"
    return True, v


def validate_patch(patch):
    """逐项校验，返回 (合法项, 错误列表)。非法项不会污染合法项。"""
    good, errors = {}, []
    for raw_key, value in (patch or {}).items():
        key = canonical(raw_key)
        if not key:
            errors.append(f"不认识的配置项「{raw_key}」")
            continue
        ok, v = coerce(key, value)
        if ok:
            good[key] = v
        else:
            errors.append(f"{key}：{v}")
    return good, errors


def apply_to_run(config, run_module):
    """把配置里属于 run.py 的项写回去。返回真正改动的项。"""
    kw = {}
    for key, spec in SETTINGS.items():
        if spec["run"] and key in config and config[key] is not None:
            kw[spec["run"]] = config[key]
    if not kw:
        return {}
    changed = run_module.configure(**kw)
    sync_from_run(config, run_module)
    return changed


def sync_from_run(config, run_module):
    """把 run.py 的现值读回配置。

    必须有这一步：改 raw_root / project_root 会**连带**重算 raw 与 out 目录，只写不读的话
    配置（以及界面上显示的路径）还停在改之前的值，下次全量写回时又把它们推平回去。
    """
    for key, spec in SETTINGS.items():
        if not spec["run"]:
            continue
        v = getattr(run_module, spec["run"], None)
        if v is not None:
            config[key] = str(v) if spec["kind"] in ("dir", "dir_out", "file") else v
    return config


def stale_stages(keys):
    """改了这些配置项后，哪些阶段的产物过期了。"""
    out = []
    for k in keys:
        for st in SETTINGS.get(k, {}).get("stages", []):
            if st not in out:
                out.append(st)
    return out


def describe(config, lang="zh"):
    """当前配置渲染成人话，给"看看现在的参数"用。"""
    lines = []
    for group, title in (("path", "路径" if lang == "zh" else "Paths"),
                         ("param", "搜库/注释参数" if lang == "zh" else "Search parameters"),
                         ("analysis", "分析参数" if lang == "zh" else "Analysis")):
        lines.append(f"**{title}**")
        for key, spec in SETTINGS.items():
            if spec["group"] != group or key not in config:
                continue
            v = config[key]
            if v is None or v == "":
                continue
            lines.append(f"  · {key} = {v}   （{spec[lang if lang in spec else 'zh']}）")
        lines.append("")
    return "\n".join(lines)


def settings_catalog(lang="zh"):
    """给大模型的可配置项清单（键 + 类型 + 含义 + 范围）。"""
    rows = []
    for key, spec in SETTINGS.items():
        r = f'"{key}" ({spec["kind"]})'
        lo, hi = spec.get("min"), spec.get("max")
        if lo is not None or hi is not None:
            r += f" range[{lo},{hi}]"
        r += f": {spec.get(lang, spec['zh'])}"
        rows.append(r)
    return "\n".join(rows)
