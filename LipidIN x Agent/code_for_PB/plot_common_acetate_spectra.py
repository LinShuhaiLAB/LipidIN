# -*- coding: utf-8 -*-
"""
在 common（未衍生化、负模式）raw 数据里，绘制 A1.xlsx 中所有脂质的
[M+CH3COO]- 乙酸加合形式的 EIC + 对应 MS2 谱图，输出 HTML。

背景 / 逻辑
-----------
A1.xlsx（sample/common/A1.xlsx）里的物质都是 **PB 衍生化后、正模式 [M+H]+** 的记录：
  - mw.            : 未衍生化脂质的中性单同位素质量
  - PB_mw.         : 衍生化后中性质量 (= mw. + C3H6O)
  - calculated m/z : 衍生化后的 [M+H]+（A1 的 adduct 列恒为 +H+）
  - rt(min)        : 保留时间

而 common 下的 raw（MSMS1.raw / MSMS2.raw）是 **未衍生化、负模式** 采集，脂质以
乙酸加合 [M+CH3COO]- 存在（与 common 库 common_PX_CH3COO.pkl 一致）。因此对每个
A1 物质，用其**未衍生化中性质量 mw.** 换算乙酸加合母离子：

    m/z([M+CH3COO]-) = mw. + 59.01385      # 加醋酸阴离子 CH3COO-

然后在 common raw 中：
  * MS1 -> 该乙酸加合 m/z 的 EIC（沿 RT 高斯平滑），并用 A1 的 rt 标注；
  * MS2 -> 母离子匹配该乙酸加合 m/z、保留时间最接近 A1 rt 的二级扫描全谱。

注意：A1 的 calculated m/z1..12 是 **PB 正模式诊断碎片**，与负模式 common MS2 无关，
因此这里不在 MS2 上高亮它们，只标出（乙酸加合）母离子。

同一脂质的多个 C=C 位置异构体共用同一乙酸加合母离子与 RT，故按 (mw., rt) 去重，
每个脂质一张 HTML（标题列出该脂质的所有 C=C 位置）。

输出 / Output
------------
sample/common/acetate_spectra_html/<脂质名>.html + index.html
"""

import os
import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# 复用 PB 绘图脚本里的 raw 读取 / EIC / 棒状谱 / 工具组件（无副作用，import 安全）
from plot_pb_spectra import (
    RawReader, RawIndex,
    gaussian_smooth_1d, add_eic, add_stick_spectrum,
    sanitize_filename, match_peak,
    EIC_PPM, EIC_SMOOTH_SIGMA_SCANS, COLOR_PEAK,
)


# ============================ 配置 / Configuration ============================

COMMON_DIR = r"D:\bio_inf\LipidIN_github\PB_demo\sample\common"
A1_PATH = os.path.join(COMMON_DIR, "A1.xlsx")
OUT_DIR = os.path.join(COMMON_DIR, "acetate_spectra_html")

# common 的两个负模式采集（作为两行子图，类比 PB1/PB2）
RAW_TAGS = [
    ("MSMS1", os.path.join(COMMON_DIR, "MSMS1.raw")),
    ("MSMS2", os.path.join(COMMON_DIR, "MSMS2.raw")),
]

# A1.xlsx 列名
CLASS_COL = "main_class"
ABBR_COL = "lipid_abbr_name"       # 如 PC(32:1)
CC_COL = "C=C location"            # 如 (n-7)
NEUTRAL_COL = "mw."                # 未衍生化脂质中性质量
RT_COL = "rt(min)"

# 乙酸加合：m/z([M+CH3COO]-) = M + 质量(CH3COO-)
#   CH3COO- = C2H3O2 + e- = 59.0133043 + 0.0005486 = 59.0138529
ACETATE_MASS = 59.0138529

# MS2 母离子匹配容差（负模式，乙酸加合）
PRECURSOR_MATCH_DA = 0.7


# ============================ MS2 选取（负模式） ============================

def best_ms2_neg(index, precursor, rt, tol_da=PRECURSOR_MATCH_DA):
    """母离子匹配 + 保留时间最近的二级扫描（不依赖 PB 脚本的全局容差）。"""
    if not index.ms2:
        return None
    for scale in (1.0, 2.0, 3.0):
        cands = [(s, r, p) for (s, r, p) in index.ms2
                 if math.isfinite(p) and abs(p - precursor) <= tol_da * scale]
        if cands:
            return min(cands, key=lambda t: abs(t[1] - rt))
    return None


# ============================ 绘图 / Plotting ============================

def build_figure(title_name, cc_list, acetate_mz, rt, per_tag):
    """per_tag: tag -> {'eic': (rt_arr,int_arr,mz), 'ms2': (mz,inten,scan,rt,prec)}。"""
    rows = len(RAW_TAGS)
    subplot_titles = []
    for tag, _ in RAW_TAGS:
        d = per_tag.get(tag, {})
        eic = d.get("eic")
        ms2 = d.get("ms2")
        s1 = "%s  EIC" % tag
        if eic:
            s1 += "  (m/z %.4f +/- %g ppm)" % (eic[2], EIC_PPM)
        s2 = "%s  MS2" % tag
        if ms2:
            s2 += "  (scan %d, RT %.2f, prec %.3f)" % (ms2[2], ms2[3], ms2[4])
        else:
            s2 += "  (no matching MS2)"
        subplot_titles += [s1, s2]

    fig = make_subplots(rows=rows, cols=2, subplot_titles=subplot_titles,
                        horizontal_spacing=0.08, vertical_spacing=0.13)

    for r, (tag, _) in enumerate(RAW_TAGS, start=1):
        d = per_tag.get(tag, {})
        eic = d.get("eic")
        ms2 = d.get("ms2")

        if eic:
            rt_arr, int_arr, mz_target = eic
            add_eic(fig, rt_arr, int_arr, r, 1, mz_target, mark_rt=rt)
        else:
            fig.add_annotation(text="no EIC", showarrow=False, x=0.5, y=0.5,
                               xref="x domain", yref="y domain", row=r, col=1,
                               font=dict(color="#999"))

        if ms2:
            mz, inten, _, _, prec = ms2
            # 负模式 common MS2：只标乙酸加合母离子，不高亮 PB 正模式诊断碎片
            add_stick_spectrum(fig, mz, inten, r, 2, COLOR_PEAK, precursor=acetate_mz)
        else:
            fig.add_annotation(text="no matching MS2", showarrow=False, x=0.5, y=0.5,
                               xref="x domain", yref="y domain", row=r, col=2,
                               font=dict(color="#999"))

    for r in range(1, rows + 1):
        fig.update_yaxes(title_text="Intensity", row=r, col=1)
        fig.update_yaxes(title_text="Intensity", row=r, col=2)
        fig.update_xaxes(title_text="Retention time (min)", row=r, col=1)
        fig.update_xaxes(title_text="m/z", row=r, col=2)

    cc_txt = ", ".join(cc_list) if cc_list else ""
    title = ("%s  |  [M+CH3COO]- m/z %.4f  |  RT %.2f min  |  C=C: %s"
             % (title_name, acetate_mz, rt, cc_txt))
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=14)),
        template="plotly_white", height=760, width=1200, bargap=0,
        margin=dict(t=90, l=70, r=40, b=60), hovermode="closest",
    )
    return fig


# ============================ 主流程 / Main ============================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_excel(A1_PATH)
    df = df.dropna(subset=[NEUTRAL_COL, RT_COL]).copy()

    # 按 (未衍生化中性质量, rt) 去重 -> 每个唯一脂质一张图；聚合其 C=C 位置
    df["_key"] = list(zip(df[NEUTRAL_COL].round(4), df[RT_COL].round(3)))
    groups = []
    for key, g in df.groupby("_key", sort=False):
        first = g.iloc[0]
        name = str(first[ABBR_COL])
        neutral = float(first[NEUTRAL_COL])
        rt = float(first[RT_COL])
        cc_list = [str(v) for v in g[CC_COL].dropna().unique()]
        groups.append((name, neutral, rt, cc_list))
    print("A1 唯一脂质数 (按 mw+rt 去重): %d" % len(groups))

    # 打开 raw 并建索引
    print("building common raw scan index ...")
    readers, indices = {}, {}
    for tag, path in RAW_TAGS:
        if not os.path.isfile(path):
            print("  [warn] missing raw: %s" % path)
            continue
        reader = RawReader(path).__enter__()
        readers[tag] = reader
        indices[tag] = RawIndex(reader)
        print("  %s: MS1=%d  MS2=%d" % (tag, len(indices[tag].ms1), len(indices[tag].ms2)))

    index_rows = []
    try:
        for name, neutral, rt, cc_list in groups:
            acetate_mz = neutral + ACETATE_MASS

            per_tag = {}
            for tag, _ in RAW_TAGS:
                if tag not in indices:
                    continue
                index = indices[tag]
                reader = readers[tag]
                entry = {}

                rt_arr, int_arr = index.eic(acetate_mz, ppm=EIC_PPM)
                int_smooth = gaussian_smooth_1d(int_arr, EIC_SMOOTH_SIGMA_SCANS)
                entry["eic"] = (rt_arr, int_smooth, acetate_mz)

                m2 = best_ms2_neg(index, acetate_mz, rt)
                if m2:
                    scan, ms2_rt, prec = m2
                    mz, inten = reader.get_peaks(scan)
                    entry["ms2"] = (mz, inten, scan, ms2_rt, prec)

                per_tag[tag] = entry

            fig = build_figure(name, cc_list, acetate_mz, rt, per_tag)
            fname = sanitize_filename(name) + ".html"
            fig.write_html(os.path.join(OUT_DIR, fname),
                           include_plotlyjs="directory", full_html=True)
            index_rows.append((name, acetate_mz, rt, cc_list, fname))
            print("[ok] %-22s acetate m/z %.4f  RT %.2f -> %s" % (name, acetate_mz, rt, fname))
    finally:
        for reader in readers.values():
            reader.__exit__(None, None, None)

    write_index(index_rows)
    print("\nDone. Output dir: %s" % OUT_DIR)


def write_index(index_rows):
    lines = [
        "<!doctype html><meta charset='utf-8'>",
        "<title>common [M+CH3COO]- spectra index</title>",
        "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:40px;}"
        "h1{font-size:20px;} li{margin:5px 0;} a{color:#2f6fed;text-decoration:none;}"
        "a:hover{text-decoration:underline;} .g{color:#888;font-size:12px;}</style>",
        "<h1>common (negative) [M+CH3COO]- EIC / MS2 &mdash; from A1.xlsx lipids</h1>",
        "<ol>",
    ]
    for name, mz, rt, cc_list, fname in index_rows:
        cc = ", ".join(cc_list)
        lines.append("<li><a href='%s'>%s</a> <span class='g'>[M+CH3COO]- %.4f, RT %.2f, C=C %s</span></li>"
                     % (fname, name, mz, rt, cc))
    lines.append("</ol>")
    with open(os.path.join(OUT_DIR, "index.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass
    main()
