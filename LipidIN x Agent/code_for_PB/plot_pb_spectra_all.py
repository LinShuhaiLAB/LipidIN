# -*- coding: utf-8 -*-
"""
为 A1.xlsx 与 B119.xlsx 里的“所有”脂质绘制一级(EIC)/二级谱图 HTML。

复用 plot_pb_spectra.py 的读取器、EIC/MS2 提取与绘图函数，遍历两个注释表的
每一行(每个 lipid_omega_name)，各生成一个 HTML(PB1+PB2 合并，2x2 子图)。

输出 / Output
------------
pos/spectra_html_all/A1/*.html      + index.html   (A1.xlsx 全部)
pos/spectra_html_all/B119/*.html    + index.html   (B119.xlsx 全部)
pos/spectra_html_all/index.html                    (总索引，带搜索框)

所有 HTML 面向文本均为纯 ASCII/英文，避免中文乱码。
"""

import os
import math
import pandas as pd

import plot_pb_spectra as P   # 复用同目录脚本的函数与配置


GROUPS = ["A1", "B119"]
OUT_ROOT = os.path.join(P.POS_DIR, "spectra_html_all")


def build_indices():
    """为两个批次的 PB1/PB2 raw 建立扫描索引(每个 raw 只遍历一次)。"""
    readers, indices = {}, {}
    print("building raw scan index ...")
    for g in GROUPS:
        for tag in P.PB_TAGS:
            raw_path = os.path.join(P.POS_DIR, "%s_%s.raw" % (tag, g))
            if not os.path.isfile(raw_path):
                print("  [warn] missing raw file: %s" % raw_path)
                continue
            reader = P.RawReader(raw_path).__enter__()
            readers[(g, tag)] = reader
            indices[(g, tag)] = P.RawIndex(reader)
            print("  %s_%s: MS1=%d  MS2=%d"
                  % (tag, g, len(indices[(g, tag)].ms1),
                     len(indices[(g, tag)].ms2)))
    return readers, indices


def process_row(row, group, readers, indices):
    """为一行注释生成 per_pb(EIC + MS2)。"""
    targets = P.extract_targets(row)
    per_pb = {}
    for tag in P.PB_TAGS:
        key = (group, tag)
        if key not in indices:
            continue
        reader, index = readers[key], indices[key]
        entry = {"targets": targets}
        mz_target = targets["prec_calc"]
        rt_arr, int_arr = index.eic(mz_target, ppm=P.EIC_PPM)
        int_smooth = P.gaussian_smooth_1d(int_arr, P.EIC_SMOOTH_SIGMA_SCANS)
        entry["eic"] = (rt_arr, int_smooth, mz_target)
        m2 = index.best_ms2(targets["prec_calc"], targets["rt"])
        if m2:
            scan, rt, prec = m2
            mz, inten = reader.get_peaks(scan)
            entry["ms2"] = (mz, inten, scan, rt, prec)
        per_pb[tag] = entry
    return per_pb


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    readers, indices = build_indices()

    master = {}   # group -> list[(name, fname)]
    try:
        for g in GROUPS:
            df = P.load_annotation(g)
            out_dir = os.path.join(OUT_ROOT, g)
            os.makedirs(out_dir, exist_ok=True)
            rows_meta = []
            seen = {}
            n = len(df)
            print("\n=== %s.xlsx : %d rows ===" % (g, n))
            for i, (_, row) in enumerate(df.iterrows(), 1):
                name = str(row[P.NAME_COL])
                if not name or name == "nan":
                    continue
                try:
                    per_pb = process_row(row, g, readers, indices)
                    fig = P.build_figure(name, g, per_pb)
                except Exception as e:
                    print("  [skip] %s: %s: %s" % (name, type(e).__name__, e))
                    continue
                # 文件名唯一化(同名行加序号后缀)
                base = P.sanitize_filename(name)
                if base in seen:
                    seen[base] += 1
                    base = "%s_%d" % (base, seen[base])
                else:
                    seen[base] = 0
                fname = base + ".html"
                fig.write_html(os.path.join(out_dir, fname),
                               include_plotlyjs="directory", full_html=True)
                rows_meta.append((name, fname))
                if i % 25 == 0 or i == n:
                    print("  [%s] %d/%d" % (g, i, n))
            master[g] = rows_meta
            write_group_index(out_dir, g, rows_meta)
            print("  -> %d html written to %s" % (len(rows_meta), out_dir))
    finally:
        for reader in readers.values():
            reader.__exit__(None, None, None)

    write_master_index(master)
    total = sum(len(v) for v in master.values())
    print("\nDone. %d html total. Root: %s" % (total, OUT_ROOT))


# ============================ 索引页 / index pages ============================

_FILTER_JS = (
    "<script>function flt(){var q=document.getElementById('q')"
    ".value.toLowerCase();var n=0;document.querySelectorAll('li.it')"
    ".forEach(function(li){var m=li.textContent.toLowerCase().indexOf(q)>=0;"
    "li.style.display=m?'':'none';if(m)n++;});"
    "document.getElementById('cnt').textContent=n;}</script>"
)
_STYLE = (
    "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:32px;}"
    "h1{font-size:20px;} h2{font-size:16px;margin-top:22px;}"
    "ul{list-style:none;padding:0;columns:3;-webkit-columns:3;}"
    "li{margin:3px 0;break-inside:avoid;} a{color:#2f6fed;text-decoration:none;}"
    "a:hover{text-decoration:underline;} .g{color:#888;font-size:12px;}"
    "#q{padding:6px 10px;width:320px;font-size:14px;margin:8px 0;}"
    "</style>"
)


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def write_group_index(out_dir, group, rows_meta):
    lines = ["<!doctype html><meta charset='utf-8'>",
             "<title>PB spectra - %s</title>" % group, _STYLE, _FILTER_JS,
             "<h1>PB MS1(EIC) / MS2 spectra - %s.xlsx (PB1 &amp; PB2)</h1>" % group,
             "<input id='q' onkeyup='flt()' placeholder='filter by name...'>",
             " <span>matches: <b id='cnt'>%d</b></span>" % len(rows_meta),
             "<ul>"]
    for name, fname in rows_meta:
        lines.append("<li class='it'><a href='%s'>%s</a></li>"
                     % (fname, _esc(name)))
    lines.append("</ul>")
    with open(os.path.join(out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def write_master_index(master):
    total = sum(len(v) for v in master.values())
    lines = ["<!doctype html><meta charset='utf-8'>",
             "<title>PB spectra - all</title>", _STYLE, _FILTER_JS,
             "<h1>PB MS1(EIC) / MS2 spectra - all lipids (PB1 &amp; PB2)</h1>",
             "<p>Total: %d lipids across %d tables. "
             "Per-table index: %s</p>"
             % (total, len(master),
                " ".join("<a href='%s/index.html'>%s</a>" % (g, g) for g in master)),
             "<input id='q' onkeyup='flt()' placeholder='filter by name...'>",
             " <span>matches: <b id='cnt'>%d</b></span>" % total]
    for g, rows_meta in master.items():
        lines.append("<h2>%s.xlsx (%d)</h2>" % (g, len(rows_meta)))
        lines.append("<ul>")
        for name, fname in rows_meta:
            lines.append("<li class='it'><a href='%s/%s'>%s</a> "
                         "<span class='g'>[%s]</span></li>"
                         % (g, fname, _esc(name), g))
        lines.append("</ul>")
    with open(os.path.join(OUT_ROOT, "index.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
