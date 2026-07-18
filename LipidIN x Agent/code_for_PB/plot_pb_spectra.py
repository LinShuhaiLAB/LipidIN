# -*- coding: utf-8 -*-
"""
绘制指定脂质分子的一级 (MS1) / 二级 (MS2) 谱图
Plot MS1 / MS2 spectra for user-specified lipids.

数据来源 / Data source
----------------------
文件夹 pos/ 下：
  - 原始文件 raw : PB1_A1.raw / PB2_A1.raw / PB1_B119.raw / PB2_B119.raw
  - 注释文件 xlsx: A1.xlsx / B119.xlsx  (sheet: result_step3_all)

注释表每一行是一个脂质异构体 (lipid_omega_name)，提供：
  - detected/calculated m/z : 母离子 (MS1) 质荷比
  - rt(min)                 : 保留时间
  - detected/calculated m/z1..12 + intensity1..12 : PB 反应产生的二级诊断碎片

思路 / Logic
-----------
同一母离子/保留时间下，不同双键位置的异构体 (如 PC(32:1)(n-7) 与 (n-9))
共用同一张二级谱，只是各自的“诊断碎片”不同。因此：
  * MS1: 取保留时间最接近的一级扫描，绘制全谱，并标出母离子。
  * MS2: 取母离子匹配、保留时间最接近的二级扫描，绘制全谱，并高亮该异构体的诊断碎片。

每个物质的 PB1、PB2 结果画在同一个 HTML 内 (2x2 子图：行=PB1/PB2，列=MS1/MS2)，
鼠标悬停即可看到具体的 m/z 与强度。

输出 / Output
------------
pos/spectra_html/<物质名>.html  以及 index.html 索引页
"""

import os
import re
import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ============================ 配置 / Configuration ============================

POS_DIR = r"D:\bio_inf\LipidIN_github\PB_demo\pos"
OUT_DIR = os.path.join(POS_DIR, "spectra_html")

# 每个 PB 数据集对应两个采集：PB1 / PB2
PB_TAGS = ["PB1", "PB2"]

# 待绘制的物质。group 指定注释表 / raw 文件所属批次：
#   A1   -> 使用 A1.xlsx 与 PB1_A1.raw / PB2_A1.raw
#   B119 -> 使用 B119.xlsx 与 PB1_B119.raw / PB2_B119.raw
# 按用户说明：除最后的 PG(46:5)(n-5) 在 B119，其余全部在 A1。
SUBSTANCES = [
    ("PC(32:1)(n-7)",        "A1"),
    ("PC(32:1)(n-9)",        "A1"),
    ("PE(36:2)(n-9)",        "A1"),
    ("PE(36:1)(n-6)",        "A1"),
    ("PI(40:3)(n-10)",       "A1"),
    ("PI(40:3)(n-13)",       "A1"),
    ("LPC(16:1)(n-10)",      "A1"),
    ("LPC(18:1)(n-6)",       "A1"),
    ("PC(O-34:1)(n-9)",      "A1"),
    ("PC(O-34:2)(n-6,n-9)",  "A1"),
    ("LPE(18:1)(n-3)",       "A1"),
    ("LPE(18:1)(n-9)",       "A1"),
    ("PE(O-34:2)(n-15)",     "A1"),
    ("PE(O-34:2)(n-6)",      "A1"),
    ("PG(46:5)(n-5)",        "B119"),
]

NAME_COL = "lipid_omega_name"        # 物质名列
PRECURSOR_CALC_COL = "calculated m/z"  # 母离子理论质荷比
PRECURSOR_DET_COL = "detected m/z"     # 母离子实测质荷比
RT_COL = "rt(min)"                     # 保留时间
ADDUCT_COL = "adduct"
MAX_FRAG = 12                          # 二级碎片列最多 12 个

# 匹配容差 / matching tolerances
PRECURSOR_MATCH_DA = 0.7   # 选二级扫描时母离子允许偏差 (Da)
FRAG_MATCH_PPM = 30.0      # 碎片/母离子峰在 raw 中匹配容差 (ppm)
FRAG_MATCH_DA = 0.02       # 碎片匹配的绝对容差下限 (Da)

# MS1 -> EIC (提取离子流色谱图) / Extracted Ion Chromatogram
# MS1 面板绘制该物质母离子 m/z 的 EIC：在每个 MS1 扫描中按 ppm 容差提取
# 母离子峰强度，作强度 vs 保留时间曲线，并沿 RT 做高斯平滑。
# 注意：本批 raw 存在约 +23 ppm 系统质量偏差，故容差放宽并设绝对下限，
# 以理论母离子 m/z 直接提取（不重定位，避免弱母离子误锁到干扰峰）。
EIC_PPM = 25.0                 # EIC 提取的质量容差 (ppm)
EIC_TOL_FLOOR_DA = 0.015       # EIC 容差绝对下限 (Da)，覆盖低 m/z 处的偏差
EIC_SMOOTH_SIGMA_SCANS = 2.0   # EIC 沿 RT 的高斯平滑标准差 (以 MS1 扫描数计)
EIC_VIEW_WINDOW_MIN = 1.5      # EIC 视图聚焦到注释 RT ± 该窗口 (min)，避免被别处干扰峰拉伸坐标

# 颜色 / colors
COLOR_PEAK = "#9aa5b1"        # 普通峰 灰
COLOR_PRECURSOR = "#2f6fed"   # 母离子 蓝
COLOR_FRAG = "#e8433a"        # 诊断碎片 红


# ==================== Thermo .raw 读取器 (fisher_py) ====================

class RawReader:
    """精简的 Thermo .raw 读取器，封装 fisher_py。"""

    def __init__(self, path):
        self.path = path
        self._raw = None
        self._scan_cls = None
        self.first_scan = None
        self.last_scan = None

    def __enter__(self):
        from fisher_py.raw_file_reader import RawFileReaderAdapter
        from fisher_py.data import Device
        from fisher_py.data.business import Scan

        raw = RawFileReaderAdapter.file_factory(self.path)
        if raw.is_error:
            raise RuntimeError("无法打开 raw 文件 / cannot open: %s" % self.path)
        raw.select_instrument(Device.MS, 1)
        header = raw.run_header_ex
        self._raw = raw
        self._scan_cls = Scan
        self.first_scan = int(header.first_spectrum)
        self.last_scan = int(header.last_spectrum)
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if self._raw is not None:
            try:
                self._raw.dispose()
            except Exception:
                pass
            self._raw = None
        return False

    def scan_numbers(self):
        return range(self.first_scan, self.last_scan + 1)

    def get_ms_order(self, scan_number):
        event = self._raw.get_scan_event_for_scan_number(scan_number)
        return int(event.ms_order.value), event

    def get_rt(self, scan_number):
        return float(self._raw.retention_time_from_scan_number(scan_number))

    def get_precursor_mz(self, scan_number, event):
        """优先取 trailer 的 monoisotopic m/z，否则取 reaction 的 precursor_mass。"""
        try:
            trailer = self._raw.get_trailer_extra_information(scan_number)
            for label, value in zip(trailer.labels, trailer.values):
                if "monoisotopic" in str(label).lower():
                    try:
                        mono = float(str(value).strip())
                    except Exception:
                        continue
                    if math.isfinite(mono) and mono > 0:
                        return mono
        except Exception:
            pass
        try:
            reaction = event.get_reaction(0)
            prec = float(reaction.precursor_mass)
            if math.isfinite(prec) and prec > 0:
                return prec
        except Exception:
            pass
        return float("nan")

    def get_peaks(self, scan_number):
        """返回 (mz_array, intensity_array)，优先 centroid。"""
        scan = self._scan_cls.from_file(self._raw, scan_number)
        if scan.has_centroid_stream and scan.centroid_scan.length > 0:
            mz = np.asarray(scan.centroid_scan.masses, dtype=np.float64)
            inten = np.asarray(scan.centroid_scan.intensities, dtype=np.float64)
        else:
            seg = scan.segmented_scan
            mz = np.asarray(seg.positions, dtype=np.float64)
            inten = np.asarray(seg.intensities, dtype=np.float64)
        return mz, inten


class RawIndex:
    """一次遍历建立扫描索引，供多个物质复用。"""

    def __init__(self, reader):
        self.reader = reader
        self.ms1 = []       # list[(scan, rt)]
        self.ms2 = []       # list[(scan, rt, precursor)]
        self.ms1_mz = []    # list[np.array]  每个 MS1 扫描按 m/z 升序的质心
        self.ms1_int = []   # list[np.array]
        for sn in reader.scan_numbers():
            try:
                order, event = reader.get_ms_order(sn)
            except Exception:
                continue
            if order == 1:
                self.ms1.append((sn, reader.get_rt(sn)))
                # 预载 MS1 质心峰，供 EIC 提取复用（按 m/z 排序便于二分）
                mz, inten = reader.get_peaks(sn)
                if mz.size and (np.diff(mz) < 0).any():
                    o = np.argsort(mz)
                    mz, inten = mz[o], inten[o]
                self.ms1_mz.append(mz)
                self.ms1_int.append(inten)
            elif order == 2:
                prec = reader.get_precursor_mz(sn, event)
                self.ms2.append((sn, reader.get_rt(sn), prec))
        self._ms1_rt = np.array([r for _, r in self.ms1]) if self.ms1 else np.array([])

    def nearest_ms1(self, rt):
        if not self.ms1:
            return None
        idx = int(np.argmin(np.abs(self._ms1_rt - rt)))
        return self.ms1[idx][0], self.ms1[idx][1]

    def eic(self, mz_target, ppm=EIC_PPM):
        """提取母离子 m/z 的 EIC：返回 (rt_array, intensity_array)。
        每个 MS1 扫描取容差内最强峰的强度，无匹配则为 0。"""
        n = len(self.ms1)
        inten_out = np.zeros(n, dtype=np.float64)
        if n == 0 or not math.isfinite(mz_target):
            return self._ms1_rt, inten_out
        tol = max(mz_target * ppm * 1e-6, EIC_TOL_FLOOR_DA)
        for i in range(n):
            mzs = self.ms1_mz[i]
            if mzs.size == 0:
                continue
            j = int(np.searchsorted(mzs, mz_target))
            best = 0.0
            for k in (j - 1, j):
                if 0 <= k < mzs.size and abs(mzs[k] - mz_target) <= tol:
                    v = self.ms1_int[i][k]
                    if v > best:
                        best = v
            inten_out[i] = best
        return self._ms1_rt, inten_out

    def best_ms2(self, precursor, rt):
        """母离子匹配 + 保留时间最近的二级扫描。"""
        if not self.ms2:
            return None
        cands = [(s, r, p) for (s, r, p) in self.ms2
                 if math.isfinite(p) and abs(p - precursor) <= PRECURSOR_MATCH_DA]
        if not cands:
            # 放宽母离子容差
            cands = [(s, r, p) for (s, r, p) in self.ms2
                     if math.isfinite(p) and abs(p - precursor) <= PRECURSOR_MATCH_DA * 2]
        if not cands:
            return None
        s, r, p = min(cands, key=lambda t: abs(t[1] - rt))
        return s, r, p


# ============================ 工具函数 / Helpers ============================

def sanitize_filename(name):
    """把物质名转为合法 Windows 文件名。"""
    return re.sub(r'[:\\/*?"<>|]', "-", name)


def match_peak(mz_array, target):
    """在峰列表中找与 target 最接近且落在容差内的峰，返回 (mz, intensity_index) 或 None。"""
    if mz_array.size == 0 or not math.isfinite(target):
        return None
    j = int(np.argmin(np.abs(mz_array - target)))
    tol = max(target * FRAG_MATCH_PPM * 1e-6, FRAG_MATCH_DA)
    if abs(mz_array[j] - target) <= tol:
        return j
    return None


def extract_targets(row):
    """从注释行提取母离子、rt、诊断碎片(理论/实测 m/z + 强度)。"""
    prec_calc = float(row[PRECURSOR_CALC_COL])
    prec_det = float(row[PRECURSOR_DET_COL])
    rt = float(row[RT_COL])
    adduct = str(row.get(ADDUCT_COL, ""))
    frags = []
    for i in range(1, MAX_FRAG + 1):
        c = row.get("calculated m/z%d" % i)
        d = row.get("detected m/z%d" % i)
        it = row.get("intensity%d" % i)
        if pd.notna(c):
            frags.append({
                "calc": float(c),
                "det": float(d) if pd.notna(d) else float("nan"),
                "intensity": float(it) if pd.notna(it) else float("nan"),
            })
    return {"prec_calc": prec_calc, "prec_det": prec_det, "rt": rt,
            "adduct": adduct, "frags": frags}


def gaussian_smooth_1d(y, sigma_points):
    """对一维序列做高斯平滑（沿 RT），sigma 以点(扫描)数计。"""
    y = np.asarray(y, dtype=np.float64)
    if y.size == 0 or sigma_points <= 0:
        return y
    radius = max(1, int(round(3.0 * sigma_points)))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma_points) ** 2)
    kernel /= kernel.sum()
    return np.convolve(y, kernel, mode="same")


# ============================ 绘图 / Plotting ============================

def add_eic(fig, rt, inten, row, col, mz_target, mark_rt=None):
    """在指定子图绘制母离子 m/z 的高斯平滑 EIC，并标出鉴定所在 RT。"""
    if rt is None or len(rt) == 0 or float(np.max(inten)) <= 0:
        fig.add_annotation(text="no EIC signal", showarrow=False,
                           x=0.5, y=0.5, xref="x domain", yref="y domain",
                           row=row, col=col, font=dict(color="#999"))
        return

    fig.add_trace(
        go.Scatter(
            x=rt, y=inten, mode="lines",
            line=dict(color=COLOR_PRECURSOR, width=1.2),
            fill="tozeroy", fillcolor="rgba(47,111,237,0.18)",
            hovertemplate="RT: %{x:.3f} min<br>intensity: %{y:,.0f}<extra></extra>",
            showlegend=False,
        ),
        row=row, col=col,
    )

    # 鉴定/取二级谱所在保留时间：红色虚线标注
    if mark_rt is not None and math.isfinite(mark_rt):
        fig.add_vline(
            x=float(mark_rt), line=dict(color=COLOR_FRAG, width=1.2, dash="dot"),
            row=row, col=col,
            annotation_text="RT %.2f" % mark_rt,
            annotation_position="top",
            annotation_font=dict(size=9, color=COLOR_FRAG),
        )
        # 视图聚焦到注释 RT ± 窗口，y 轴按该窗口内峰高自适应，
        # 避免别处同质异构干扰峰把坐标撑大而看不清目标峰
        w = EIC_VIEW_WINDOW_MIN
        rt_np = np.asarray(rt, dtype=np.float64)
        int_np = np.asarray(inten, dtype=np.float64)
        m = (rt_np >= mark_rt - w) & (rt_np <= mark_rt + w)
        if m.any() and float(int_np[m].max()) > 0:
            fig.update_xaxes(range=[mark_rt - w, mark_rt + w], row=row, col=col)
            fig.update_yaxes(range=[0, float(int_np[m].max()) * 1.12],
                             row=row, col=col)


# ============================ 绘图 / Plotting ============================

def add_stick_spectrum(fig, mz, inten, row, col, base_color,
                       highlights=None, precursor=None):
    """
    在指定子图绘制棒状谱 (stick spectrum)。
    highlights: [{"mz": float, "label": str}]  诊断碎片，高亮成红色并标注。
    precursor : float  母离子 m/z，高亮成蓝色。
    """
    if mz.size == 0:
        fig.add_annotation(text="no spectrum", showarrow=False,
                           x=0.5, y=0.5, xref="x domain", yref="y domain",
                           row=row, col=col, font=dict(color="#999"))
        return

    # 每个峰上色：默认灰，母离子蓝，诊断碎片红
    colors = np.full(mz.size, base_color, dtype=object)

    special_idx = {}  # index -> (color, label)
    if precursor is not None:
        j = match_peak(mz, precursor)
        if j is not None:
            special_idx[j] = (COLOR_PRECURSOR, "precursor %.4f" % mz[j])

    if highlights:
        for h in highlights:
            j = match_peak(mz, h["mz"])
            if j is not None:
                special_idx[j] = (COLOR_FRAG, h["label"])

    for j, (c, _lab) in special_idx.items():
        colors[j] = c

    # 竖线（棒状谱本体）：分三组绘制 灰/蓝/红，保证高亮在上层且颜色正确
    _add_lines(fig, mz, inten, colors, row, col)

    # 峰顶 marker，承载 hover
    fig.add_trace(
        go.Scatter(
            x=mz, y=inten, mode="markers",
            marker=dict(size=4, color=list(colors)),
            hovertemplate="m/z: %{x:.4f}<br>intensity: %{y:,.0f}<extra></extra>",
            showlegend=False,
        ),
        row=row, col=col,
    )

    # 诊断碎片 / 母离子 文本标注
    for j, (c, lab) in special_idx.items():
        fig.add_annotation(
            x=mz[j], y=inten[j], text="%.4f" % mz[j], showarrow=True,
            arrowhead=0, arrowwidth=1, arrowcolor=c, ax=0, ay=-18,
            font=dict(size=9, color=c),
            row=row, col=col,
        )


def _add_lines(fig, mz, inten, colors, row, col):
    """按颜色分组画竖线，避免单条 trace 混色。"""
    for target_color in (COLOR_PEAK, COLOR_PRECURSOR, COLOR_FRAG):
        idx = [j for j in range(mz.size) if colors[j] == target_color]
        if not idx:
            continue
        xs, ys = [], []
        for j in idx:
            xs += [mz[j], mz[j], None]
            ys += [0, inten[j], None]
        lw = 2.4 if target_color != COLOR_PEAK else 1.0
        fig.add_trace(
            go.Scatter(x=xs, y=ys, mode="lines",
                       line=dict(color=target_color, width=lw),
                       hoverinfo="skip", showlegend=False),
            row=row, col=col,
        )


def _axref(row, col):
    """subplot (row,col) -> plotly 轴编号 (2 列布局)。"""
    return (row - 1) * 2 + col


def build_figure(name, group, per_pb):
    """
    per_pb: dict tag -> {
        "targets": targets,
        "eic": (rt_array, intensity_array, mz_target) 或 None,   # MS1 -> EIC
        "ms2": (mz, inten, ms2_scan, ms2_rt, ms2_prec) 或 None,
    }
    """
    rows = len(PB_TAGS)
    subplot_titles = []
    for tag in PB_TAGS:
        d = per_pb.get(tag, {})
        eic = d.get("eic")
        ms2 = d.get("ms2")
        s1 = "%s  EIC" % tag
        if eic:
            s1 += "  (m/z %.4f +/- %g ppm)" % (eic[2], EIC_PPM)
        s2 = "%s  MS2" % tag
        if ms2:
            s2 += "  (scan %d, RT %.2f min, prec %.3f)" % (ms2[2], ms2[3], ms2[4])
        subplot_titles += [s1, s2]

    fig = make_subplots(
        rows=rows, cols=2,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.08, vertical_spacing=0.13,
    )

    any_targets = None
    for r, tag in enumerate(PB_TAGS, start=1):
        d = per_pb.get(tag, {})
        targets = d.get("targets")
        any_targets = any_targets or targets
        eic = d.get("eic")
        ms2 = d.get("ms2")

        # ---- MS1 (左列)：母离子 EIC，沿 RT 高斯平滑 ----
        if eic:
            rt_arr, int_arr, mz_target = eic
            mark_rt = targets["rt"] if targets else None
            add_eic(fig, rt_arr, int_arr, r, 1, mz_target, mark_rt=mark_rt)
        else:
            _empty(fig, r, 1)

        # ---- MS2 (右列) ----
        if ms2:
            mz, inten, _, _, _ = ms2
            highlights = []
            if targets:
                for fr in targets["frags"]:
                    highlights.append({"mz": fr["calc"], "label": "%.4f" % fr["calc"]})
            add_stick_spectrum(fig, mz, inten, r, 2, COLOR_PEAK, highlights=highlights)
        else:
            _empty(fig, r, 2)

    # 轴标题：左列 EIC 的 x 轴为保留时间，右列 MS2 的 x 轴为 m/z
    for r in range(1, rows + 1):
        fig.update_yaxes(title_text="Intensity", row=r, col=1)
        fig.update_yaxes(title_text="Intensity", row=r, col=2)
        fig.update_xaxes(title_text="Retention time (min)", row=r, col=1)
        fig.update_xaxes(title_text="m/z", row=r, col=2)

    tt = any_targets or {}
    title = ("%s  |  group %s  |  precursor calc %.4f / det %.4f  |  RT %.2f min  |  adduct %s"
             % (name, group,
                tt.get("prec_calc", float("nan")),
                tt.get("prec_det", float("nan")),
                tt.get("rt", float("nan")),
                tt.get("adduct", "")))
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=15)),
        template="plotly_white",
        height=760, width=1200,
        bargap=0,
        margin=dict(t=90, l=70, r=40, b=60),
        hovermode="closest",
    )
    return fig


def _empty(fig, row, col):
    fig.add_annotation(text="no matching scan", showarrow=False,
                       x=0.5, y=0.5, xref="x domain", yref="y domain",
                       row=row, col=col, font=dict(color="#999"))


# ============================ 主流程 / Main ============================

def load_annotation(group):
    path = os.path.join(POS_DIR, "%s.xlsx" % group)
    return pd.read_excel(path)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 需要用到哪些批次
    groups = sorted({g for _, g in SUBSTANCES})
    annotations = {g: load_annotation(g) for g in groups}

    # 为每个 (批次, PB) 建立 raw 扫描索引，避免重复遍历
    print("building raw scan index ...")
    raw_indices = {}   # (group, tag) -> RawIndex
    raw_readers = {}   # (group, tag) -> RawReader (保持打开)
    for g in groups:
        for tag in PB_TAGS:
            raw_path = os.path.join(POS_DIR, "%s_%s.raw" % (tag, g))
            if not os.path.isfile(raw_path):
                print("  [warn] missing raw file: %s" % raw_path)
                continue
            reader = RawReader(raw_path).__enter__()
            raw_readers[(g, tag)] = reader
            raw_indices[(g, tag)] = RawIndex(reader)
            print("  %s_%s: MS1=%d  MS2=%d"
                  % (tag, g, len(raw_indices[(g, tag)].ms1),
                     len(raw_indices[(g, tag)].ms2)))

    index_rows = []
    try:
        for name, group in SUBSTANCES:
            df = annotations[group]
            sub = df[df[NAME_COL].astype(str) == name]
            if sub.empty:
                print("[skip] %s not found in %s.xlsx" % (name, group))
                continue
            row = sub.iloc[0]
            targets = extract_targets(row)

            per_pb = {}
            for tag in PB_TAGS:
                key = (group, tag)
                if key not in raw_indices:
                    continue
                reader = raw_readers[key]
                index = raw_indices[key]
                entry = {"targets": targets}

                # MS1 -> 母离子 EIC (沿 RT 高斯平滑)
                # 以理论母离子 m/z 直接提取（容差已放宽覆盖系统质量偏差）
                mz_target = targets["prec_calc"]
                rt_arr, int_arr = index.eic(mz_target, ppm=EIC_PPM)
                int_smooth = gaussian_smooth_1d(int_arr, EIC_SMOOTH_SIGMA_SCANS)
                entry["eic"] = (rt_arr, int_smooth, mz_target)

                # MS2
                m2 = index.best_ms2(targets["prec_calc"], targets["rt"])
                if m2:
                    scan, rt, prec = m2
                    mz, inten = reader.get_peaks(scan)
                    entry["ms2"] = (mz, inten, scan, rt, prec)

                per_pb[tag] = entry

            fig = build_figure(name, group, per_pb)
            fname = sanitize_filename(name) + ".html"
            out_path = os.path.join(OUT_DIR, fname)
            fig.write_html(out_path, include_plotlyjs="directory",
                           full_html=True)
            index_rows.append((name, group, fname))
            print("[ok] %s -> %s" % (name, fname))
    finally:
        for reader in raw_readers.values():
            reader.__exit__(None, None, None)

    # 索引页
    write_index(index_rows)
    print("\nDone. Output dir: %s" % OUT_DIR)


def write_index(index_rows):
    lines = [
        "<!doctype html><meta charset='utf-8'>",
        "<title>PB spectra index</title>",
        "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:40px;}"
        "h1{font-size:20px;} li{margin:6px 0;} a{color:#2f6fed;text-decoration:none;}"
        "a:hover{text-decoration:underline;} .g{color:#888;font-size:12px;}</style>",
        "<h1>PB MS1 / MS2 spectra (PB1 &amp; PB2)</h1>",
        "<ol>",
    ]
    for name, group, fname in index_rows:
        lines.append("<li><a href='%s'>%s</a> <span class='g'>[%s]</span></li>"
                     % (fname, name, group))
    lines.append("</ol>")
    with open(os.path.join(OUT_DIR, "index.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
