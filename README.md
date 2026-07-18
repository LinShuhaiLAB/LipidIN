# LipidIN

Lipid annotation, reverse lipidomics, and carbon–carbon double-bond localization — from raw MS data to a biological report.

*Nature Communications* 2025, 16:4566 · Xu, H. *et al.* — [LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics](https://doi.org/10.1038/s41467-025-59683-5)

> 中文说明见文末 [**中文简介**](#中文简介)。

> ### Maintenance status
> **The Python implementation — `LipidIN × Agent` — is the actively developed version and will keep being updated.**
> The earlier **R version is kept for reproducibility only and will no longer be updated.** New users should use the Python track below.

---

## Figures (as reported in the paper)

| Metric | Value |
|---|---|
| Theoretical library size | 168.6 million fragmentation entries (all chain compositions + C=C double-bond positions) |
| Query throughput | ~70 billion spectral comparisons in < 1 s |
| Annotation false discovery rate | 5.7 % (estimated), across 8,923 lipids |
| WMYn recall change | +20 % (estimated) |
| Agent general library | 24,631,702 spectra; ~66 (pos) / ~95 (neg) subclasses; 91 / 108 subclass×adduct combinations |
| Agent PB library | 2,289,913 spectra; 11 glycerophospholipid / lysophospholipid subclasses |
| Minimum RAM | 16 GB |
| Input formats | Thermo `.raw`, SCIEX `.wiff`, Bruker/Agilent `.d`, `.mzML` |

## What it does (Python · LipidIN × Agent)

The Python version runs the workflow as a web app driven from a chat box:

- **Annotation** — MS2 library search against the theoretical library (chain compositions + C=C double-bond positions).
- **End-to-end analysis** — annotation → quantification → QC (CV / PCA / PLS-DA) → differential analysis → pathway enrichment → LLM biological conclusions → 4-axis verification → bilingual report.
- **Double-bond localization** — Paternò–Büchi (PB) derivatization localizes C=C positions in fatty-acyl chains to omega notation (`n-9`, `n-6`, …), with position-level quantification and differential analysis.
- **Parameters editable mid-run** — change any parameter and re-run only the affected step; structured commands run offline, natural-language requests use an optional LLM.

## Repository layout

```
LipidIN-main/
├─ LipidIN x Agent/                 # Python · actively developed
│   ├─ code/                        #   general untargeted lipidomics   (agent_server.py, port 7860)
│   │   └─ README_EN.pdf / README_ZH.pdf
│   └─ code_for_PB/                 #   PB carbon–carbon double-bond localization (pb_agent_server.py, port 8790)
│       └─ README_EN.pdf / README_ZH.pdf
├─ LipidIN for R (20250601)/        # R framework — reproducibility only, no longer updated
├─ README_LipidIN_for_R.md          # R-track module usage (archived)
└─ README.md                        # this overview
```

## Quick start (Python)

1. Install Python **3.11 (Windows x64)** and the packages listed in the sub-manual (`flask pandas numpy scipy scikit-learn matplotlib polars fisher_py …`).
2. Download the spectral libraries from **[Zenodo → https://zenodo.org/records/21421866](https://zenodo.org/records/21421866)** (see below).
3. Launch the web agent and open it in a browser:
   - General lipidomics: `python agent_server.py` → http://127.0.0.1:7860
   - PB double-bond localization: `python pb_agent_server.py --no-open --port 8790` → http://127.0.0.1:8790
4. Point it at your data folder, type `run all` / `全流程`, and confirm.

Each folder also has a one-shot `run.py` (edit the path lines at the top and run).

**Manuals (self-contained PDF, with a worked demo):**
- General — [English (PDF)](LipidIN%20x%20Agent/code/README_EN.pdf) · [中文 (PDF)](LipidIN%20x%20Agent/code/README_ZH.pdf)
- PB double-bond — [English (PDF)](LipidIN%20x%20Agent/code_for_PB/README_EN.pdf) · [中文 (PDF)](LipidIN%20x%20Agent/code_for_PB/README_ZH.pdf)

## Spectral libraries — where & what

The libraries are too large to ship in this repository and are hosted on Zenodo.

### For LipidIN × Agent (Python) → https://zenodo.org/records/21421866

Three delivery forms of the same theoretical library:

| Library | Format | Spectra | Coverage | Use |
|---|---|---|---|---|
| General (exclude PB) | Python `.pkl` (`pos_common.pkl`, `neg_common.pkl`) | 24,631,702 | ~66 pos / ~95 neg subclasses; 91 / 108 subclass×adduct combos | conventional LC-MS/MS lipidomics |
| PB (only PB) | Python `.pkl` (`PB_PX_H.pkl`, `common_PX_CH3COO.pkl` + MS1 index) | 2,289,913 | 11 glycerophospholipid / lysophospholipid subclasses | PB double-bond localization |
| R-native | `.rda` (`pos_ALL.rda`, `neg_ALL.rda`) | (= general) | same as general | R pipeline |

Subclasses covered by the general library: glycerophospholipids (PC/PE/PS/PG/PI/PA + ether & lyso forms), glycerolipids (TG/DG/MG, OxTG), sphingolipids (Cer/SM/HexCer/gangliosides), sterol esters & bile-acid derivatives, fatty acids & oxidized lipids, and acylated glycolipids (MGDG/DGDG/SQDG). Requires ≥ 16 GB RAM.

(R libraries: the original hierarchical / MS-DIAL libraries are at https://zenodo.org/records/13350719; additional oxidized-lipid spectra at https://zenodo.org/records/13735639.)

## Requirements (summary)

- **OS:** Windows. **RAM:** ≥ 16 GB for the `.pkl` libraries.
- **Python track:** Python 3.11 (x64); Thermo `.raw` via `fisher_py` (.NET); `.wiff` / `.d` via ProteoWizard **msconvert**.

## LipidIN for R (legacy)

The R framework (EQ / LCI / WMYn modules, a GUI, and demos) is in [`LipidIN for R (20250601)/`](LipidIN%20for%20R%20%2820250601%29/) with archived instructions in [README_LipidIN_for_R.md](README_LipidIN_for_R.md). It is kept for reproducing published results and will not be updated; use the Python track for new work.

## Citation & license

Xu, H. *et al.* LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics. *Nat. Commun.* 16, 4566 (2025). https://doi.org/10.1038/s41467-025-59683-5

See [LICENSE](LICENSE).

---

## 中文简介

LipidIN：脂质注释、逆向脂质组学与碳碳双键定位——从原始质谱数据到生物学报告。

> ### 维护状态
> **Python 版实现——`LipidIN × Agent`——是当前持续开发的版本，后续会继续更新。**
> 早期的 **R 版仅作复现之用，后续不再更新。** 新用户请使用下面的 Python 版。

### 指标（引自论文）

| 指标 | 数值 |
|---|---|
| 理论库规模 | 1.686 亿条碎裂谱（覆盖所有链组成 + C=C 双键位置） |
| 检索吞吐 | 1 秒内约 700 亿次谱比对 |
| 注释假发现率 FDR | 5.7 %（估计），覆盖 8,923 种脂质 |
| WMYn 召回变化 | +20 %（估计） |
| Agent 常规库 | 24,631,702 条；约 66（正）/95（负）个亚类；91/108 种「亚类×加合」组合 |
| Agent PB 库 | 2,289,913 条；11 种甘油磷脂/溶血磷脂亚类 |
| 最低内存 | 16 GB |
| 输入格式 | Thermo `.raw`、SCIEX `.wiff`、Bruker/Agilent `.d`、`.mzML` |

### 它能做什么（Python · LipidIN × Agent）

Python 版把流程做成对话框驱动的网页应用：

- **注释**：对理论库（链组成 + C=C 双键位置）做 MS2 搜库。
- **一站式分析**：注释 → 定量 → 质控（CV/PCA/PLS-DA）→ 差异分析 → 通路富集 → 大模型生物学结论 → 四维验证 → 中英双语报告。
- **双键定位**：通过 Paternò–Büchi（PB）衍生化，把脂肪酰链上的 C=C 定位到 omega 位置（`n-9`、`n-6`……），并做位置级定量与差异分析。
- **运行中可改参数**：改任意参数后只重跑受影响的步骤；结构化命令离线可用，自然语言请求可接入可选的大模型。

### 如何快速上手（Python）

1. 安装 Python **3.11（Windows x64）** 及子手册列出的依赖包。
2. 从 **[Zenodo → https://zenodo.org/records/21421866](https://zenodo.org/records/21421866)** 下载谱图库（见下）。
3. 启动网页助手并在浏览器打开：
   - 常规脂质组学：`python agent_server.py` → http://127.0.0.1:7860
   - PB 双键定位：`python pb_agent_server.py --no-open --port 8790` → http://127.0.0.1:8790
4. 把它指向你的数据目录，输入 `全流程` 并确认即可。（各目录也有一键脚本 `run.py`。）

手册（自包含 PDF，含 demo 演示）：
- 常规：[中文 PDF](LipidIN%20x%20Agent/code/README_ZH.pdf) · [英文 PDF](LipidIN%20x%20Agent/code/README_EN.pdf)
- PB 双键：[中文 PDF](LipidIN%20x%20Agent/code_for_PB/README_ZH.pdf) · [英文 PDF](LipidIN%20x%20Agent/code_for_PB/README_EN.pdf)

### 谱图库在哪里、包含什么

谱库体积过大，不随仓库分发，托管在 Zenodo：

- **LipidIN × Agent（Python）用库 → https://zenodo.org/records/21421866**，同一套理论库的三种形态：
  - 常规库（`pos_common.pkl`+`neg_common.pkl`）：24,631,702 条，约 66（正）/95（负）个亚类，91/108 种「亚类×加合」组合；
  - PB 库（`PB_PX_H.pkl`+`common_PX_CH3COO.pkl`+MS1 索引）：2,289,913 条，11 种甘油磷脂/溶血磷脂亚类；
  - R 原生（`pos_ALL.rda`、`neg_ALL.rda`）：与常规库同内容。
  - 常规库覆盖：甘油磷脂（PC/PE/PS/PG/PI/PA + 醚型与溶血型）、甘油酯（TG/DG/MG、OxTG）、鞘脂（Cer/SM/HexCer/神经节苷脂）、固醇酯与胆汁酸衍生物、脂肪酸及氧化脂、酰基化糖脂（MGDG/DGDG/SQDG）。需 ≥ 16 GB 内存。
- （R 库）原始层次库 / MS-DIAL 库见 https://zenodo.org/records/13350719，氧化脂质补充谱见 https://zenodo.org/records/13735639。

### LipidIN for R（旧版）

原始 R 框架（EQ / LCI / WMYn 模块、GUI 与 demo）保留在 [`LipidIN for R (20250601)/`](LipidIN%20for%20R%20%2820250601%29/)，说明见 [README_LipidIN_for_R.md](README_LipidIN_for_R.md)。它用于复现已发表结果，后续不再更新——新项目请使用 Python 版。

### 引用

Xu, H. *et al.* LipidIN: a comprehensive repository for flash platform-independent annotation and reverse lipidomics. *Nat. Commun.* 16, 4566 (2025). https://doi.org/10.1038/s41467-025-59683-5
