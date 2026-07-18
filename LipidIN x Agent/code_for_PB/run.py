import os
import sys
import time
from pathlib import Path

# 管线里大量 print 含 emoji/特殊符号；Windows 默认 gbk 控制台无法编码会抛
# UnicodeEncodeError。入口处统一把标准输出改成 utf-8，覆盖所有下游模块的打印。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8")
    except Exception:
        pass

from tm_raw2pkl import *
from EQ2 import *
from LCI2 import *
from master_table import *
from quantify import *
from pb_search import run_pb_search
from pb_localize import localize_all_pb_diagnostic
from pb_quantify import quantify_pb_dbpositions

# ============================================================
# 路径配置（PB_demo）
# ============================================================
CODE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CODE_DIR.parent                       # D:\bio_inf\LipidIN_github\PB_demo
SAMPLE_ROOT = PROJECT_ROOT / "sample"
LIBRARY_DIR = PROJECT_ROOT / "PB library"
RAW_ROOT = SAMPLE_ROOT      # "原始文件都在这"：其下 common/ 与 PB/ 两个子目录分别是两种采集

# 第一阶段：常规注释（sample/common，负模式 CH3COO 库）
COMMON_RAW_DIR = SAMPLE_ROOT / "common"
COMMON_OUT_DIR = COMMON_RAW_DIR / "pkl"
COMMON_LIBRARY = LIBRARY_DIR / "common_PX_CH3COO.pkl"

# 第二阶段：PB 衍生化检索（sample/PB，正模式 [M+H]+ PB 库）
PB_RAW_DIR = SAMPLE_ROOT / "PB"
PB_OUT_DIR = PB_RAW_DIR / "pkl"
PB_LIBRARY = LIBRARY_DIR / "PB_PX_H.pkl"

CPP_FILE = str(CODE_DIR / "eq3_similarity.cpp")

# ============================================================
# 可调参数
# ============================================================
RUN_STAGE1 = True      # 常规注释
RUN_STAGE2 = True      # PB 检索
RUN_STAGE3 = True      # PB 双键位置重建
RUN_STAGE4 = True      # PB 链×双键位置级定量 + 总量归一
STAGE1_QUANTIFY = True  # 常规注释后是否定量
MASTER_MS2_ONLY = True  # 主表只保留有 MS2 证据的注释（丢弃 MS1_only/unannotated 的 MS1 推断）

PB_LOC_RT_TOL = 0.5    # 第三阶段：PB 峰与 master 就近匹配的 RT 容差(min)。
                       # master 只用于细化链划分、不做闸门，故放宽只影响“能否用上
                       # master 的链划分”，不会因对不上而丢峰。PB 衍生化会让 RT 相对
                       # common 漂移，0.3 会漏掉真实对应(如 PE 36:1 漂移 0.34)。
PB_LOC_ISOMER_TOP_N = 20   # 每个 (峰,物质) 最多写出多少个竞争异构体行(按 位置数/置信度 降序)
PB_LOC_L12_MIN_COV = 0.5   # 第三阶段：common 无背书的峰，需 LEVEL1(母离子对)与 LEVEL2(头基对)
                           # 各自碎片覆盖率 ≥ 此值才认物质成立。每个 level 各 2 根碎片，
                           # 故 0.5 = "两根里至少命中一根"。对上 master 的峰有 common 背书，
                           # 不走这道门。注:非胆碱类(PE/PE-O/PI/PG/LPE)库里 L2 与 L1 是同一
                           # 对离子，该判据对它们等价于单一条件。

PB_QUANT_RT_CLUSTER_TOL = 0.3  # 第四阶段：跨样本对齐同一 (物质,异构体) 时的 RT 聚类容差(min)。
                               # 相差超过它视为不同色谱峰、分成两行，不强行合并。

EQ3_WORKERS = 1        # 搜库进程数（16GB 内存对应 1）
MS2_filter = 0.0005    # 过滤低丰度二级碎片
ppm_err1_input = 20    # MS1 ppm 误差（common 与 PB 一致，MS1/MS2 都 20ppm）
da_err1_input = 0.01   # MS1 道尔顿误差
ppm_err2_input = 20    # MS2 ppm 误差
da_err2_input = 0.05   # MS2 道尔顿误差（第一阶段 common 用）
pb_ppm_err1 = 60       # 第二阶段 PB 诊断搜库单独用更松的【母离子】容差(common 保持 20)。
pb_da_err1 = 0.05      # 原因:这批 raw 的 DDA 母离子被量化到 0.1 网格(PB1_A1 的 1983 个 MS2
                       # 只有 82 个唯一母离子值、100% 落在 0.1 网格上)，raw 本身给不出优于
                       # ±0.05 的母离子。tm_raw2pkl 的精修只能在隔离窗内取最强 MS1 峰去猜，
                       # 目标若非窗内最丰度物种就会被共洗脱物吸走——同一化合物的精修母离子
                       # 在扫描间从 -144 飘到 +84 ppm。用 20ppm 卡这个"猜出来的"母离子会把
                       # 真注释误杀:实测携带完整诊断离子对的扫描母离子偏差全在 20-58ppm，
                       # 无一落在 20ppm 内(LPC 17:1 的 n-8 两根离子实测只差 -7.4/-3.6ppm)。
                       # 故母离子只用于粗筛(≈隔离窗量级)，真正的鉴别交给高精度诊断离子对。
pb_da_err2 = 0.02      # 第二阶段 PB 诊断搜库单独用更紧的碎片容差。PB 诊断离子往往近噪声，
                       # 0.05 太松会把偏 ~30ppm 的噪声峰当成诊断离子(如把 PC42:4 误定到
                       # n-6,9,12 盖过真实 n-13)。收到 0.02(~25ppm@800) 保住亚ppm真离子、砍噪声。
raw_workers = 6        # raw2pkl 并发进程数。None→约 cpu*0.85(≈20+)，在 16GB 机器上
                       # 会因 24 个 worker 各持大 MS1 池同时 pickle.dump(ms1_only) 而 OOM
                       # (报错 MemoryError、文件被标 failed)。限到 6 兼顾速度与内存安全。

# 派生路径
common_aligned_table = str(COMMON_OUT_DIR / "rt_ccs_aligned_table.csv")
common_rt_checked_table = str(COMMON_OUT_DIR / "rt_ccs_aligned_table_RT_checked.csv")
common_master_table = str(COMMON_OUT_DIR / "master_table.csv")
common_master_quant = str(COMMON_OUT_DIR / "master_table_quant.csv")


# ============================================================
# 运行时改配置（给 agent 用；直接跑 run.py 的人不受影响）
# ============================================================
# 上面这些名字都是模块级全局，而各 stage 函数是在**被调用时**才读它们的，所以运行时
# 改了立刻生效。唯一的坑是派生路径：common_master_table 这类是 import 时算好的字符串，
# 改了 COMMON_OUT_DIR 而不重算，它们还会指着老目录——这正是 _recompute_derived 存在的理由。
_PATH_KEYS = {"PROJECT_ROOT", "SAMPLE_ROOT", "LIBRARY_DIR", "RAW_ROOT", "COMMON_RAW_DIR",
              "COMMON_OUT_DIR", "COMMON_LIBRARY", "PB_RAW_DIR", "PB_OUT_DIR", "PB_LIBRARY"}


def raw_files(folder):
    """目录**顶层**的 .raw 列表（raw2pkl 也是 recursive=0，这里必须与它口径一致）。"""
    folder = Path(folder)
    if not folder.is_dir():
        return []
    return sorted(p for p in folder.iterdir() if p.is_file() and p.suffix.lower() == ".raw")


def subdirs_with_raw(folder):
    """{子目录名: 路径}，只列真的装着 .raw 的子目录。用于把"指错一层"说清楚。"""
    folder = Path(folder)
    if not folder.is_dir():
        return {}
    return {p.name: p for p in sorted(folder.iterdir()) if p.is_dir() and raw_files(p)}


def require_raw(folder, what):
    """转换前的硬闸门：没有 raw 就当场报错，并指出 raw 到底在哪一层。

    以前这里只是让 raw2pkl 打个"未找到 .raw"的警告然后继续，结果三步之后炸在
    find_ms1_msms_files 上，报的是"没有 pkl"——真因（raw 目录指错一层）被掩盖了。
    """
    files = raw_files(folder)
    if files:
        return files
    subs = subdirs_with_raw(folder)
    hint = ""
    if subs:
        hint = ("。不过它的这些子目录里有 .raw：" + "、".join(subs) +
                f"——你多半是把目录指高了一层。可以直接说「raw 根目录 {folder}」"
                "让我按 common/PB 自动分，或分别指定 common raw 目录 / PB raw 目录")
    raise FileNotFoundError(f"{what} 的 raw 目录里没有 .raw 文件：{folder}{hint}")


def set_raw_root(root):
    """告诉管线"我的原始文件都在这"：<root>/common 是常规采集、<root>/PB 是 PB 衍生化采集。

    两批 raw 必须分在两个子目录里——它们是不同的采集模式（负模式 CH3COO / 正模式 [M+H]+）、
    走不同的库，混在一个目录里没法自动分辨。pkl 输出各自落在 <raw 目录>/pkl。
    """
    global RAW_ROOT, COMMON_RAW_DIR, COMMON_OUT_DIR, PB_RAW_DIR, PB_OUT_DIR
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"raw 根目录不存在 {root}")
    subs = subdirs_with_raw(root)
    common = next((p for n, p in subs.items() if "common" in n.lower() or "常规" in n), None)
    pb = next((p for n, p in subs.items() if n.lower().startswith("pb")), None)
    if common is None or pb is None:
        missing = "common" if common is None else "PB"
        raise FileNotFoundError(
            f"{root} 下没找到装着 .raw 的 {missing} 子目录。约定布局是 "
            f"<raw 根目录>/common（常规采集）+ <raw 根目录>/PB（PB 衍生化采集）；"
            f"当前含 .raw 的子目录：{'、'.join(subs) or '一个都没有'}。"
            f"不按这个布局的话，用 common_raw_dir / pb_raw_dir 分别指定。")
    RAW_ROOT = root
    COMMON_RAW_DIR, COMMON_OUT_DIR = common, common / "pkl"
    PB_RAW_DIR, PB_OUT_DIR = pb, pb / "pkl"
    _recompute_derived()


def _recompute_derived():
    global common_aligned_table, common_rt_checked_table, common_master_table, common_master_quant
    common_aligned_table = str(COMMON_OUT_DIR / "rt_ccs_aligned_table.csv")
    common_rt_checked_table = str(COMMON_OUT_DIR / "rt_ccs_aligned_table_RT_checked.csv")
    common_master_table = str(COMMON_OUT_DIR / "master_table.csv")
    common_master_quant = str(COMMON_OUT_DIR / "master_table_quant.csv")


def set_project_root(root):
    """把整套路径挪到新的项目根目录（按 PB_demo 的标准布局重建）。

    布局约定：<root>/sample/common、<root>/sample/PB、<root>/PB library。
    个别目录不按约定的，再用 configure() 单独覆盖即可。
    """
    global PROJECT_ROOT, SAMPLE_ROOT, LIBRARY_DIR, RAW_ROOT, COMMON_RAW_DIR, COMMON_OUT_DIR
    global COMMON_LIBRARY, PB_RAW_DIR, PB_OUT_DIR, PB_LIBRARY
    PROJECT_ROOT = Path(root)
    SAMPLE_ROOT = PROJECT_ROOT / "sample"
    # 新根目录下没有库就留用原来的：用户挪的往往只是数据目录，库还在老地方，
    # 跟着挪过去只会在下一步报"常规库不存在"。
    LIBRARY_DIR = (PROJECT_ROOT / "PB library") if (PROJECT_ROOT / "PB library").is_dir() \
        else LIBRARY_DIR
    # raw 通常在 <root>/sample 下；但用户常常直接把装着 common/ 与 PB/ 的那个目录当根目录给我，
    # 这时硬套 sample/ 会指到一个不存在的目录上。谁下面真有 common+PB 就用谁。
    RAW_ROOT = SAMPLE_ROOT if (SAMPLE_ROOT / "common").is_dir() else PROJECT_ROOT
    COMMON_RAW_DIR = RAW_ROOT / "common"
    COMMON_OUT_DIR = COMMON_RAW_DIR / "pkl"
    COMMON_LIBRARY = LIBRARY_DIR / "common_PX_CH3COO.pkl"
    PB_RAW_DIR = RAW_ROOT / "PB"
    PB_OUT_DIR = PB_RAW_DIR / "pkl"
    PB_LIBRARY = LIBRARY_DIR / "PB_PX_H.pkl"
    _recompute_derived()


def configure(**kwargs):
    """运行时设置任意模块级配置项，返回真正改动了的 {名字: 新值}。

    未知的名字直接抛 KeyError——宁可报错，也不要让一个拼错的参数名被静默忽略、
    让人以为改生效了。
    """
    g = globals()
    changed = {}
    # 两个"根目录"是宏：一句话重排下面一整组路径。调用方（agent）每次都把整份配置全量写回，
    # 所以只在**真的变了**的时候才展开，否则每次调用都会把用户单独改过的子路径推平。
    root = kwargs.pop("PROJECT_ROOT", None)
    if root and Path(root) != PROJECT_ROOT:
        set_project_root(root)
        for k in ("PROJECT_ROOT", "RAW_ROOT", "COMMON_RAW_DIR", "COMMON_OUT_DIR",
                  "PB_RAW_DIR", "PB_OUT_DIR", "LIBRARY_DIR", "COMMON_LIBRARY", "PB_LIBRARY"):
            changed[k] = str(g[k])
            kwargs.pop(k, None)   # 根目录已经重排过这些了，别让同一次调用捎带的旧值覆盖回去
    raw_root = kwargs.pop("RAW_ROOT", None)
    if raw_root and Path(raw_root) != RAW_ROOT:
        set_raw_root(raw_root)
        for k in ("RAW_ROOT", "COMMON_RAW_DIR", "COMMON_OUT_DIR", "PB_RAW_DIR", "PB_OUT_DIR"):
            changed[k] = str(g[k])
            kwargs.pop(k, None)   # 别让同一次调用里带着的老 raw/out 值把刚算出来的覆盖回去

    before = {k: g[k] for k in _PATH_KEYS}
    for k, v in kwargs.items():
        if k not in g:
            raise KeyError(f"未知配置项 {k}")
        if k in _PATH_KEYS:
            v = Path(v)
        if g[k] != v:
            g[k] = v
            changed[k] = str(v)
    # raw 目录换了、而输出目录还停在**老 raw 目录下的默认位置**，说明用户没特意指定过输出目录，
    # 让它跟着 raw 走。不这样的话 pkl 写到新地方、却仍从老目录去读，报的还是"没有 pkl"。
    for raw_key, out_key in (("COMMON_RAW_DIR", "COMMON_OUT_DIR"), ("PB_RAW_DIR", "PB_OUT_DIR")):
        if g[raw_key] != before[raw_key] and g[out_key] == before[raw_key] / "pkl":
            g[out_key] = g[raw_key] / "pkl"
            changed[out_key] = str(g[out_key])
    _recompute_derived()
    return changed


def current_paths():
    """当前所有路径，给界面展示 / 校验用。"""
    return {k: str(globals()[k]) for k in sorted(_PATH_KEYS)}


def check_inputs():
    if RUN_STAGE1:
        if not COMMON_RAW_DIR.is_dir():
            raise FileNotFoundError(f"common raw 目录不存在 {COMMON_RAW_DIR}")
        if not COMMON_LIBRARY.exists():
            raise FileNotFoundError(f"常规库不存在 {COMMON_LIBRARY}")
        print(f"[信息] 第一阶段 common raw {len(require_raw(COMMON_RAW_DIR, '常规注释'))} 个，"
              f"库 {COMMON_LIBRARY.name}")
    if RUN_STAGE2:
        if not PB_RAW_DIR.is_dir():
            raise FileNotFoundError(f"PB raw 目录不存在 {PB_RAW_DIR}")
        if not PB_LIBRARY.exists():
            raise FileNotFoundError(f"PB 库不存在 {PB_LIBRARY}")
        print(f"[信息] 第二阶段 PB raw {len(require_raw(PB_RAW_DIR, 'PB 衍生化检索'))} 个，"
              f"库 {PB_LIBRARY.name}")
    if not os.path.exists(CPP_FILE):
        raise FileNotFoundError(f"C++ 相似度源文件不存在 {CPP_FILE}")


def stage1_common_annotation():
    """第一阶段：sample/common 常规注释（流程与既有逻辑一致，仅路径/库改为 PB_demo）。"""
    print("\n" + "=" * 60)
    print("第一阶段 | 常规注释 sample/common")
    print("=" * 60)

    n_raw = len(require_raw(COMMON_RAW_DIR, "常规注释"))
    print(f"[信息] raw {n_raw} 个：{COMMON_RAW_DIR} -> pkl：{COMMON_OUT_DIR}")

    # ThermoRawToPklAgent 只在构造时保存参数，需显式调用 run() 才会执行 raw->pkl。
    ThermoRawToPklAgent(
        folder_path=str(COMMON_RAW_DIR), MS2_filter=MS2_filter, max_workers=raw_workers,
        recursive=0, output_folder=str(COMMON_OUT_DIR),
    ).run()

    EQ3_batch(
        cpp_file_input=CPP_FILE,
        ppm_err1_input=ppm_err1_input, da_err1_input=da_err1_input,
        ppm_err2_input=ppm_err2_input, da_err2_input=da_err2_input,
        sample_inputs=find_ms1_msms_files(str(COMMON_OUT_DIR)),
        library_input=str(COMMON_LIBRARY), max_workers=EQ3_WORKERS,
    )

    retention_time_inference(INPUT_CSV=common_aligned_table, overwrite=False)

    build_master_table(
        pkl_dir=str(COMMON_OUT_DIR), rt_checked_csv=common_rt_checked_table,
        library_pkl=str(COMMON_LIBRARY), ppm_err1=ppm_err1_input, da_err1=da_err1_input,
        ms2_only=MASTER_MS2_ONLY,
    )

    if STAGE1_QUANTIFY:
        quantify_master_table(
            master_table_csv=common_master_table, raw_dir=str(COMMON_RAW_DIR), max_workers=None,
        )


def stage2_pb_search():
    """第二阶段：用常规注释主表过滤 PB 库后，对 sample/PB 做 PB 检索。"""
    print("\n" + "=" * 60)
    print("第二阶段 | PB 衍生化检索 sample/PB")
    print("=" * 60)

    if not os.path.exists(common_master_table):
        raise FileNotFoundError(
            f"未找到第一阶段主表 {common_master_table}，请先运行第一阶段(RUN_STAGE1=True)"
        )

    n_raw = len(require_raw(PB_RAW_DIR, "PB 衍生化检索"))
    print(f"[信息] raw {n_raw} 个：{PB_RAW_DIR} -> pkl：{PB_OUT_DIR}")

    run_pb_search(
        pb_raw_dir=str(PB_RAW_DIR),
        pb_library_pkl=str(PB_LIBRARY),
        cpp_file=CPP_FILE,
        ms2_filter=MS2_filter,
        ppm_err1=pb_ppm_err1, da_err1=pb_da_err1,
        ppm_err2=ppm_err2_input, da_err2=pb_da_err2,
        eq3_workers=EQ3_WORKERS, raw_workers=raw_workers,
        pb_output_dir=str(PB_OUT_DIR),
    )


def stage3_pb_localize():
    """第三阶段：在 *_PB_diagnostic.csv 上做双键位置重建，输出 *_PB_isomers.csv。"""
    print("\n" + "=" * 60)
    print("第三阶段 | PB 双键位置重建 sample/PB")
    print("=" * 60)

    if not os.path.exists(common_master_quant):
        raise FileNotFoundError(
            f"未找到第一阶段定量主表 {common_master_quant}，请先运行第一阶段(STAGE1_QUANTIFY=True)"
        )

    localize_all_pb_diagnostic(
        pb_output_dir=str(PB_OUT_DIR),
        master_quant_csv=common_master_quant,
        pb_library_pkl=str(PB_LIBRARY),
        l12_min_coverage=PB_LOC_L12_MIN_COV,
        ppm_err2=ppm_err2_input, da_err2=pb_da_err2,
        rt_tol=PB_LOC_RT_TOL, isomer_top_n=PB_LOC_ISOMER_TOP_N,
    )


def stage4_pb_quantify():
    """第四阶段：在 *_PB_dbpos.csv 上定量到 (物质,链,双键位置)，出 PB_dbposition_quant.csv。"""
    print("\n" + "=" * 60)
    print("第四阶段 | PB 链×双键位置级定量 sample/PB")
    print("=" * 60)

    quantify_pb_dbpositions(
        pb_output_dir=str(PB_OUT_DIR),
        pb_raw_dir=str(PB_RAW_DIR),
        rt_tol=PB_QUANT_RT_CLUSTER_TOL,
        max_workers=raw_workers,
    )


if __name__ == "__main__":
    start = time.perf_counter()
    check_inputs()
    if RUN_STAGE1:
        stage1_common_annotation()
    if RUN_STAGE2:
        stage2_pb_search()
    if RUN_STAGE3:
        stage3_pb_localize()
    if RUN_STAGE4:
        stage4_pb_quantify()
    print(f"\n全流程耗时: {time.perf_counter() - start:.4f} 秒")
