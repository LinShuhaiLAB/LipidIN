import pickle
import pandas as pd
import numpy as np
import time
import os
import re
import gc
import concurrent.futures
from bisect import bisect_left
from functools import lru_cache
from pathlib import Path
from processEQ2v2 import *
def compile_cpp_module(cpp_file_input):
    try:
        import eq3_search
        return True
    except ImportError:
        print("Cpp modules need to be compiled...\nCpp 模块需要编译...")
        if not os.path.exists(cpp_file_input):
            print(f"✗ Cpp source file not found: {cpp_file_input}\n✗ 未找到 Cpp 源文件: {cpp_file_input}")
            return False
        try:
            from compile_cpp_fixed import compile_with_fixed_path
            success = compile_with_fixed_path()
            if not success:
                return False
            import eq3_search
            return True
        except Exception as e:
            print(f"✗ Cpp module compilation failed: {e}\n✗ Cpp 模块编译失败: {e}")
            return False
def load_and_prepare_data(sample_input, library_input):
    library_start_time = time.time()
    print("=" * 60)
    print("Load Dynamic Hierarchical Spectral Library...\n加载动态分层谱库...")
    pkl_start_time = time.time()
    with open(sample_input, 'rb') as f:
        standard_data = pickle.load(f)
    with open(library_input, 'rb') as f:
        library_df = pickle.load(f)
    pkl_time = time.time() - pkl_start_time
    standard_data_sorted = sorted(standard_data, key=lambda x: x[1])
    library_time = time.time() - library_start_time
    print(f"✓ Sample spectrum: {len(standard_data_sorted)}\n✓ 样本谱图: {len(standard_data_sorted)}")
    print(f"✓ Database records: {len(library_df)}\n✓ 数据库记录数: {len(library_df)}")
    print(f"✓ Load pkl data time: {pkl_time:.4f} s\n✓ 加载 pkl 数据耗时: {pkl_time:.4f} s")
    print(f"✓ Load Dynamic Hierarchical Spectral Library time: {library_time:.4f} s\n✓ 加载动态分层谱库耗时: {library_time:.4f} s")
    return standard_data_sorted, library_df
def prepare_cpp_data(library_df):
    db_pmz = library_df['pmz'].to_numpy(dtype=np.float64)
    db_ms2_mz = library_df['mz_value'].to_numpy(dtype=np.float64)
    db_ms2_int = library_df['ms2_int'].to_numpy(dtype=np.float64)
    db_labels = library_df['lab_mz'].tolist()
    return (
        db_pmz.tolist(),
        db_ms2_mz.tolist(),
        db_ms2_int.tolist(),
        db_labels
    )
def prepare_sample_data(standard_data_sorted):
    peak_ids = []
    sample_pmz_list = []
    sample_rt_list = []
    sample_ms2_list = []
    sample_intensity_list = []
    for spectrum in standard_data_sorted:
        peak_ids.append(int(spectrum[0]))
        sample_pmz_list.append(float(spectrum[1]))
        sample_rt_list.append(float(spectrum[2]))
        sample_ms2_list.append([float(x) for x in spectrum[3]])
        sample_intensity_list.append([float(x) for x in spectrum[4]])
    return peak_ids, sample_pmz_list, sample_rt_list, sample_ms2_list, sample_intensity_list
def initialize_cpp_search_engine(db_pmz, db_ms2_mz, db_ms2_int, db_labels):
    import eq3_search
    print("=" * 60)
    search_engine = eq3_search.FastSpectralSearchEQ4()
    search_engine.load_database(db_pmz, db_ms2_mz, db_labels)
    return search_engine
def run_single_search(search_engine, peak_ids, sample_pmz_list, sample_rt_list,
                      sample_ms2_list, sample_intensity_list,
                      ppm_err1=20, da_err1=0.01, ppm_err2=20, da_err2=0.01):
    start_time = time.time()
    matches = search_engine.search_batch_eq4(
        peak_ids,
        sample_pmz_list,
        sample_rt_list,
        sample_ms2_list,
        sample_intensity_list,
        ppm_err1,
        da_err1,
        ppm_err2,
        da_err2
    )
    search_time = time.time() - start_time
    print(f"✓ Library search time: {search_time:.4f} s\n✓ 谱库搜索耗时: {search_time:.4f} s")
    return matches
def convert_matches_to_dataframe(matches, search_engine=None):
    output_columns = [
        'peak_id',
        'label',
        'sample_pmz',
        'sample_rt',
        'sample_pmz_da_error',
        'sample_ms2_mz',
        'sample_ms2_int'
    ]
    if not matches:
        return pd.DataFrame(columns=output_columns)

    n_match = len(matches)
    peak_id = np.empty(n_match, dtype=np.int64)
    da_error = np.empty(n_match, dtype=np.float64)
    label = [None] * n_match

    for i, match in enumerate(matches):
        peak_id[i] = match.peak_id
        da_error[i] = match.sample_pmz_da_error
        label[i] = match.label

    
    
    unique_peak, first_row = np.unique(peak_id, return_index=True)
    codes = np.searchsorted(unique_peak, peak_id).tolist()

    pmz_by_peak = np.array([round(matches[i].sample_pmz, 4) for i in first_row], dtype=np.float64)
    rt_by_peak = np.array([matches[i].sample_rt for i in first_row], dtype=np.float64)
    ms2_mz_by_peak = [matches[i].sample_ms2_mz for i in first_row]
    ms2_int_by_peak = [matches[i].sample_ms2_int for i in first_row]

    results_df = pd.DataFrame({
        'peak_id': peak_id,
        'label': label,
        'sample_pmz': pmz_by_peak[codes],
        'sample_rt': rt_by_peak[codes],
        'sample_pmz_da_error': da_error,
        'sample_ms2_mz': [ms2_mz_by_peak[c] for c in codes],
        'sample_ms2_int': [ms2_int_by_peak[c] for c in codes]
    })
    results_df = results_df[output_columns]
    
    
    results_df.drop_duplicates(
        subset=['peak_id', 'label'],
        keep='first',
        inplace=True
    )
    results_df.sort_values(
        ['peak_id', 'label'],
        inplace=True,
        kind='mergesort'
    )
    results_df.reset_index(drop=True, inplace=True)
    print(f"✓ EQ3 matches retained: {len(results_df)}\n✓ 保留的 EQ3 匹配数: {len(results_df)}")
    return results_df
def _to_excel_text(value):
    if pd.isna(value):
        return value

    text = str(value)

    if text == "":
        return text

    return f'="{text}"'
def _from_excel_text(value):
    if pd.isna(value):
        return value

    text = str(value)

    if text.startswith('="') and text.endswith('"'):
        return text[2:-1]

    return text
def save_results(results_df, sample_input, filename=None, ppm_err2=20, da_err2=0.02):
    if results_df is None or results_df.empty:
        print("No repeated-label results need to be saved.\n没有需要保存的重复标签结果。")
        return

    results_df = split_label_columns_before_output(results_df)
    results_df = annotate_subclass_and_fragment_checks(results_df, ppm_err2, da_err2)

    if filename is None:
        filename = f'{sample_input}_EQ2.csv'

    start_time = time.time()

    
    results_df = results_df.copy()
    results_df["fa_composition"] = results_df["fa_composition"].map(_to_excel_text)

    
    results_df.to_csv(filename, index=False, encoding="utf-8-sig")

    output_time = time.time() - start_time

    print(f"✓ Result saved: {filename}\n✓ 结果已保存: {filename}")
    print(f"✓ CSV output time: {output_time:.4f} s\n✓ CSV 输出耗时: {output_time:.4f} s")
_EQ3_ENGINES = {}


def _eq3_worker_init(library_input, cpp_file_input):
    compile_cpp_module(cpp_file_input)
    compile_cpp_similarity_module("eq3_similarity.cpp")
    search_engine, _, similarity_engine = load_library_once(library_input, quiet=True)
    _EQ3_ENGINES["search"] = search_engine
    _EQ3_ENGINES["similarity"] = similarity_engine


def _eq3_process_one_sample(task_args):
    import io
    from contextlib import redirect_stdout

    sample_input, ppm_err1, da_err1, ppm_err2, da_err2 = task_args
    buffer = io.StringIO()

    try:
        with redirect_stdout(buffer):
            results_df = run_one_eq3_sample(
                _EQ3_ENGINES["search"],
                _EQ3_ENGINES["similarity"],
                sample_input,
                ppm_err1,
                da_err1,
                ppm_err2,
                da_err2,
            )
    except Exception as e:
        import traceback
        return sample_input, None, buffer.getvalue() + "\n" + traceback.format_exc()

    return sample_input, results_df, buffer.getvalue()


def run_one_eq3_sample(search_engine, similarity_engine, sample_input,
                       ppm_err1_input, da_err1_input, ppm_err2_input, da_err2_input):
    print("=" * 60)
    print(f"Process sample: {sample_input}\n处理样本: {sample_input}")

    sample_total_start_time = time.time()

    pkl_start_time = time.time()
    with open(sample_input, 'rb') as f:
        standard_data = pickle.load(f)

    standard_data_sorted = sorted(standard_data, key=lambda x: x[1])

    print(f"✓ Load sample pkl data time: {time.time() - pkl_start_time:.4f} s\n✓ 加载样本 pkl 数据耗时: {time.time() - pkl_start_time:.4f} s")
    print(f"✓ Sample spectrum: {len(standard_data_sorted)}\n✓ 样本谱图: {len(standard_data_sorted)}")

    prepare_start_time = time.time()
    peak_ids, sample_pmz_list, sample_rt_list, sample_ms2_list, sample_intensity_list = prepare_sample_data(
        standard_data_sorted
    )
    print(f"✓ Prepare sample data time: {time.time() - prepare_start_time:.4f} s\n✓ 准备样本数据耗时: {time.time() - prepare_start_time:.4f} s")

    matches = run_single_search(
        search_engine,
        peak_ids,
        sample_pmz_list,
        sample_rt_list,
        sample_ms2_list,
        sample_intensity_list,
        ppm_err1_input,
        da_err1_input,
        ppm_err2_input,
        da_err2_input
    )

    results_df = convert_matches_to_dataframe(matches, search_engine)

    
    del matches
    gc.collect()

    results_df = postprocess_results_before_output(
        results_df,
        similarity_engine,
        ppm_err2=ppm_err2_input,
        da_err2=da_err2_input,
        jcr_cutoff=0.2
    )

    annotated_peak_ids = set()
    if results_df is not None and not results_df.empty and "peak_id" in results_df.columns:
        annotated_peak_ids = set(results_df["peak_id"].astype(int).tolist())

    unannotated_records = []
    for peak_id, sample_pmz, sample_rt, sample_ms2_mz, sample_ms2_int in zip(
        peak_ids,
        sample_pmz_list,
        sample_rt_list,
        sample_ms2_list,
        sample_intensity_list
    ):
        if int(peak_id) not in annotated_peak_ids:
            unannotated_records.append({
                "peak_id": int(peak_id),
                "sample_pmz": float(sample_pmz),
                "sample_rt": float(sample_rt),
                "sample_ms2_mz": "_".join(map(str, sample_ms2_mz)),
                "sample_ms2_int": "_".join(map(str, sample_ms2_int))
            })

    unannotated_df = pd.DataFrame(
        unannotated_records,
        columns=[
            "peak_id",
            "sample_pmz",
            "sample_rt",
            "sample_ms2_mz",
            "sample_ms2_int"
        ]
    )

    unannotated_output = f"{sample_input.replace('.pkl', '')}_unannotated.csv"
    unannotated_df.to_csv(unannotated_output, index=False)
    print(f"✓ Unannotated result saved: {unannotated_output}\n✓ 未注释结果已保存: {unannotated_output}")
    print(f"✓ Unannotated peaks: {len(unannotated_df)}\n✓ 未注释峰数: {len(unannotated_df)}")

    save_results(
        results_df,
        sample_input.replace('.pkl', ''),
        ppm_err2=ppm_err2_input,
        da_err2=da_err2_input
    )

    print(f"✓ Sample total processing time: {time.time() - sample_total_start_time:.4f} s\n✓ 样本总处理耗时: {time.time() - sample_total_start_time:.4f} s")

    return results_df




EQ3_BYTES_PER_WORKER = 2.2e9
EQ3_MAX_WORKERS = 4


def resolve_eq3_workers(n_samples, max_workers=None, bytes_per_worker=EQ3_BYTES_PER_WORKER):
    if max_workers is None:
        gc.collect()  
        try:
            import psutil
            available = psutil.virtual_memory().available
        except Exception:
            available = bytes_per_worker  
        max_workers = int((available * 0.85) // bytes_per_worker)

    cpu_limit = max(1, (os.cpu_count() or 2) - 1)
    return max(1, min(int(max_workers), n_samples, cpu_limit, EQ3_MAX_WORKERS))


def EQ3_batch(cpp_file_input, ppm_err1_input, da_err1_input,
              ppm_err2_input, da_err2_input,
              sample_inputs, library_input, max_workers=None):
    if isinstance(sample_inputs, (str, Path)):
        raise TypeError(
            f"sample_inputs 需要传入 pkl 文件路径的列表，不能传单个路径字符串，"
            f"否则会被逐字符遍历。当前收到：{sample_inputs!r}"
        )
    print("🚀 Start batch EQ3 (LipidIN 2.0)\n🚀 开始批量 EQ3 (LipidIN 2.0)")
    print("=" * 60)
    print(f"C++ source path: {cpp_file_input}\nC++ 源码路径: {cpp_file_input}")

    if not compile_cpp_module(cpp_file_input):
        print("❌ Unable to use Cpp optimization, the program exits.\n❌ 无法使用 Cpp 优化，程序退出。")
        return None

    if not compile_cpp_similarity_module("eq3_similarity.cpp"):
        print("❌ Unable to use Cpp MS2 similarity optimization, the program exits.\n❌ 无法使用 Cpp MS2 相似度优化，程序退出。")
        return None

    try:
        workers = resolve_eq3_workers(len(sample_inputs), max_workers)
        all_results = {}

        if workers == 1:
            search_engine, library_df, similarity_engine = load_library_once(library_input)
            for sample_input in sample_inputs:
                all_results[sample_input] = run_one_eq3_sample(
                    search_engine, similarity_engine, sample_input,
                    ppm_err1_input, da_err1_input, ppm_err2_input, da_err2_input
                )
        else:
            print("=" * 60)
            print(f"Process {len(sample_inputs)} samples with {workers} worker processes\n使用 {workers} 个工作进程处理 {len(sample_inputs)} 个样本")

            tasks = [
                (sample_input, ppm_err1_input, da_err1_input, ppm_err2_input, da_err2_input)
                for sample_input in sample_inputs
            ]

            failures = []
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=workers,
                initializer=_eq3_worker_init,
                initargs=(library_input, cpp_file_input),
            ) as executor:
                
                for sample_input, results_df, log_text in executor.map(_eq3_process_one_sample, tasks):
                    print(log_text, end="")
                    if results_df is None:
                        failures.append(sample_input)
                    all_results[sample_input] = results_df

            if failures:
                raise RuntimeError(f"以下样本处理失败：{failures}")

        print("=" * 60)
        print("🎉 Batch search complete!\n🎉 批量搜索完成！")
        merge_root_dir = Path(sample_inputs[0]).parent

        merge_eq2_to_rt_ccs_aligned_table(
            root_dir=merge_root_dir,
            file_suffix="_EQ2.csv",
            mz_tol=da_err1_input,
            rt_tol=0.2,
            output_file="rt_ccs_aligned_table.csv"
        )
        return all_results

    except Exception as e:
        print(f"❌ Errors in batch program execution: {str(e)}\n❌ 批量程序执行出错: {str(e)}")
        import traceback
        traceback.print_exc()
        return None
def load_library_once(library_input, quiet=False):
    def emit(message):
        if not quiet:
            print(message)

    emit("=" * 60)
    emit("Load Dynamic Hierarchical Spectral Library...\n加载动态分层谱库...")

    total_start_time = time.time()

    pkl_start_time = time.time()
    with open(library_input, 'rb') as f:
        library_df = pickle.load(f)
    emit(f"✓ Load library pkl data time: {time.time() - pkl_start_time:.4f} s\n✓ 加载谱库 pkl 数据耗时: {time.time() - pkl_start_time:.4f} s")

    prepare_start_time = time.time()
    db_pmz, db_ms2_mz, db_ms2_int, db_labels = prepare_cpp_data(library_df)
    emit(f"✓ Prepare library data time: {time.time() - prepare_start_time:.4f} s\n✓ 准备谱库数据耗时: {time.time() - prepare_start_time:.4f} s")

    engine_start_time = time.time()
    search_engine = initialize_cpp_search_engine(
        db_pmz,
        db_ms2_mz,
        db_ms2_int,
        db_labels
    )
    emit(f"✓ Load database into search engine time: {time.time() - engine_start_time:.4f} s\n✓ 将数据库载入搜索引擎耗时: {time.time() - engine_start_time:.4f} s")

    
    del db_pmz, db_ms2_mz, db_ms2_int, db_labels
    gc.collect()

    similarity_start_time = time.time()
    similarity_engine = initialize_cpp_similarity_engine(library_df)
    emit(f"✓ Load database into similarity engine time: {time.time() - similarity_start_time:.4f} s\n✓ 将数据库载入相似度引擎耗时: {time.time() - similarity_start_time:.4f} s")

    emit(f"✓ Database records: {len(library_df)}\n✓ 数据库记录数: {len(library_df)}")
    emit(f"✓ Load Dynamic Hierarchical Spectral Library time: {time.time() - total_start_time:.4f} s\n✓ 加载动态分层谱库耗时: {time.time() - total_start_time:.4f} s")

    if quiet:
        del library_df
        gc.collect()
        return search_engine, None, similarity_engine

    return search_engine, library_df, similarity_engine
def strip_dtl_prefix(label):
    text = str(label)
    if "_DTL_" in text:
        return text.split("_DTL_", 1)[1]
    return text
def postprocess_results_before_output(results_df, similarity_engine,
                                      ppm_err2=20, da_err2=0.01,
                                      jcr_cutoff=0.4):
    if results_df is None or results_df.empty:
        return results_df

    start_time = time.time()

    results_df = results_df.copy()
    results_df["label"] = results_df["label"].map(strip_dtl_prefix)

    
    
    
    
    label_subclass = results_df["label"].astype(str).str.split(" ", n=1).str[0]
    single_peak_mask = label_subclass.isin(SINGLE_PEAK_SUBCLASSES)
    repeated_mask = results_df.duplicated(subset=["peak_id", "label"], keep=False)

    results_df = results_df.loc[repeated_mask | single_peak_mask].copy()

    if results_df.empty:
        print("✓ Results retained after repeated-label filtering: 0\n✓ 重复标签过滤后保留的结果数: 0")
        return results_df

    results_df.drop_duplicates(
        subset=["peak_id", "label"],
        keep="first",
        inplace=True
    )
    results_df.reset_index(drop=True, inplace=True)

    similarity_results = similarity_engine.calculate_batch(
        results_df["label"].astype(str).tolist(),
        results_df["sample_pmz"].astype(np.float64).tolist(),
        results_df["sample_pmz_da_error"].astype(np.float64).tolist(),
        results_df["sample_ms2_mz"].astype(str).tolist(),
        results_df["sample_ms2_int"].astype(str).tolist(),
        float(ppm_err2),
        float(da_err2),
        float(jcr_cutoff)
    )

    results_df["jcr_similarity"] = [
        item.jcr_similarity for item in similarity_results
    ]
    results_df["cosine_similarity"] = [
        item.cosine_similarity for item in similarity_results
    ]
    results_df["intersection_over_union"] = [
        item.intersection_over_union for item in similarity_results
    ]

    n_before_jcr_filter = len(results_df)

    results_df["jcr_similarity"] = pd.to_numeric(
        results_df["jcr_similarity"],
        errors="coerce"
    )

    results_df = results_df.loc[
        results_df["jcr_similarity"] > float(jcr_cutoff)
    ].copy()

    results_df.reset_index(drop=True, inplace=True)

    n_removed_by_jcr = n_before_jcr_filter - len(results_df)

    print(f"✓ Results retained after repeated-label deduplication: {n_before_jcr_filter}\n✓ 重复标签去重后保留的结果数: {n_before_jcr_filter}")
    print(f"✓ Results removed by JCR <= {jcr_cutoff}: {n_removed_by_jcr}\n✓ 被 JCR <= {jcr_cutoff} 移除的结果数: {n_removed_by_jcr}")
    print(f"✓ Results retained after JCR filtering: {len(results_df)}\n✓ JCR 过滤后保留的结果数: {len(results_df)}")
    print(f"✓ C++ MS2 similarity calculation time: {time.time() - start_time:.4f} s\n✓ C++ MS2 相似度计算耗时: {time.time() - start_time:.4f} s")

    return results_df
def merge_eq2_to_rt_ccs_aligned_table(
    root_dir,
    file_suffix="_EQ2.csv",
    mz_col="sample_pmz",
    rt_col="sample_rt",
    name_col="total_chain_name",
    sub_name_col="sub_chain_name",
    adduct_col="adduct_form",
    source_ms2mz_col="sample_ms2_mz",
    source_ms2int_col="sample_ms2_int",
    output_ms2mz_col="ms2mz",
    output_ms2int_col="ms2intensity",
    mz_tol=0.01,
    rt_tol=0.2,
    output_file="rt_ccs_aligned_table.csv"
):
    root_dir = Path(root_dir)

    similarity_cols = [
        "jcr_similarity",
        "cosine_similarity",
        "intersection_over_union"
    ]

    csv_files = [
        p for p in root_dir.rglob("*.csv")
        if p.name.endswith(file_suffix)
    ]

    if len(csv_files) == 0:
        raise FileNotFoundError(f"No files ending with {file_suffix} found under {root_dir}")

    def _clean_ms2_text(value):
        if pd.isna(value):
            return ""

        text = str(value).strip()

        if text == "" or text.lower() == "nan":
            return ""

        return text

    def _split_unique_values(values, sep=";"):
        out = []
        seen = set()

        for value in values:
            if pd.isna(value):
                continue

            for item in str(value).split(sep):
                item = item.strip()

                if item and item not in seen:
                    seen.add(item)
                    out.append(item)

        return out

    
    sub_name_sep = "|"

    def _add_ms2_to_peak(peak, row):
        ms2mz_text = _clean_ms2_text(row.get(source_ms2mz_col, np.nan))
        ms2int_text = _clean_ms2_text(row.get(source_ms2int_col, np.nan))

        if ms2mz_text == "" and ms2int_text == "":
            return

        pair = (ms2mz_text, ms2int_text)

        if pair not in peak["ms2_pair_seen"]:
            peak["ms2_pair_seen"].add(pair)
            peak["ms2_pairs"].append(pair)

    def _align_table(df, id_prefix):
        aligned_rows = []
        peak_index = 0

        df = df.sort_values(
            by=[name_col, rt_col, mz_col],
            ascending=[True, True, True]
        ).reset_index(drop=True)

        for lipid_name, group in df.groupby(name_col, sort=False):
            group = group.reset_index(drop=True)

            peaks = []
            rt_bins = {}

            for _, row in group.iterrows():
                mz = row[mz_col]
                rt = row[rt_col]
                rt_bin = int(np.floor(rt / rt_tol))

                candidate_peak_ids = set()

                for b in range(rt_bin - 1, rt_bin + 2):
                    candidate_peak_ids.update(rt_bins.get(b, set()))

                candidates = []

                for pid in candidate_peak_ids:
                    peak = peaks[pid]

                    mz_err = abs(mz - peak["aligned_mz"])
                    rt_err = abs(rt - peak["aligned_rt"])

                    if mz_err <= mz_tol and rt_err <= rt_tol:
                        candidates.append((pid, rt_err, mz_err))

                if len(candidates) > 0:
                    candidates.sort(key=lambda x: (x[1], x[2]))
                    best_pid = candidates[0][0]
                    peak = peaks[best_pid]

                    old_bin = peak["rt_bin"]

                    row_n = int(row.get("n_records", 1))
                    old_n = peak["n_records"]
                    new_n = old_n + row_n

                    peak["mz_sum"] += mz * row_n
                    peak["rt_sum"] += rt * row_n
                    peak["n_records"] = new_n
                    peak["aligned_mz"] = peak["mz_sum"] / new_n
                    peak["aligned_rt"] = peak["rt_sum"] / new_n

                    peak["mz_min"] = min(peak["mz_min"], row.get("mz_min", mz))
                    peak["mz_max"] = max(peak["mz_max"], row.get("mz_max", mz))
                    peak["rt_min"] = min(peak["rt_min"], row.get("rt_min", rt))
                    peak["rt_max"] = max(peak["rt_max"], row.get("rt_max", rt))

                    for adduct in _split_unique_values([row.get(adduct_col, np.nan)]):
                        peak["adduct_forms"].add(adduct)

                    for sub_name in _split_unique_values([row.get(sub_name_col, np.nan)], sep=sub_name_sep):
                        peak["sub_chain_names"].add(sub_name)

                    _add_ms2_to_peak(peak, row)

                    for sim_col in similarity_cols:
                        value = row.get(sim_col, np.nan)

                        if pd.notna(value):
                            if pd.isna(peak[sim_col]):
                                peak[sim_col] = value
                            else:
                                peak[sim_col] = max(peak[sim_col], value)

                    new_bin = int(np.floor(peak["aligned_rt"] / rt_tol))

                    if new_bin != old_bin:
                        rt_bins[old_bin].discard(best_pid)
                        rt_bins.setdefault(new_bin, set()).add(best_pid)
                        peak["rt_bin"] = new_bin

                else:
                    pid = len(peaks)
                    row_n = int(row.get("n_records", 1))

                    peak = {
                        name_col: lipid_name,
                        "aligned_mz": mz,
                        "aligned_rt": rt,
                        "mz_sum": mz * row_n,
                        "rt_sum": rt * row_n,
                        "n_records": row_n,
                        "mz_min": row.get("mz_min", mz),
                        "mz_max": row.get("mz_max", mz),
                        "rt_min": row.get("rt_min", rt),
                        "rt_max": row.get("rt_max", rt),
                        "adduct_forms": set(_split_unique_values([row.get(adduct_col, np.nan)])),
                        "sub_chain_names": set(_split_unique_values([row.get(sub_name_col, np.nan)], sep=sub_name_sep)),
                        "ms2_pairs": [],
                        "ms2_pair_seen": set(),
                        "rt_bin": rt_bin
                    }

                    _add_ms2_to_peak(peak, row)

                    for sim_col in similarity_cols:
                        peak[sim_col] = row.get(sim_col, np.nan)

                    peaks.append(peak)
                    rt_bins.setdefault(rt_bin, set()).add(pid)

            for peak in peaks:
                peak_index += 1

                output_row = {
                    "peak_id": f"{id_prefix}_{peak_index}",
                    name_col: peak[name_col],
                    sub_name_col: sub_name_sep.join(sorted(peak["sub_chain_names"])),
                    "aligned_mz": peak["aligned_mz"],
                    "aligned_rt": peak["aligned_rt"],
                    "n_records": peak["n_records"],
                    "mz_min": peak["mz_min"],
                    "mz_max": peak["mz_max"],
                    "rt_min": peak["rt_min"],
                    "rt_max": peak["rt_max"],
                    adduct_col: ";".join(sorted(peak["adduct_forms"])),
                    output_ms2mz_col: "||".join([pair[0] for pair in peak["ms2_pairs"]]),
                    output_ms2int_col: "||".join([pair[1] for pair in peak["ms2_pairs"]])
                }

                for sim_col in similarity_cols:
                    output_row[sim_col] = peak[sim_col]

                aligned_rows.append(output_row)

        return pd.DataFrame(aligned_rows)

    per_file_aligned = []

    for csv_file in csv_files:
        df = pd.read_csv(csv_file, encoding="utf-8-sig")

        if "fa_composition" in df.columns:
            df["fa_composition"] = df["fa_composition"].map(_from_excel_text)

        required_cols = [
            mz_col,
            rt_col,
            name_col,
            adduct_col,
            source_ms2mz_col,
            source_ms2int_col
        ]

        missing_cols = [
            col for col in required_cols
            if col not in df.columns
        ]

        if missing_cols:
            raise ValueError(f"{csv_file} missing required columns: {missing_cols}")

        if sub_name_col not in df.columns:
            print(f"⚠ {csv_file.name} is missing the {sub_name_col} column; sub-chain names for this file will be left empty\n⚠ {csv_file.name} 缺少 {sub_name_col} 列，该文件的亚链名将留空")
            df[sub_name_col] = ""

        check_cols = ["headgroup_check", "fa_chain_check"]
        missing_check_cols = [col for col in check_cols if col not in df.columns]

        if missing_check_cols:
            print(f"⚠ {csv_file.name} is missing {missing_check_cols}; skipping headgroup/fatty-acid-chain filtering\n⚠ {csv_file.name} 缺少 {missing_check_cols}，跳过头基/脂肪链过滤")
        else:
            n_before_check_filter = len(df)

            df = df.loc[
                df["headgroup_check"].isin(CHECK_ACCEPTED)
                & df["fa_chain_check"].isin(CHECK_ACCEPTED)
            ].copy()

            print(
                f"✓ {csv_file.name} headgroup/fatty-acid-chain filtering: {n_before_check_filter} -> {len(df)}\n✓ {csv_file.name} 头基/脂肪链过滤：{n_before_check_filter} -> {len(df)}"
            )

        for col in [mz_col, rt_col] + similarity_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        for sim_col in similarity_cols:
            if sim_col not in df.columns:
                df[sim_col] = np.nan

        df = df.dropna(
            subset=[
                mz_col,
                rt_col,
                name_col
            ]
        ).copy()

        if df.empty:
            print(f"⚠ {csv_file.name} has no remaining records after filtering; skipping\n⚠ {csv_file.name} 过滤后无剩余记录，跳过")
            continue

        df["n_records"] = 1
        df["mz_min"] = df[mz_col]
        df["mz_max"] = df[mz_col]
        df["rt_min"] = df[rt_col]
        df["rt_max"] = df[rt_col]

        file_result = _align_table(
            df=df,
            id_prefix=csv_file.stem
        )

        per_file_aligned.append(file_result)

    if not per_file_aligned:
        raise ValueError(f"所有 {file_suffix} 文件在头基/脂肪链过滤后均无剩余记录")

    data = pd.concat(per_file_aligned, ignore_index=True)

    data = data.rename(
        columns={
            "aligned_mz": mz_col,
            "aligned_rt": rt_col,
            output_ms2mz_col: source_ms2mz_col,
            output_ms2int_col: source_ms2int_col
        }
    )

    result = _align_table(
        df=data,
        id_prefix="peak"
    )

    result = result.sort_values(
        by=[
            name_col,
            "aligned_rt",
            "aligned_mz"
        ],
        ascending=[
            True,
            True,
            True
        ]
    ).reset_index(drop=True)

    result["peak_id"] = [
        f"peak_{i + 1}"
        for i in range(len(result))
    ]

    first_columns = [
        "peak_id",
        "aligned_mz",
        "aligned_rt",
        output_ms2mz_col,
        output_ms2int_col,
        name_col,
        sub_name_col,
        adduct_col,
        "n_records",
        "mz_min",
        "mz_max",
        "rt_min",
        "rt_max",
        "jcr_similarity",
        "cosine_similarity",
        "intersection_over_union"
    ]

    first_columns = [
        col for col in first_columns
        if col in result.columns
    ]

    remaining_columns = [
        col for col in result.columns
        if col not in first_columns
    ]

    result = result[first_columns + remaining_columns]

    output_path = root_dir / output_file
    result.to_csv(output_path, index=False, encoding="utf-8-sig")

    print(f"✓ Final aligned table saved: {output_path}\n✓ 最终对齐表已保存: {output_path}")
    print(f"✓ Final aligned records: {len(result)}\n✓ 最终对齐记录数: {len(result)}")

    return result
def split_label_columns_before_output(results_df):
    if results_df is None or results_df.empty:
        return results_df

    results_df = results_df.copy()

    split_pattern = r"_CAH_|_ADD_|_FML_"

    label_parts = results_df["label"].astype(str).str.split(
        split_pattern,
        n=3,
        expand=True
    )

    results_df["total_chain_name"] = label_parts[0]
    results_df["sub_chain_name"] = label_parts[1]
    results_df["adduct_form"] = label_parts[2]
    results_df["chemical_formula"] = label_parts[3]

    output_columns = [
        "peak_id",
        "sample_pmz",
        "sample_rt",
        "total_chain_name",
        "sub_chain_name",
        "adduct_form",
        "chemical_formula",
        "sample_pmz_da_error",
        "jcr_similarity",
        "cosine_similarity",
        "intersection_over_union",
        "sample_ms2_mz",
        "sample_ms2_int"
    ]

    output_columns = [
        col for col in output_columns
        if col in results_df.columns
    ]

    return results_df[output_columns]


CHECK_PASS = "通过"
CHECK_MISSING = "缺失"
CHECK_NOT_APPLICABLE = "无需判断"
CHECK_SKIPPED_BY_IOU = "无需搜索(IoU=1)"
CHECK_SKIPPED_BY_TG_IOU = "无需搜索(TG IoU>0.75)"
CHECK_TG_IOU_TOO_LOW = "TG IoU<=0.75"


CHECK_ACCEPTED = (CHECK_PASS, CHECK_SKIPPED_BY_IOU, CHECK_SKIPPED_BY_TG_IOU)


TG_SUBCLASSES = ("TG", "TG-O")
TG_IOU_CUTOFF = 0.75




SINGLE_PEAK_SUBCLASSES = ("FA",)

CODE_DIR = Path(__file__).resolve().parent
SPB_FA_REFERENCE_FILE = CODE_DIR / "SPB_FA_intact_ion_mz.csv"
HG_REFERENCE_FILE = CODE_DIR / "HG.csv"
def _adduct_to_ionmode(adduct_form):
    if pd.isna(adduct_form):
        return ""

    text = str(adduct_form).strip()

    if text.endswith("+"):
        return "pos"

    if text.endswith("-"):
        return "neg"

    return ""
@lru_cache(maxsize=None)
def _load_spb_fa_reference(csv_path):
    df = pd.read_csv(csv_path)
    reference = {}

    for name, mz, ionmode in zip(df["name"], df["mz"], df["ionmode"]):
        if pd.isna(name) or pd.isna(mz) or pd.isna(ionmode):
            continue

        reference[(str(ionmode).strip(), str(name).strip())] = float(mz)

    return reference
@lru_cache(maxsize=None)
def _load_hg_reference(csv_path):
    df = pd.read_csv(csv_path)
    reference = {}

    for name, mz, ionmode in zip(df["name"], df["mz"], df["ionmode"]):
        if pd.isna(name) or pd.isna(mz) or pd.isna(ionmode):
            continue

        subclasses = [item.strip() for item in str(name).split("|") if item.strip()]
        fragment_mzs = [float(item.strip()) for item in str(mz).split("|") if item.strip()]
        mode = str(ionmode).strip()

        for subclass in subclasses:
            reference.setdefault((mode, subclass), []).extend(fragment_mzs)

    return reference
def _parse_sample_ms2_mz(sample_ms2_mz):
    if pd.isna(sample_ms2_mz):
        return []

    values = []

    for item in str(sample_ms2_mz).split("_"):
        item = item.strip()

        if not item:
            continue

        try:
            values.append(float(item))
        except ValueError:
            continue

    values.sort()

    return values
def _has_fragment(sorted_sample_mz, target_mz, ppm_err2, da_err2):
    if not sorted_sample_mz:
        return False

    pos = bisect_left(sorted_sample_mz, target_mz)

    for idx in (pos - 1, pos):
        if idx < 0 or idx >= len(sorted_sample_mz):
            continue

        sample_mz = sorted_sample_mz[idx]
        tol = max(sample_mz * ppm_err2 * 1e-6, da_err2)

        if abs(sample_mz - target_mz) <= tol:
            return True

    return False
def _split_chain_components(fa_composition):
    if not fa_composition:
        return []

    text = str(fa_composition).replace("/", "_")

    return [item.strip() for item in text.split("_") if item.strip()]
def _lookup_chain_mz(reference, ionmode, component):
    for prefix in ("FA", "SPB"):
        mz = reference.get((ionmode, f"{prefix} {component}"))

        if mz is not None:
            return mz

    return None
def _check_fa_chains(fa_composition, ionmode, sorted_sample_mz, reference, ppm_err2, da_err2):
    components = _split_chain_components(fa_composition)

    if not components or not ionmode:
        return CHECK_NOT_APPLICABLE

    known_mzs = [
        mz for mz in (
            _lookup_chain_mz(reference, ionmode, component)
            for component in components
        )
        if mz is not None
    ]

    
    if not known_mzs:
        return CHECK_NOT_APPLICABLE

    
    for mz in known_mzs:
        if _has_fragment(sorted_sample_mz, mz, ppm_err2, da_err2):
            return CHECK_PASS

    return CHECK_MISSING
def _check_headgroup(subclass, ionmode, sorted_sample_mz, reference, ppm_err2, da_err2):
    if not subclass or not ionmode:
        return CHECK_NOT_APPLICABLE

    fragment_mzs = reference.get((ionmode, subclass))

    if not fragment_mzs:
        return CHECK_NOT_APPLICABLE

    
    for mz in fragment_mzs:
        if _has_fragment(sorted_sample_mz, mz, ppm_err2, da_err2):
            return CHECK_PASS

    return CHECK_MISSING
def annotate_subclass_and_fragment_checks(results_df, ppm_err2, da_err2):

    if results_df is None or results_df.empty:
        return results_df

    results_df = results_df.copy()

    name_parts = results_df["sub_chain_name"].astype(str).str.split(" ", n=1, expand=True)

    results_df["subclass"] = name_parts[0].fillna("")
    results_df["fa_composition"] = (
        name_parts[1].fillna("") if name_parts.shape[1] > 1 else ""
    )

    spb_fa_reference = _load_spb_fa_reference(str(SPB_FA_REFERENCE_FILE))
    hg_reference = _load_hg_reference(str(HG_REFERENCE_FILE))

    iou_values = pd.to_numeric(results_df["intersection_over_union"], errors="coerce")

    headgroup_checks = []
    fa_chain_checks = []

    for iou, adduct_form, subclass, fa_composition, sample_ms2_mz in zip(
        iou_values,
        results_df["adduct_form"],
        results_df["subclass"],
        results_df["fa_composition"],
        results_df["sample_ms2_mz"]
    ):
        
        if subclass in TG_SUBCLASSES:
            if pd.notna(iou) and iou > TG_IOU_CUTOFF:
                status = CHECK_SKIPPED_BY_TG_IOU
            else:
                status = CHECK_TG_IOU_TOO_LOW

            headgroup_checks.append(status)
            fa_chain_checks.append(status)
            continue

        
        
        if subclass in SINGLE_PEAK_SUBCLASSES:
            ionmode = _adduct_to_ionmode(adduct_form)
            sorted_sample_mz = _parse_sample_ms2_mz(sample_ms2_mz)
            fa_status = _check_fa_chains(
                fa_composition, ionmode, sorted_sample_mz, spb_fa_reference, ppm_err2, da_err2
            )
            headgroup_checks.append(fa_status)
            fa_chain_checks.append(fa_status)
            continue

        
        if pd.notna(iou) and iou >= 1.0:
            headgroup_checks.append(CHECK_SKIPPED_BY_IOU)
            fa_chain_checks.append(CHECK_SKIPPED_BY_IOU)
            continue

        ionmode = _adduct_to_ionmode(adduct_form)
        sorted_sample_mz = _parse_sample_ms2_mz(sample_ms2_mz)

        headgroup_checks.append(
            _check_headgroup(subclass, ionmode, sorted_sample_mz, hg_reference, ppm_err2, da_err2)
        )
        fa_chain_checks.append(
            _check_fa_chains(fa_composition, ionmode, sorted_sample_mz, spb_fa_reference, ppm_err2, da_err2)
        )

    results_df["headgroup_check"] = headgroup_checks
    results_df["fa_chain_check"] = fa_chain_checks

    insert_after = list(results_df.columns).index("sub_chain_name") + 1
    new_columns = ["subclass", "fa_composition", "headgroup_check", "fa_chain_check"]
    other_columns = [col for col in results_df.columns if col not in new_columns]

    ordered_columns = (
        other_columns[:insert_after] + new_columns + other_columns[insert_after:]
    )

    return results_df[ordered_columns]
def compile_cpp_similarity_module(cpp_similarity_file_input="eq3_similarity.cpp"):
    try:
        import eq3_similarity
        return True
    except ImportError:
        print("Cpp similarity module needs to be compiled...\nCpp 相似度模块需要编译...")

        if not os.path.exists(cpp_similarity_file_input):
            print(f"✗ Cpp similarity source file not found: {cpp_similarity_file_input}\n✗ 未找到 Cpp 相似度源文件: {cpp_similarity_file_input}")
            return False

        try:
            from compile_cpp_similarity import compile_similarity_cpp
            success = compile_similarity_cpp(cpp_similarity_file_input)
            if not success:
                return False

            import eq3_similarity
            return True

        except Exception as e:
            print(f"✗ Cpp similarity module compilation failed: {e}\n✗ Cpp 相似度模块编译失败: {e}")
            return False
def initialize_cpp_similarity_engine(library_df):
    import eq3_similarity

    db_pmz = library_df["pmz"].to_numpy(dtype=np.float64).tolist()
    db_ms2_mz = library_df["mz_value"].to_numpy(dtype=np.float64).tolist()
    db_ms2_int = library_df["ms2_int"].to_numpy(dtype=np.float64).tolist()
    db_labels = library_df["lab_mz"].astype(str).tolist()

    similarity_engine = eq3_similarity.FastMS2Similarity()
    similarity_engine.load_library(
        db_pmz,
        db_ms2_mz,
        db_ms2_int,
        db_labels
    )

    return similarity_engine
def EQ3(cpp_file_input, ppm_err1_input, da_err1_input,
        ppm_err2_input, da_err2_input,
        sample_input, library_input):
    print("🚀 Start EQ3 (LipidIN 2.0)\n🚀 开始 EQ3 (LipidIN 2.0)")
    print("=" * 60)
    print(f"C++ source path: {cpp_file_input}\nC++ 源码路径: {cpp_file_input}")

    if not compile_cpp_module(cpp_file_input):
        print("❌ Unable to use Cpp optimization, the program exits.\n❌ 无法使用 Cpp 优化，程序退出。")
        return None

    if not compile_cpp_similarity_module("eq3_similarity.cpp"):
        print("❌ Unable to use Cpp MS2 similarity optimization, the program exits.\n❌ 无法使用 Cpp MS2 相似度优化，程序退出。")
        return None

    try:
        standard_data_sorted, library_df = load_and_prepare_data(sample_input, library_input)

        start_time = time.time()
        db_pmz, db_ms2_mz, db_ms2_int, db_labels = prepare_cpp_data(library_df)
        print(f"✓ Prepare library data time: {time.time() - start_time:.4f} s\n✓ 准备谱库数据耗时: {time.time() - start_time:.4f} s")

        start_time = time.time()
        peak_ids, sample_pmz_list, sample_rt_list, sample_ms2_list, sample_intensity_list = prepare_sample_data(
            standard_data_sorted
        )
        print(f"✓ Prepare sample data time: {time.time() - start_time:.4f} s\n✓ 准备样本数据耗时: {time.time() - start_time:.4f} s")

        start_time = time.time()
        search_engine = initialize_cpp_search_engine(
            db_pmz,
            db_ms2_mz,
            db_ms2_int,
            db_labels
        )
        print(f"✓ Load database into search engine time: {time.time() - start_time:.4f} s\n✓ 将数据库载入搜索引擎耗时: {time.time() - start_time:.4f} s")

        start_time = time.time()
        similarity_engine = initialize_cpp_similarity_engine(library_df)
        print(f"✓ Load database into similarity engine time: {time.time() - start_time:.4f} s\n✓ 将数据库载入相似度引擎耗时: {time.time() - start_time:.4f} s")

        matches = run_single_search(
            search_engine,
            peak_ids,
            sample_pmz_list,
            sample_rt_list,
            sample_ms2_list,
            sample_intensity_list,
            ppm_err1_input,
            da_err1_input,
            ppm_err2_input,
            da_err2_input
        )

        results_df = convert_matches_to_dataframe(matches, search_engine)

        results_df = postprocess_results_before_output(
            results_df,
            similarity_engine,
            ppm_err2=ppm_err2_input,
            da_err2=da_err2_input,
            jcr_cutoff=0.4
        )

        annotated_peak_ids = set()
        if results_df is not None and not results_df.empty and "peak_id" in results_df.columns:
            annotated_peak_ids = set(results_df["peak_id"].astype(int).tolist())

        unannotated_records = []
        for peak_id, sample_pmz, sample_rt, sample_ms2_mz, sample_ms2_int in zip(
            peak_ids,
            sample_pmz_list,
            sample_rt_list,
            sample_ms2_list,
            sample_intensity_list
        ):
            if int(peak_id) not in annotated_peak_ids:
                unannotated_records.append({
                    "peak_id": int(peak_id),
                    "sample_pmz": float(sample_pmz),
                    "sample_rt": float(sample_rt),
                    "sample_ms2_mz": "_".join(map(str, sample_ms2_mz)),
                    "sample_ms2_int": "_".join(map(str, sample_ms2_int))
                })

        unannotated_df = pd.DataFrame(
            unannotated_records,
            columns=[
                "peak_id",
                "sample_pmz",
                "sample_rt",
                "sample_ms2_mz",
                "sample_ms2_int"
            ]
        )

        unannotated_output = f"{sample_input.replace('.pkl', '')}_unannotated.csv"
        unannotated_df.to_csv(unannotated_output, index=False)
        print(f"✓ Unannotated result saved: {unannotated_output}\n✓ 未注释结果已保存: {unannotated_output}")
        print(f"✓ Unannotated peaks: {len(unannotated_df)}\n✓ 未注释峰数: {len(unannotated_df)}")

        save_results(results_df, sample_input.replace('.pkl', ''))

        print("🎉 Search complete!\n🎉 搜索完成！")
        print("=" * 60)

        return results_df

    except Exception as e:
        print(f"❌ Errors in program execution: {str(e)}\n❌ 程序执行出错: {str(e)}")
        import traceback
        traceback.print_exc()
        return None