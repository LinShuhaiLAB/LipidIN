import pandas as pd
import numpy as np
import os
import re
from functools import lru_cache
import polars as pl
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
from typing import Union, List
import warnings
import shutil
import time
import subprocess
warnings.filterwarnings('ignore')

def split_dataframe_for_parallel(df: pl.DataFrame, n_chunks: int) -> List[pl.DataFrame]:
    total_rows = len(df)
    chunk_size = total_rows // n_chunks
    remainder = total_rows % n_chunks

    chunks = []
    start_idx = 0

    for i in range(n_chunks):
        current_chunk_size = chunk_size + (1 if i < remainder else 0)
        chunks.append(df.slice(start_idx, current_chunk_size))
        start_idx += current_chunk_size

    return chunks
def process_single_chunk_polars(chunk: pl.DataFrame, chunk_id: int) -> pl.DataFrame:
    try:
        
        chunk = chunk.with_columns([
            pl.col("Label").str.extract(r"(.+?)_DTL_(.+?)_CAH_(.+?)_ADD_(.+?)_FML_(.+)", 1).alias("level"),
            pl.col("Label").str.extract(r"(.+?)_DTL_(.+?)_CAH_(.+?)_ADD_(.+?)_FML_(.+)", 2).alias("Title"),
            pl.col("Label").str.extract(r"(.+?)_DTL_(.+?)_CAH_(.+?)_ADD_(.+?)_FML_(.+)", 3).alias("annotation"),
            pl.col("Label").str.extract(r"(.+?)_DTL_(.+?)_CAH_(.+?)_ADD_(.+?)_FML_(.+)", 4).alias("Adduct"),
            pl.col("Label").str.extract(r"(.+?)_DTL_(.+?)_CAH_(.+?)_ADD_(.+?)_FML_(.+)", 5).alias("Formula"),
        ])

        
        grouped = (chunk
                   .group_by(['PeakID', 'annotation', 'level'])
                   .agg(pl.col('Matched_Intensity').sum())
                   )

        
        pivot_df = grouped.pivot(
            values='Matched_Intensity',
            index=['PeakID', 'annotation'],
            columns='level',
            aggregate_function='sum'
        ).fill_null(0)

        
        lvl_cols = [col for col in pivot_df.columns if 'LVL' in col]
        if lvl_cols:
            pivot_df = pivot_df.with_columns(
                pl.sum_horizontal(lvl_cols).alias('Total_Intensity')
            )
        else:
            pivot_df = pivot_df.with_columns(pl.lit(0.0).alias('Total_Intensity'))

        
        meta_df = (chunk
        .group_by(['PeakID', 'annotation'])
        .agg([
            pl.col('Title').first(),
            pl.col('Adduct').first(),
            pl.col('Formula').first(),
            pl.col('Sample_RT').first(),
            pl.col('Sample_PMZ').first()
        ])
        )

        
        result = meta_df.join(pivot_df, on=['PeakID', 'annotation'], how='inner')

        return result

    except Exception as e:
        print(f"Chunk {chunk_id} processing error: {e}\n分块 {chunk_id} 处理出错: {e}")
        import traceback
        traceback.print_exc()
        return None
def merge_chunks_polars(chunks: List[pl.DataFrame]) -> pl.DataFrame:
    valid_chunks = [c for c in chunks if c is not None]

    if not valid_chunks:
        raise ValueError("All chunks processing failed")

    
    combined = pl.concat(valid_chunks, how='diagonal')

    lvl_cols = [col for col in combined.columns if 'LVL' in col]

    
    agg_exprs = [
        pl.col('Title').first(),
        pl.col('Adduct').first(),
        pl.col('Formula').first(),
        pl.col('Sample_RT').first(),
        pl.col('Sample_PMZ').first(),
    ]

    
    for col in lvl_cols:
        agg_exprs.append(pl.col(col).sum())

    
    final_result = (combined
                    .group_by(['PeakID', 'annotation'])
                    .agg(agg_exprs)
                    )

    
    if lvl_cols:
        final_result = final_result.with_columns(
            pl.sum_horizontal(lvl_cols).alias('Total_Intensity')
        )
    else:
        final_result = final_result.with_columns(pl.lit(0.0).alias('Total_Intensity'))

    return final_result
def process_spectral_data_parallel_polars(input_data: Union[str, pd.DataFrame],
                                          n_workers: int = None) -> pl.DataFrame:


    if n_workers is None:
        n_workers = max(1, mp.cpu_count() - 1)

    
    if isinstance(input_data, str):
        df = pl.read_csv(input_data, low_memory=False)
    elif isinstance(input_data, pd.DataFrame):
        df = pl.from_pandas(input_data)
    else:
        df = input_data

    chunks = split_dataframe_for_parallel(df, n_workers)

    results = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_single_chunk_polars, chunk, i): i
                   for i, chunk in enumerate(chunks)}

        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
            except Exception as e:
                print(f"Chunk {chunk_id} failed: {e}\n分块 {chunk_id} 失败: {e}")

    final_result = merge_chunks_polars(results)
    return final_result
@lru_cache(maxsize=1)
def load_dynamic_dict_polars(dynamic_file_path: str) -> dict:

    try:
        da_D = pl.read_csv(dynamic_file_path)
        keys = list(zip(da_D['c2'].to_list(), da_D['c1'].to_list()))
        values = da_D.select(pl.col("*").exclude(['c1', 'c2'])).to_numpy().flatten()
        return dict(zip(keys, values))
    except Exception as e:
        print(f"Failed to load dynamic dictionary: {e}\n加载动态字典失败: {e}")
        return {}
def calculate_jcs_polars(df: pl.DataFrame,
                         dynamic_file: str,
                         work_dir: str = "") -> pl.DataFrame:

    dynamic_file_path = os.path.join(work_dir, dynamic_file) if work_dir else dynamic_file

    if not os.path.exists(dynamic_file_path):
        print(f"Dynamic level file not found: {dynamic_file_path}\n未找到动态等级文件: {dynamic_file_path}")
        return df.with_columns(pl.lit(None).alias('Total_JCS'))

    dynamic_dict = load_dynamic_dict_polars(dynamic_file_path)

    if not dynamic_dict:
        return df.with_columns(pl.lit(None).alias('Total_JCS'))

    lvl_cols = [col for col in df.columns if 'LVL' in col]

    if not lvl_cols:
        return df.with_columns(pl.lit(None).alias('Total_JCS'))

    
    df = df.with_columns([
        pl.col('Title').str.split(' ').list.first().alias('subclass')
    ])

    
    df = df.with_columns([
        pl.sum_horizontal([
            (pl.col(col) != 0).cast(pl.Int32)
            for col in lvl_cols
        ]).alias('non_zero_count')
    ])

    
    try:
        dynamic_df = pl.DataFrame({
            'Adduct': [k[0] for k in dynamic_dict.keys()],
            'subclass': [k[1] for k in dynamic_dict.keys()],
            'dynamic_value': list(dynamic_dict.values())
        })

        df = df.join(dynamic_df, on=['Adduct', 'subclass'], how='left')
    except Exception:
        
        pandas_df = df.select(['Adduct', 'subclass']).to_pandas()
        dynamic_values = np.array([
            dynamic_dict.get((row['Adduct'], row['subclass']), np.nan)
            for _, row in pandas_df.iterrows()
        ])
        df = df.with_columns(pl.Series('dynamic_value', dynamic_values))

    
    df = df.with_columns([
        pl.when((pl.col('dynamic_value') > 0) & pl.col('dynamic_value').is_not_null())
        .then(pl.col('non_zero_count').cast(pl.Float64) / pl.col('dynamic_value'))
        .otherwise(None)
        .alias('Total_JCS')
    ])

    df = df.drop(['subclass', 'non_zero_count', 'dynamic_value'])

    return df
def filter_by_jcs_polars(df: pl.DataFrame,
                         jcs_threshold: float = 0.5,
                         top_n: int = 10) -> pl.DataFrame:

    
    high_jcs = df.filter(pl.col('Total_JCS') >= jcs_threshold)

    top_n_records = (df
                     .sort(['PeakID', 'Total_JCS'], descending=[False, True])
                     .group_by('PeakID', maintain_order=True)
                     .head(top_n)
                     )

    
    combined = pl.concat([high_jcs, top_n_records], how='diagonal').unique(
        subset=['PeakID', 'annotation'], maintain_order=False
    )

    return combined
def step1_process_EQ2(
        input_file: Union[str, pd.DataFrame],
        dynamic_level_file: str = None,
        jcs_threshold: float = 0.5,
        top_n: int = 10,
        n_workers: int = None,
        all_Dynamic_level=r"D:\bio_inf\02 LipidIN 3.0\pycode\all_Dynamic_level.csv"
) -> pd.DataFrame:

    if dynamic_level_file is None:
        dynamic_level_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            all_Dynamic_level
        )

    try:
        result_df = process_spectral_data_parallel_polars(input_file, n_workers)

        if result_df is None or len(result_df) == 0:
            print("Spectral processing failed\n谱图处理失败")
            return None

        result_df = calculate_jcs_polars(result_df, dynamic_level_file, "")
        filtered_df = filter_by_jcs_polars(result_df, jcs_threshold, top_n)
        filtered_df = filtered_df.filter(pl.col("Total_JCS") >= 0.33)

        
        filtered_df = filtered_df.with_columns([
            (pl.col('Total_Intensity') * pl.col('Total_JCS').fill_null(0)).alias('Combined_Score')
        ])
        filtered_df = (filtered_df
                    .sort(['PeakID', 'Combined_Score'], descending=[False, True])
                    .group_by('PeakID', maintain_order=True)
                    .head(top_n)
                    )

        final_df = (filtered_df
                    .sort(['annotation', 'Combined_Score'], descending=[False, True])
                    .group_by('annotation', maintain_order=True)
                    .head(top_n)
                    .drop('Combined_Score')
                    )

        return final_df.to_pandas()

    except Exception as e:
        print(f"Processing error: {e}\n处理出错: {e}")
        import traceback
        traceback.print_exc()
        return None