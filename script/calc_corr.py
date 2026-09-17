#!/share/data1/software/miniconda3/bin/python3
# -*- encoding:utf-8 -*-

from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from multiprocessing import Manager
import threading
import time
import logging
from scipy import stats
import numpy as np
import pandas as pd
import argparse
import sqlite3
import os
import asyncio
import math

__doc__ = '''
对于输入matrix格式：
    Name    col1    col2    col3
    row1    0.1     0.3     0.45
    row2    0.38    0.9     0.38
    row3    0.12    0.01    0.47
    row4    0.43    0.03    0.39

注意：第一行，小心错位
'''


def get_args():

    corr_methods = ['spearman','pearson','fisher','Chi-Square','Mutual Information','Correlation Ratio','Maximal Information Coefficient','euclidean', 'Cosine Similarity','Kendalls','Point-Biserial']

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="matrix_1", required=True, type=str, help="Input data1")
    parser.add_argument("-I", metavar="matrix_2", required=False, type=str, help="Input data2")
    parser.add_argument("-o", metavar="out", required=True, type=str, help="Output")
    parser.add_argument("-t1", metavar="title", default='T', choices=['T','F'], type=str, help="data1是否有列名(F/T:default)")
    parser.add_argument("-t2", metavar="title", default='T', choices=['T','F'], type=str, help="data2是否有列名(F/T:default)")
    parser.add_argument("-M",metavar="Method",required=False, default="spearman",
                        choices = corr_methods,
                        type=str, help=f"!!!!,暂时只能用spearman,pearson,fisher. 相关性算法. {corr_methods}")
    parser.add_argument("-s", required=False, type=int,  default=1, choices=[0,1,2,3], help="Split Str fo matrix.[default: 1]\n0 -> \\s+\n1 -> \\t\n2 -> |\n3 -> ,")
    parser.add_argument("--skip1", type=int, default=0, help="data1 忽略前多少行(default=0).")
    parser.add_argument("--skip2", type=int, default=0, help="data2 忽略前多少行(default=0).")
    parser.add_argument("--threads", metavar="threads", type=int, default=10, help="Number of threads.")
    parser.add_argument("--axis", metavar="", default="row", type=str, help="by col/row:default.")
    args = parser.parse_args()
    return args

### 异步写入结果到数据库
async def async_write_to_db(df, db_name):
    loop = asyncio.get_running_loop()
    await loop.run_in_executor(None, lambda: df.to_sql('results', sqlite3.connect(db_name), if_exists='append', index=False))


def calcute_pearson(df1, start1, end1, df2, db_name=None, progress=None, index=None):
    try:
        names1 = df1.index[start1:end1]

        nrows1 = end1-start1

        names2 = df2.index
        nrows2 = len(names2)
        st = time.time()

        ## init results data
        nrow = nrows2 * (end1 - start1)

        row_num = 0
        corrs = np.full(nrow, np.nan)
        pvals = np.full(nrow, np.nan)
        for i in range(start1, end1):
            x = df1.iloc[i,:]
            for j in range(nrows2):
                y = df2.iloc[j,:]
                nas = np.logical_or(x.isna(), y.isna())
                corrs[row_num],pvals[row_num] = stats.pearsonr(x[~nas], y[~nas])
                row_num  += 1
            progress[index] += 1

        res = pd.DataFrame({
            'A':  np.repeat(names1, nrows2), 
            'B':  np.tile(names2, nrows1),
            'corr': corrs,
            'pvalue': pvals
            })

    except Exception as e:
        logging.error(f"Error in spearman calculation: {e}")
    
    # 异步执行写入文件
    asyncio.run(async_write_to_db(res, db_name))
    logging.info(f"# {start1} - {end1}\t{time.time() - st} Sec.")


def calcute_fisher(df1, start1, end1, df2, db_name=None, progress=None, index=None):
     
    names1 = df1.index[start1:end1]
    nrow1 = end1-start1

    names2 = df2.index
    nrow2 = len(names2)
    st = time.time()

    ## init results data
    nrow = nrow2 * (end1 - start1)

    row_num = 0
    corrs = np.full(nrow, np.nan)
    pvals = np.full(nrow, np.nan)

    x_values = df1.values[start1:end1]
    y_values = df2.values

    for i in range(x_values.shape[0]):
        x = x_values[i, :]
        for j in range(y_values.shape[0]):
            y = y_values[j, :]
            ## 计算每个组合的 0 和 1 数量
            contingency_table = np.array([
                        [ np.sum((x == 0) & (y == 0)), np.sum((x == 0) & (y == 1))],
                        [ np.sum((x == 1) & (y == 0)), np.sum((x == 1) & (y == 1))]
                    ])

            # corrs[row_num],pvals[row_num] = stats.fisher_exact(pd.crosstab(x,y), alternative='two-sided')
            corrs[row_num],pvals[row_num] = stats.fisher_exact(contingency_table, alternative='two-sided')
            row_num  += 1
        progress[index] += 1

    res = pd.DataFrame({
        'A':  [x for x in names1 for _ in range(0,nrow2) ], 
        'B':  names2.to_list() * nrow1,
        'corr': corrs,
        'pvalue': pvals
        })
    
    # 异步执行写入文件
    asyncio.run(async_write_to_db(res, db_name))
    logging.info(f"# {start1} - {end1}\t{time.time() - st} Sec.")


def read_file(inf, skp, sstr, title, method, axis):
    inf = os.path.abspath(f"{inf}") 
    logging.info(f"read file: {inf}")
    if title == "F":
        title=None
    else:
        title=0
    df = pd.read_csv(inf, skiprows=skp, header=title, sep=sstr, index_col=0)
    
    ## 判断是否按行
    if axis == "col":
        df = df.T
    ## 如果是spearman，那就先计算秩序（必须先转置，再计算，不然就错了）
    if method == "spearman":
        df = df.rank(axis=1)
    ##  如果是fisher，先将矩阵转换成binary
    elif method == "fisher":
        df = (df>0) + 0
    return df


def get_input_matrix(args):
    sstr = {0:"\s+", 1:"\t", 2:"\|", 3:","}[args.s]
    logging.info("## Step1: read data")
    logging.info(f"axis: {args.axis}")
    ## 读取文件
    df1 = read_file(args.i, args.skip1, sstr, args.t1, method=args.M, axis=args.axis)
    nf,ns = df1.shape
    if not args.I:
        df2 = df1
        logging.info("-"*50)
        logging.info(f"\tN_row\tN_col")
        logging.info(f"data\t{nf}\t\t{ns}")
    else:
        df2 = read_file(args.I, args.skip2, sstr, args.t2, method=args.M, axis=args.axis)
        df2 = df2.reindex(columns=df1.columns)
        nf2,ns2 = df2.shape
        logging.info("-"*50)
        logging.info(f"\tN_row\tN_col")
        logging.info(f"data1\t{nf}\t\t{ns}")
        logging.info(f"data2\t{nf2}\t\t{ns2}" )
    logging.info("-"*50)
    return df1, df2

def merge_results(db_name):
    # 从SQLite数据库中读取所有结果
    conn = sqlite3.connect(db_name)
    final_df = pd.read_sql_query("SELECT * FROM results", conn)
    conn.close()
    return final_df


def print_progress(progress, total):
    while True:
        current_progress = sum(progress)
        print(f"Total Progress: {current_progress}/{total}", end="\r")
        if current_progress >= total:
            print("\nMerged Data...")
            break
        time.sleep(1)

def main(args):
    
    outf = args.o
    db_name = f"{outf}.tmp.db"
    with sqlite3.connect(db_name) as conn:
         conn.execute('CREATE TABLE IF NOT EXISTS results (A TEXT, B TEXT, corr REAL, pvalue REAL)')
    
    ## get data1 and data2
    df1, df2 = get_input_matrix(args)
    
    ## 因为如果是计算spearman的话，会在读取数据时就计算数据rank，所以本质还是计算pearson
    my_method = {'spearman':calcute_pearson, 'pearson':calcute_pearson, 'fisher':calcute_fisher}.get(args.M)
    
    ## calc corr
    logging.info(f"method: {args.M}")

    ## 多线程计算
    num_threads = args.threads

    ## 将data1分割成多份
    nf1 = df1.shape[0]
    step_size = math.ceil(nf1 / num_threads)
    loop_size = nf1 // step_size + (nf1 % step_size > 0)

    names1 = df1.index
    logging.info(f"actual: {num_threads}")

    with Manager() as manager:
        progress = manager.list([0]*loop_size)

        progress_thread = threading.Thread(target=print_progress, args=(progress, nf1))
        progress_thread.start()

        ## 为了减少上下文交换，每次循环区间
        with ProcessPoolExecutor(max_workers = num_threads ) as exector:
            futures = []
            for i in range(loop_size):
                start = i * step_size
                end = min( (i + 1) * step_size, nf1)
                future = exector.submit(my_method, df1, start, end, df2, db_name, progress, i)
                futures.append(future)

            logging.info("Start subprocess...")
            for completed_future in as_completed(futures):
                try:
                    completed_future.result()  # 确保任务完成
                except Exception as e:
                    logging.error(f"Error in future: {e}")

        progress_thread.join()  # 等待进度打印线程结束

    # 关闭数据库连接
    if os.path.exists(db_name):
        ## 合并结果
        final_result = merge_results(db_name)
        os.remove(db_name)
        final_result.to_csv(f"{outf}", sep="\t", index=False)
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = get_args()
    main(args)
    
    
