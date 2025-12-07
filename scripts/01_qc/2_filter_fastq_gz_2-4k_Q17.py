#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
filter_fastq_gz_2-4k_Q17.py
---------------------------
递归扫描 /Volumes/OWC/Analysis/trimmed_out 下所有 .fastq.gz 文件，
按如下条件过滤并输出到 /Volumes/OWC/Analysis/filter_out_2-4k_Q17：
- read 长度在 [2000, 4000] 区间
- 平均 Q score > 17

为最大化利用本机性能，采用多进程并行（按文件级别并行），
进程数自动等于 min(可用CPU核心数, 待处理文件数)。

用法：
    python3 /Volumes/OWC/Analysis/filter_fastq_gz_2-4k_Q17.py
"""

import os
import gzip
import multiprocessing as mp
from typing import List, Tuple
import argparse


# 配置路径
INPUT_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_dorado"
OUTPUT_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_filter_2"

# 过滤阈值
MIN_LENGTH = 2000
MAX_LENGTH = 4000
MIN_QSCORE = 17.0


def filter_single_file(in_fastq_gz_path: str,
                       out_fastq_gz_path: str,
                       min_len: int = MIN_LENGTH,
                       max_len: int = MAX_LENGTH,
                       min_q: float = MIN_QSCORE) -> None:
    """对单个 .fastq.gz 文件进行过滤，并写入对应的输出 .fastq.gz 文件。

    具备容错：若文件不是有效 gzip/损坏/解码失败，打印警告并跳过。
    """
    try:
        out_dir = os.path.dirname(out_fastq_gz_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        total_reads = 0
        passed_reads = 0

        # 为提高速度，将 gzip 压缩级别设为 1（更快，文件稍大）
        with gzip.open(in_fastq_gz_path, mode="rt", encoding="ascii", errors="strict") as fin, \
             gzip.open(out_fastq_gz_path, mode="wt", encoding="ascii", compresslevel=1) as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                sequence = fin.readline()
                plus_line = fin.readline()
                quality = fin.readline()

                # 防止非完整 read（意外截断）
                if not quality:
                    break

                header = header.rstrip("\n")
                sequence = sequence.rstrip("\n")
                plus_line = plus_line.rstrip("\n")
                quality = quality.rstrip("\n")

                total_reads += 1

                read_length = len(sequence)
                if read_length < min_len or read_length > max_len:
                    continue

                # 质量行长度应与序列长度一致
                if len(quality) != read_length:
                    continue

                # 计算平均 Q 值（Phred+33）
                total_qscore = 0
                for qual_char in quality:
                    total_qscore += (ord(qual_char) - 33)
                average_qscore = total_qscore / float(read_length) if read_length else 0.0

                if average_qscore > min_q:
                    fout.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
                    passed_reads += 1

        print(f"[{os.path.basename(in_fastq_gz_path)}] Total reads: {total_reads}, Passed: {passed_reads}, Output: {os.path.basename(out_fastq_gz_path)}")
    except Exception as exc:
        print(f"[SKIP] {in_fastq_gz_path} -> 跳过原因: {type(exc).__name__}: {exc}")


def _is_probably_gzip(file_path: str) -> bool:
    """快速判断文件是否为 gzip（检查魔数 0x1f 0x8b）。"""
    try:
        with open(file_path, "rb") as f:
            magic = f.read(2)
        return magic == b"\x1f\x8b"
    except Exception:
        return False


def _is_valid_gzip(file_path: str) -> bool:
    """尝试打开并读取少量数据以确认 gzip 文件有效。"""
    if not _is_probably_gzip(file_path):
        return False
    try:
        with gzip.open(file_path, "rb") as f:
            _ = f.read(1)
        return True
    except Exception:
        return False


def collect_file_pairs(root_dir: str, output_root_dir: str) -> List[Tuple[str, str]]:
    """递归收集 (输入文件, 输出文件) 对，并保持子目录结构。

    - 跳过以 `._` 开头的伪文件
    - 跳过非 gzip/空文件
    """
    pairs: List[Tuple[str, str]] = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if not filename.endswith(".fastq.gz"):
                continue
            if filename.startswith("._"):
                print(f"[SKIP] {os.path.join(dirpath, filename)} -> 伪文件（前缀 ._）")
                continue
            in_path = os.path.join(dirpath, filename)
            # 跳过非 gzip 或空文件
            try:
                if os.path.getsize(in_path) == 0:
                    print(f"[SKIP] {in_path} -> 空文件")
                    continue
            except Exception:
                print(f"[SKIP] {in_path} -> 无法获取文件大小")
                continue
            if not _is_probably_gzip(in_path):
                print(f"[SKIP] {in_path} -> 非 gzip 文件")
                continue

            rel_dir = os.path.relpath(dirpath, root_dir)
            target_dir = os.path.join(output_root_dir, rel_dir) if rel_dir != "." else output_root_dir
            out_filename = filename[:-9] + "_filtered.fastq.gz"  # 去掉 ".fastq.gz" 共 9 个字符
            out_path = os.path.join(target_dir, out_filename)
            pairs.append((in_path, out_path))
    return pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter Nanopore fastq.gz reads by length and Q score")
    parser.add_argument("--only-missing", action="store_true",
                        help="仅处理输出缺失或损坏（无效 gzip）的文件")
    parser.add_argument("--workers", type=int, default=0,
                        help="并行进程数（默认自动=CPU核数与任务数的较小值）")
    return parser.parse_args()


def main():
    args = parse_args()
    file_pairs = collect_file_pairs(INPUT_DIR, OUTPUT_DIR)
    if not file_pairs:
        print(f"在目录 {INPUT_DIR} 下未找到任何 .fastq.gz 文件。请检查。")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 仅重跑缺失或损坏的输出
    if args.only_missing:
        rerun_pairs: List[Tuple[str, str]] = []
        for in_path, out_path in file_pairs:
            # 输出不存在或不是有效 gzip 即需要重跑
            if (not os.path.exists(out_path)) or (not _is_valid_gzip(out_path)):
                rerun_pairs.append((in_path, out_path))
        file_pairs = rerun_pairs
        if not file_pairs:
            print("没有需要重跑的文件，已全部完成或输出有效。")
            return

    available_cpus = mp.cpu_count()
    requested_workers = args.workers if args.workers and args.workers > 0 else None
    num_workers = min(available_cpus, len(file_pairs)) if requested_workers is None else min(requested_workers, len(file_pairs))
    print(f"检测到 CPU 核心数: {available_cpus}，将并行处理 {num_workers} 个文件。")

    with mp.Pool(processes=num_workers) as pool:
        results = []
        for in_path, out_path in file_pairs:
            r = pool.apply_async(filter_single_file, args=(in_path, out_path))
            results.append(r)

        # 等待所有子任务完成
        for r in results:
            r.get()

    print("过滤完成。所有处理结果已输出到:", OUTPUT_DIR)


if __name__ == "__main__":
    main()


