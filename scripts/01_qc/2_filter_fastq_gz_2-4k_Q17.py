#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter FASTQ reads by length and quality score.

Filters reads by:
- Length: [2000, 4000] bp (HBV genome ~3.2kb)
- Mean Q score: > 17

Usage:
    python 2_filter_fastq_gz_2-4k_Q17.py --input /path/to/input --output /path/to/output
"""

import os
import gzip
import multiprocessing as mp
from typing import List, Tuple
import argparse

# Default paths (relative to project)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))

INPUT_DIR = os.environ.get("HBV_FILTER_INPUT", os.path.join(_PROJECT_DIR, "fastq_dorado"))
OUTPUT_DIR = os.environ.get("HBV_FILTER_OUTPUT", os.path.join(_PROJECT_DIR, "fastq_filter_2"))

# Filter thresholds
MIN_LENGTH = 2000
MAX_LENGTH = 4000
MIN_QSCORE = 17.0


def filter_single_file(in_fastq_gz_path: str,
                       out_fastq_gz_path: str,
                       min_len: int = MIN_LENGTH,
                       max_len: int = MAX_LENGTH,
                       min_q: float = MIN_QSCORE) -> None:
    """Filter a single FASTQ.gz file by length and quality."""
    try:
        out_dir = os.path.dirname(out_fastq_gz_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        total_reads = 0
        passed_reads = 0

        with gzip.open(in_fastq_gz_path, mode="rt", encoding="ascii", errors="strict") as fin, \
             gzip.open(out_fastq_gz_path, mode="wt", encoding="ascii", compresslevel=1) as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                sequence = fin.readline()
                plus_line = fin.readline()
                quality = fin.readline()

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

                if len(quality) != read_length:
                    continue

                # Calculate mean Q score (Phred+33)
                total_qscore = sum(ord(c) - 33 for c in quality)
                average_qscore = total_qscore / float(read_length) if read_length else 0.0

                if average_qscore > min_q:
                    fout.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
                    passed_reads += 1

        print(f"[{os.path.basename(in_fastq_gz_path)}] Total: {total_reads}, Passed: {passed_reads}")
    except Exception as exc:
        print(f"[SKIP] {in_fastq_gz_path} -> {type(exc).__name__}: {exc}")


def _is_probably_gzip(file_path: str) -> bool:
    """Check if file has gzip magic number."""
    try:
        with open(file_path, "rb") as f:
            magic = f.read(2)
        return magic == b"\x1f\x8b"
    except Exception:
        return False


def _is_valid_gzip(file_path: str) -> bool:
    """Verify gzip file is readable."""
    if not _is_probably_gzip(file_path):
        return False
    try:
        with gzip.open(file_path, "rb") as f:
            _ = f.read(1)
        return True
    except Exception:
        return False


def collect_file_pairs(root_dir: str, output_root_dir: str) -> List[Tuple[str, str]]:
    """Collect (input, output) file pairs, preserving directory structure."""
    pairs: List[Tuple[str, str]] = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if not filename.endswith(".fastq.gz"):
                continue
            if filename.startswith("._"):
                print(f"[SKIP] {os.path.join(dirpath, filename)} -> pseudo file")
                continue
            in_path = os.path.join(dirpath, filename)
            try:
                if os.path.getsize(in_path) == 0:
                    print(f"[SKIP] {in_path} -> empty file")
                    continue
            except Exception:
                print(f"[SKIP] {in_path} -> cannot get file size")
                continue
            if not _is_probably_gzip(in_path):
                print(f"[SKIP] {in_path} -> not gzip")
                continue

            rel_dir = os.path.relpath(dirpath, root_dir)
            target_dir = os.path.join(output_root_dir, rel_dir) if rel_dir != "." else output_root_dir
            out_filename = filename[:-9] + "_filtered.fastq.gz"
            out_path = os.path.join(target_dir, out_filename)
            pairs.append((in_path, out_path))
    return pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter Nanopore reads by length and Q score")
    parser.add_argument("--input", "-i", default=INPUT_DIR, help="Input directory")
    parser.add_argument("--output", "-o", default=OUTPUT_DIR, help="Output directory")
    parser.add_argument("--min-len", type=int, default=MIN_LENGTH, help="Minimum read length")
    parser.add_argument("--max-len", type=int, default=MAX_LENGTH, help="Maximum read length")
    parser.add_argument("--min-quality", type=float, default=MIN_QSCORE, help="Minimum mean Q score")
    parser.add_argument("--only-missing", action="store_true", help="Only process missing/invalid outputs")
    parser.add_argument("--workers", type=int, default=0, help="Number of parallel workers (0=auto)")
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = args.input
    output_dir = args.output
    
    file_pairs = collect_file_pairs(input_dir, output_dir)
    if not file_pairs:
        print(f"No .fastq.gz files found in {input_dir}")
        return

    os.makedirs(output_dir, exist_ok=True)

    if args.only_missing:
        rerun_pairs: List[Tuple[str, str]] = []
        for in_path, out_path in file_pairs:
            if (not os.path.exists(out_path)) or (not _is_valid_gzip(out_path)):
                rerun_pairs.append((in_path, out_path))
        file_pairs = rerun_pairs
        if not file_pairs:
            print("All outputs already exist and are valid.")
            return

    available_cpus = mp.cpu_count()
    requested_workers = args.workers if args.workers > 0 else None
    num_workers = min(available_cpus, len(file_pairs)) if requested_workers is None else min(requested_workers, len(file_pairs))
    print(f"CPU cores: {available_cpus}, using {num_workers} workers")

    with mp.Pool(processes=num_workers) as pool:
        results = []
        for in_path, out_path in file_pairs:
            r = pool.apply_async(filter_single_file, args=(in_path, out_path, args.min_len, args.max_len, args.min_quality))
            results.append(r)
        for r in results:
            r.get()

    print(f"Done. Output: {output_dir}")


if __name__ == "__main__":
    main()
