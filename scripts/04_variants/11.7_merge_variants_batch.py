#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge Accepted and Candidate Variants Script
=============================================

Purpose:
    Merge accepted.tsv and candidates.tsv files from the transformed coordinate results.
    Combines both files into a single unified output with standardized columns and
    source tracking.

Required Environment:
    - pandas (data manipulation)
    
    Install via conda/mamba:
        conda activate ivar  # or any environment with pandas

Input:
    - {WORK_DIR}/variants/transformed/{sample_id}/
      * accepted.tsv - High-confidence variants (transformed coordinates)
      * candidates.tsv - Candidate variants (transformed coordinates)

Output:
    - {WORK_DIR}/variants/combined/{sample_id}_combined.tsv
      * Merged file with all variants
      * Includes file_source column to track origin
    - {WORK_DIR}/variants/combined/
      * merge_summary.csv - Per-sample merge statistics
      * merge_statistics.txt - Overall statistics

Adjustable Parameters:
    --test SAMPLE: Test mode for single sample
    --list: List available samples
    --jobs N: Number of parallel jobs (default: 8)

Usage:
    # Test single sample
    python 11.7_merge_variants_batch.py --test 10090
    
    # List available samples
    python 11.7_merge_variants_batch.py --list
    
    # Process all samples (with confirmation)
    python 11.7_merge_variants_batch.py --jobs 8

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import sys
import csv
import argparse
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# ========== PATH CONFIGURATION ==========
# Working directory
WORK_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/variants_call_10"

# Paths
INPUT_DIR = os.path.join(WORK_DIR, "variants", "transformed")
OUTPUT_DIR = os.path.join(WORK_DIR, "variants", "combined")
# ========================================

# Target files
TARGET_FILES = ['accepted.tsv', 'candidates.tsv']

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_sample_directories():
    """Get all sample directories with target files"""
    input_path = Path(INPUT_DIR)
    sample_dirs = {}
    
    if not input_path.exists():
        print(f"[ERROR] Input directory not found: {INPUT_DIR}")
        return sample_dirs
    
    for sample_dir in input_path.iterdir():
        if not sample_dir.is_dir():
            continue
        
        sample_files = {}
        files_found = 0
        
        for target_file in TARGET_FILES:
            file_path = sample_dir / target_file
            if file_path.exists():
                sample_files[target_file] = file_path
                files_found += 1
        
        if files_found > 0:
            sample_dirs[sample_dir.name] = {
                'input_dir': sample_dir,
                'files': sample_files
            }
    
    return sample_dirs

def standardize_columns(df, file_type):
    """Standardize column names with accepted.tsv as baseline"""
    # Standard column structure
    standard_columns = [
        'sample_id', 'chrom', 'pos', 'ref', 'alt', 'var_type', 'struct_len',
        'AF_primary', 'DP_primary', 'AF_clair3', 'DP_clair3', 'AD_clair3',
        'AF_ivar', 'DP_ivar', 'AD_ivar', 'source_label', 'tier', 'decision', 'reason',
        'file_source'
    ]
    
    # Add file_source column
    df['file_source'] = file_type
    
    if file_type == 'accepted':
        # accepted.tsv should have correct structure, just rename if needed
        if 'AF_primary' not in df.columns and 'AF' in df.columns:
            df = df.rename(columns={'AF': 'AF_primary', 'DP': 'DP_primary'})
    
    elif file_type == 'candidates':
        # candidates.tsv needs column mapping
        column_mapping = {}
        if 'AF' in df.columns:
            column_mapping['AF'] = 'AF_primary'
        if 'DP' in df.columns:
            column_mapping['DP'] = 'DP_primary'
        if 'AD' in df.columns:
            column_mapping['AD'] = 'AD_primary'
        
        df = df.rename(columns=column_mapping)
    
    # Ensure all standard columns exist
    for col in standard_columns:
        if col not in df.columns:
            df[col] = ''
    
    # Reorder columns
    df = df[standard_columns]
    return df

def merge_sample_files(sample_id, sample_info):
    """Merge accepted.tsv and candidates.tsv for a single sample"""
    print(f"[Merging] {sample_id}")
    
    dataframes = []
    merge_info = {
        'sample_id': sample_id,
        'files_processed': {},
        'total_variants': 0,
        'processing_time': datetime.now().isoformat()
    }
    
    try:
        for file_type in ['accepted', 'candidates']:
            file_key = f"{file_type}.tsv"
            if file_key in sample_info['files']:
                file_path = sample_info['files'][file_key]
                
                df = pd.read_csv(file_path, sep='\t', dtype=str, na_filter=False)
                df_standardized = standardize_columns(df, file_type)
                dataframes.append(df_standardized)
                
                merge_info['files_processed'][file_type] = {
                    'input_file': str(file_path),
                    'variants_count': len(df_standardized)
                }
                merge_info['total_variants'] += len(df_standardized)
        
        if not dataframes:
            print(f"[WARN] {sample_id}: No files to merge")
            return None
        
        merged_df = pd.concat(dataframes, ignore_index=True)
        output_file = Path(OUTPUT_DIR) / f"{sample_id}_combined.tsv"
        merged_df.to_csv(output_file, sep='\t', index=False, na_rep='')
        
        merge_info['output_file'] = str(output_file)
        merge_info['merged_variants'] = len(merged_df)
        
        print(f"[OK] {sample_id}: {merge_info['total_variants']} variants → {merge_info['merged_variants']} rows")
        return merge_info
        
    except Exception as e:
        print(f"[ERROR] {sample_id}: {e}")
        return None

def generate_merge_summary(results):
    """Generate merge summary report"""
    summary_file = Path(OUTPUT_DIR) / "merge_summary.csv"
    
    try:
        with open(summary_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['# Merge Summary'])
            writer.writerow(['Merge time', datetime.now().isoformat()])
            writer.writerow([])
            writer.writerow(['Sample ID', 'accepted variants', 'candidates variants', 'Total', 'Output file'])
            
            for result in results:
                if result:
                    accepted_count = result['files_processed'].get('accepted', {}).get('variants_count', 0)
                    candidates_count = result['files_processed'].get('candidates', {}).get('variants_count', 0)
                    writer.writerow([
                        result['sample_id'],
                        accepted_count,
                        candidates_count,
                        result['total_variants'],
                        result['output_file']
                    ])
        
        print(f"[INFO] Summary saved to: {summary_file}")
        
        # Statistics file
        successful = [r for r in results if r is not None]
        total_accepted = sum(r['files_processed'].get('accepted', {}).get('variants_count', 0) for r in successful)
        total_candidates = sum(r['files_processed'].get('candidates', {}).get('variants_count', 0) for r in successful)
        total_merged = sum(r['merged_variants'] for r in successful)
        
        stats_file = Path(OUTPUT_DIR) / "merge_statistics.txt"
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("Merge Statistics\n")
            f.write("=" * 40 + "\n")
            f.write(f"Processed samples: {len(successful)}\n")
            f.write(f"Total accepted: {total_accepted:,}\n")
            f.write(f"Total candidates: {total_candidates:,}\n")
            f.write(f"Total merged rows: {total_merged:,}\n")
            f.write(f"Average per sample: {total_merged / len(successful):.1f}\n" if successful else "Average: 0\n")
        
        print(f"[INFO] Statistics saved to: {stats_file}")
        
    except Exception as e:
        print(f"[ERROR] Generating summary: {e}")

def main():
    parser = argparse.ArgumentParser(description='Merge accepted.tsv and candidates.tsv files')
    parser.add_argument('--test', type=str, help='Test single sample (sample ID)')
    parser.add_argument('--jobs', type=int, default=8, help='Parallel jobs (default: 8)')
    parser.add_argument('--list', action='store_true', help='List all available samples')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Variant Merge Script v2.0 for variants_call_10")
    print(f"Input: {INPUT_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 60)
    
    sample_dirs = get_sample_directories()
    
    if args.list:
        print(f"Found {len(sample_dirs)} samples:")
        for i, (sample_id, info) in enumerate(sorted(sample_dirs.items()), 1):
            files_str = ', '.join(info['files'].keys())
            print(f"{i:3d}. {sample_id} ({files_str})")
        return
    
    if not sample_dirs:
        print("[ERROR] No samples found")
        return
    
    if args.test:
        if args.test not in sample_dirs:
            print(f"[ERROR] Sample {args.test} not found")
            return
        
        result = merge_sample_files(args.test, sample_dirs[args.test])
        if result:
            print(f"\nSample {args.test} merge complete:")
            print(f"  Total variants: {result['total_variants']}")
            print(f"  Merged rows: {result['merged_variants']}")
            print(f"  Output: {result['output_file']}")
    else:
        # Batch processing
        print(f"\nMerging {len(sample_dirs)} samples (parallel: {args.jobs})")
        
        results = []
        sample_items = list(sample_dirs.items())
        
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {
                executor.submit(merge_sample_files, sample_id, sample_info): sample_id
                for sample_id, sample_info in sample_items
            }
            
            completed = 0
            for future in as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1
                    print(f"Progress: {completed}/{len(sample_items)}")
                except Exception as e:
                    print(f"[ERROR] {sample_id}: {e}")
        
        generate_merge_summary(results)
        
        successful = [r for r in results if r is not None]
        total_input = sum(r['total_variants'] for r in successful)
        total_output = sum(r['merged_variants'] for r in successful)
        
        print(f"\nMerge complete!")
        print(f"  Successful: {len(successful)}/{len(sample_items)} samples")
        print(f"  Input variants: {total_input:,}")
        print(f"  Merged rows: {total_output:,}")
        print(f"  Output: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
