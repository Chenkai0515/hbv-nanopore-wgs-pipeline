#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coordinate Transformation Script for Variant Results
=====================================================

Purpose:
    Transform variant coordinates from rotated reference (starting from DR1) back to
    original reference coordinates. Processes accepted.tsv and candidates.tsv files
    from the summarization step.

Required Environment:
    - pandas (data manipulation)
    
    Install via conda/mamba:
        conda activate ivar  # or any environment with pandas

Input:
    - {WORK_DIR}/variants/filtered/{sample_id}/
      * accepted.tsv - High-confidence variants
      * candidates.tsv - Candidate variants

Output:
    - {WORK_DIR}/variants/transformed/{sample_id}/
      * accepted.tsv - Coordinates transformed to original reference
      * candidates.tsv - Coordinates transformed to original reference
    - {WORK_DIR}/variants/transformed/
      * transform_summary.csv - Transformation summary
      * transform_detailed.csv - Detailed transformation log

Adjustable Parameters:
    --test SAMPLE: Test mode for single sample
    --list: List available samples
    --jobs N: Number of parallel jobs (default: 4)

Coordinate Transformation:
    - HBV genome length: 3221 bp
    - Unified offset: 1823
    - Formula: original_pos = ((current_pos - 1) + 1823) % 3221 + 1

Usage:
    # Test single sample
    python 11.6_coordinate_transform_batch.py --test 10090
    
    # List available samples
    python 11.6_coordinate_transform_batch.py --list
    
    # Process all samples
    python 11.6_coordinate_transform_batch.py --jobs 8

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import sys
import csv
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# ========== PATH CONFIGURATION ==========
# Working directory (can be set via HBV_WORK_DIR environment variable)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
WORK_DIR = os.environ.get("HBV_WORK_DIR", os.path.join(PROJECT_DIR, "variants_call_10"))

# Paths
INPUT_DIR = os.path.join(WORK_DIR, "variants", "filtered")
OUTPUT_DIR = os.path.join(WORK_DIR, "variants", "transformed")
# ========================================

# HBV genome parameters
GENOME_LENGTH = 3221
DEFAULT_OFFSET = 1823

# Target files
TARGET_FILES = ['accepted.tsv', 'candidates.tsv']

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def transform_coordinate(current_pos, offset=DEFAULT_OFFSET, genome_length=GENOME_LENGTH):
    """Transform coordinate from rotated to original reference"""
    original_pos = ((current_pos - 1) + offset) % genome_length + 1
    return original_pos

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
        if sample_dir.name == "cohort_summary":
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

def transform_single_file(input_file, output_file, sample_id, file_type):
    """Transform coordinates in a single TSV file"""
    variants_transformed = 0
    
    try:
        with open(input_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8', newline='') as outfile:
            
            reader = csv.DictReader(infile, delimiter='\t')
            fieldnames = reader.fieldnames
            
            if 'pos' not in fieldnames:
                print(f"[WARN] {sample_id} {file_type}: No 'pos' column found")
                return 0
            
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for row in reader:
                try:
                    current_pos = int(row['pos'])
                    original_pos = transform_coordinate(current_pos)
                    row['pos'] = str(original_pos)
                    writer.writerow(row)
                    variants_transformed += 1
                except (ValueError, KeyError) as e:
                    print(f"[WARN] {sample_id} {file_type}: Cannot process pos '{row.get('pos', 'N/A')}': {e}")
                    continue
        
        return variants_transformed
        
    except Exception as e:
        print(f"[ERROR] Transforming {input_file}: {e}")
        return 0

def process_single_sample(sample_id, sample_info):
    """Process coordinate transformation for a single sample"""
    print(f"[Processing] {sample_id}")
    
    sample_output_dir = Path(OUTPUT_DIR) / sample_id
    sample_output_dir.mkdir(exist_ok=True)
    
    results = {
        'sample_id': sample_id,
        'files_processed': {},
        'total_variants': 0,
        'offset_used': DEFAULT_OFFSET,
        'processing_time': datetime.now().isoformat()
    }
    
    for file_type, input_file in sample_info['files'].items():
        output_file = sample_output_dir / file_type
        variants_count = transform_single_file(input_file, output_file, sample_id, file_type)
        
        results['files_processed'][file_type] = {
            'input_file': str(input_file),
            'output_file': str(output_file),
            'variants_count': variants_count
        }
        results['total_variants'] += variants_count
    
    print(f"[OK] {sample_id}: Transformed {results['total_variants']} variants")
    return results

def generate_summary(results):
    """Generate transformation summary report"""
    summary_file = Path(OUTPUT_DIR) / "transform_summary.csv"
    
    try:
        with open(summary_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['# Coordinate Transformation Summary'])
            writer.writerow(['Transformation time', datetime.now().isoformat()])
            writer.writerow(['Offset used', DEFAULT_OFFSET])
            writer.writerow(['Formula', f'original_pos = ((current_pos - 1) + {DEFAULT_OFFSET}) % {GENOME_LENGTH} + 1'])
            writer.writerow([])
            writer.writerow(['Sample ID', 'accepted.tsv variants', 'candidates.tsv variants', 'Total variants'])
            
            for result in results:
                if result:
                    accepted_count = result['files_processed'].get('accepted.tsv', {}).get('variants_count', 0)
                    candidates_count = result['files_processed'].get('candidates.tsv', {}).get('variants_count', 0)
                    writer.writerow([result['sample_id'], accepted_count, candidates_count, result['total_variants']])
        
        print(f"[INFO] Summary saved to: {summary_file}")
        
    except Exception as e:
        print(f"[ERROR] Generating summary: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='Transform variant coordinates from rotated to original reference'
    )
    parser.add_argument('--test', type=str, help='Test single sample (sample ID)')
    parser.add_argument('--jobs', type=int, default=4, help='Parallel jobs (default: 4)')
    parser.add_argument('--list', action='store_true', help='List all available samples')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Coordinate Transformation Script v2.0 for variants_call_10")
    print(f"Offset: {DEFAULT_OFFSET}, Genome length: {GENOME_LENGTH}")
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
        
        result = process_single_sample(args.test, sample_dirs[args.test])
        if result:
            print(f"\nSample {args.test} transformation complete:")
            print(f"  Total variants: {result['total_variants']}")
            print(f"  Output: {Path(OUTPUT_DIR) / args.test}")
            for file_type, file_info in result['files_processed'].items():
                print(f"  - {file_type}: {file_info['variants_count']} variants")
    else:
        # Batch processing
        print(f"\nTransforming {len(sample_dirs)} samples (parallel: {args.jobs})")
        
        results = []
        sample_items = list(sample_dirs.items())
        
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {
                executor.submit(process_single_sample, sample_id, sample_info): sample_id
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
        
        generate_summary(results)
        
        successful = [r for r in results if r is not None]
        total_transformed = sum(r['total_variants'] for r in successful)
        
        print(f"\nTransformation complete!")
        print(f"  Successful: {len(successful)}/{len(sample_items)} samples")
        print(f"  Total variants: {total_transformed}")
        print(f"  Output: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
