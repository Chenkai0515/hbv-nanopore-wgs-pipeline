#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combined BAM Filter Script - Batch Processing Mode
===================================================

Purpose:
    Filter BAM files to remove low-quality alignments and reads.
    Supports both single sample and batch processing modes.
    
    Filtering steps:
    1. Remove unmapped and secondary alignments
    2. Remove reads with large soft clips or mismatches at start/end
    3. Remove reads shorter than minimum length
    
    Note: No PCR duplicate removal is performed.

Required Environment:
    - samtools (BAM file manipulation)
    - pysam (Python BAM file reading)
    
    Install via conda/mamba:
        conda create -n filtering samtools pysam
        conda activate filtering

Input:
    - BAM files from mapping step: {WORK_DIR}/mapping/{sample_id}/{sample_id}.sorted.bam

Output:
    - Filtered BAM files: {WORK_DIR}/filtering/{sample_id}/
      * {sample_id}.final.filtered.sorted.bam - Final filtered and sorted BAM
      * {sample_id}.final.filtered.sorted.bam.bai - BAM index
      * {sample_id}.final.filtered.sorted.bam.flagstat.txt - Statistics

Adjustable Parameters:
    --mismatch_thresh: Maximum allowed soft clip or mismatch length at start/end (default: 50)
    --filter_end: Also filter reads with mismatches at the end (default: False)
    --min_length: Minimum read length in bp (default: 500)

Usage:
    # Single sample
    python 11.2combined_bam_filter_batch.py --sample 10090
    
    # All samples
    python 11.2combined_bam_filter_batch.py --all
    
    # With custom parameters
    python 11.2combined_bam_filter_batch.py --all --mismatch_thresh 60 --min_length 600

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import sys
import subprocess
import re
import argparse
import pysam
import glob
from functools import partial
from concurrent.futures import ProcessPoolExecutor

# ========== PATH CONFIGURATION ==========
# Working directory (same as mapping output)
WORK_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/variants_call_10"

# Input directory (mapping results)
INPUT_DIR = os.path.join(WORK_DIR, "mapping")

# Output directory (filtered results)
OUTPUT_DIR = os.path.join(WORK_DIR, "filtering")
# ========================================

def filter_standard(input_bam, output_dir):
    """
    Standard BAM filtering: filter out unmapped/secondary alignments.
    No PCR duplicate removal.
    
    Args:
        input_bam: Path to input BAM file
        output_dir: Output directory
    
    Returns:
        Path to filtered BAM file
    """
    base_name = os.path.basename(input_bam)
    sample_id = base_name.replace(".sorted.bam", "").replace(".bam", "")

    # Create sample output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Output files
    final_bam = os.path.join(sample_output_dir, f"{sample_id}.filtered.bam")
    final_index = final_bam + ".bai"
    stats_txt = os.path.join(sample_output_dir, f"{sample_id}.filtered.bam.flagstat.txt")

    # Filter -F 4(unmapped) -F 256(secondary) then sort and index
    print(f"[Filtering] Removing unmapped and secondary alignments for {sample_id}...")
    cmd_view = ["samtools", "view", "-b", "-F", "4", "-F", "256", input_bam]
    cmd_sort = ["samtools", "sort", "-o", final_bam, "-"]
    p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd_sort, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()

    # Index
    print(f"[Filtering] Creating index for filtered BAM {sample_id}...")
    subprocess.run(["samtools", "index", final_bam], check=True)

    # Generate flagstat statistics
    print(f"[Filtering] Generating flagstat for {sample_id}...")
    with open(stats_txt, "w") as fstats:
        subprocess.run(["samtools", "flagstat", final_bam], stdout=fstats, check=True)

    print(f"[Filtering] Completed: {sample_id}.filtered.bam")
    return final_bam

def analyze_soft_clip_at_start(cigar_string):
    """
    Analyze the length of soft clip at the start of the CIGAR string.
    
    Args:
        cigar_string: CIGAR string from BAM alignment
    
    Returns:
        Length of soft clip at start (0 if none)
    """
    if not cigar_string:
        return 0
        
    # Extract operations from CIGAR
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    # Check if the first operation is soft clip (S)
    if cigar_tuples and cigar_tuples[0][1] == 'S':
        return int(cigar_tuples[0][0])
    return 0

def analyze_soft_clip_at_end(cigar_string):
    """
    Analyze the length of soft clip at the end of the CIGAR string.
    
    Args:
        cigar_string: CIGAR string from BAM alignment
    
    Returns:
        Length of soft clip at end (0 if none)
    """
    if not cigar_string:
        return 0
        
    # Extract operations from CIGAR
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    # Check if the last operation is soft clip (S)
    if cigar_tuples and cigar_tuples[-1][1] == 'S':
        return int(cigar_tuples[-1][0])
    return 0

def analyze_md_tag(md_tag, read_seq, head_length=50, tail_length=50, check_tail=False):
    """
    Analyze MD tag to determine if start/end of read has many mismatches.
    
    Args:
        md_tag: MD tag from BAM alignment
        read_seq: Read sequence
        head_length: Length to check at start
        tail_length: Length to check at end
        check_tail: Whether to check tail instead of head
    
    Returns:
        Number of mismatches in the specified region
    """
    if not md_tag or not read_seq:
        return 0
        
    # Parse MD tag, identify match and mismatch positions
    md_parts = re.findall(r'(\d+|\^[ACGTN]+|[ACGTN])', md_tag)
    
    # Track reference sequence position
    ref_pos = 0
    # Track mismatch positions in read sequence
    mismatches = []
    
    for part in md_parts:
        if part.isdigit():  # Match segment
            ref_pos += int(part)
        elif part.startswith('^'):  # Deletion
            # Deletions don't count as mismatches
            ref_pos += len(part) - 1  # Subtract '^'
        else:  # Single base mismatch
            mismatches.append(ref_pos)
            ref_pos += 1
    
    if check_tail:
        # Count mismatches in the tail of tail_length
        read_length = len(read_seq)
        tail_mismatches = sum(1 for pos in mismatches if pos >= read_length - tail_length)
        return tail_mismatches
    else:
        # Count mismatches in the head of head_length
        head_mismatches = sum(1 for pos in mismatches if pos < head_length)
        return head_mismatches

def filter_bam_by_mismatch(input_bam, output_bam, mismatch_thresh=50, filter_end=False, min_length=500):
    """
    Filter BAM file, discarding:
    1. Reads with long unmapped segments at start/end
    2. Reads shorter than min_length
    
    Args:
        input_bam: Input BAM file path
        output_bam: Output BAM file path
        mismatch_thresh: Maximum allowed soft clip or mismatch length
        filter_end: Whether to also filter reads with mismatches at the end
        min_length: Minimum read length
    
    Returns:
        Dictionary with filtering statistics
    """
    # Open input BAM file
    infile = pysam.AlignmentFile(input_bam, "rb")
    
    # Create output BAM file
    outfile = pysam.AlignmentFile(output_bam, "wb", template=infile)
    
    # Statistics counters
    total_reads = 0
    kept_reads = 0
    filtered_by_soft_clip_start = 0
    filtered_by_md_start = 0
    filtered_by_soft_clip_end = 0
    filtered_by_md_end = 0
    filtered_by_length = 0
    
    # Process each read
    for read in infile:
        total_reads += 1
        
        # Skip unmapped reads
        if read.is_unmapped:
            continue
        
        # 0. Check read length
        read_length = len(read.query_sequence) if read.query_sequence else 0
        if read_length < min_length:
            filtered_by_length += 1
            continue
        
        # 1. Check soft clip length at start
        soft_clip_start = analyze_soft_clip_at_start(read.cigarstring)
        if soft_clip_start >= mismatch_thresh:
            filtered_by_soft_clip_start += 1
            continue
        
        # 2. Check MD tag for mismatches at start
        md_tag = read.get_tag('MD') if read.has_tag('MD') else ""
        read_seq = read.query_sequence or ""
        
        if md_tag and read_seq:
            head_mismatches = analyze_md_tag(md_tag, read_seq, mismatch_thresh)
            # Filter if most of the start is mismatches
            if head_mismatches > mismatch_thresh * 0.8:  # >80% are mismatches
                filtered_by_md_start += 1
                continue
        
        # If end checking is requested
        if filter_end:
            # 3. Check soft clip length at end
            soft_clip_end = analyze_soft_clip_at_end(read.cigarstring)
            if soft_clip_end >= mismatch_thresh:
                filtered_by_soft_clip_end += 1
                continue
            
            # 4. Check MD tag for mismatches at end
            if md_tag and read_seq:
                tail_mismatches = analyze_md_tag(md_tag, read_seq, mismatch_thresh, mismatch_thresh, True)
                # Filter if most of the end is mismatches
                if tail_mismatches > mismatch_thresh * 0.8:  # >80% are mismatches
                    filtered_by_md_end += 1
                    continue
        
        # Passed all checks, keep this read
        outfile.write(read)
        kept_reads += 1
    
    # Close files
    infile.close()
    outfile.close()
    
    return {
        "total": total_reads,
        "kept": kept_reads,
        "filtered_by_length": filtered_by_length,
        "filtered_by_soft_clip_start": filtered_by_soft_clip_start,
        "filtered_by_md_start": filtered_by_md_start,
        "filtered_by_soft_clip_end": filtered_by_soft_clip_end if filter_end else 0,
        "filtered_by_md_end": filtered_by_md_end if filter_end else 0
    }

def process_single_bam(input_bam, mismatch_thresh=50, filter_end=False, min_length=500):
    """
    Process a single BAM file through all filtering steps.
    
    Args:
        input_bam: Path to input BAM file
        mismatch_thresh: Mismatch threshold
        filter_end: Filter end regions
        min_length: Minimum read length
    
    Returns:
        True if successful, False otherwise
    """
    try:
        base_name = os.path.basename(input_bam)
        sample_id = base_name.replace(".sorted.bam", "").replace(".bam", "")
        
        print(f"\n{'='*60}")
        print(f"Processing sample: {sample_id}")
        print(f"Input BAM: {input_bam}")
        print(f"Filter parameters: mismatch_thresh={mismatch_thresh}, filter_end={filter_end}, min_length={min_length}bp")
        print(f"{'='*60}")
        
        # Phase 1: Standard filtering
        print("\n=== Phase 1: Standard Filtering ===")
        standard_filtered_bam = filter_standard(input_bam, OUTPUT_DIR)
        
        # Phase 2: Mismatch/soft clip and length-based filtering
        print("\n=== Phase 2: Mismatch/Soft Clip and Length-based Filtering ===")
        sample_output_dir = os.path.join(OUTPUT_DIR, sample_id)
        final_bam = os.path.join(sample_output_dir, f"{sample_id}.final.filtered.bam")
        
        filter_description = "start" if not filter_end else "start and end"
        print(f"Filtering {standard_filtered_bam}, removing reads with long mismatches at {filter_description} (threshold: {mismatch_thresh} bp)...")
        print(f"Also filtering out reads shorter than {min_length}bp...")
        
        stats = filter_bam_by_mismatch(standard_filtered_bam, final_bam, mismatch_thresh, filter_end, min_length)
        
        print(f"Total reads: {stats['total']}")
        print(f"Kept reads: {stats['kept']}")
        print(f"Reads filtered by length < {min_length}bp: {stats['filtered_by_length']}")
        print(f"Reads filtered by start soft clip: {stats['filtered_by_soft_clip_start']}")
        print(f"Reads filtered by start MD tag mismatches: {stats['filtered_by_md_start']}")
        
        if filter_end:
            print(f"Reads filtered by end soft clip: {stats['filtered_by_soft_clip_end']}")
            print(f"Reads filtered by end MD tag mismatches: {stats['filtered_by_md_end']}")
        
        print(f"Total filtered reads: {stats['total'] - stats['kept']}")
        if stats['total'] > 0:
            print(f"Retention rate: {stats['kept']/stats['total']*100:.2f}%")
        else:
            print("Retention rate: N/A (no reads found)")
        
        # Sort and index final output
        print(f"\n=== Creating Index for Final Output File ===")
        final_sorted_bam = os.path.join(sample_output_dir, f"{sample_id}.final.filtered.sorted.bam")
        
        subprocess.run(["samtools", "sort", "-o", final_sorted_bam, final_bam], check=True)
        subprocess.run(["samtools", "index", final_sorted_bam], check=True)
        
        # Generate flagstat statistics
        final_stats = os.path.join(sample_output_dir, f"{sample_id}.final.filtered.sorted.bam.flagstat.txt")
        with open(final_stats, "w") as fstats:
            subprocess.run(["samtools", "flagstat", final_sorted_bam], stdout=fstats, check=True)
        
        print("\n=== Cleaning Temporary Files ===")
        # Delete intermediate files
        if os.path.exists(final_bam):
            os.remove(final_bam)
        if os.path.exists(standard_filtered_bam):
            os.remove(standard_filtered_bam)
        bai = standard_filtered_bam + ".bai"
        if os.path.exists(bai):
            os.remove(bai)
        flagstat = standard_filtered_bam + ".flagstat.txt"
        if os.path.exists(flagstat):
            os.remove(flagstat)
        
        print(f"\nProcessing complete! Final results saved to: {final_sorted_bam}")
        return True
        
    except Exception as e:
        print(f"Error processing {input_bam}: {str(e)}")
        return False

def find_all_bam_files():
    """
    Find all BAM files in the mapping directory.
    
    Returns:
        List of BAM file paths
    """
    bam_files = []
    
    if not os.path.exists(INPUT_DIR):
        print(f"Error: Input directory not found: {INPUT_DIR}")
        return bam_files
    
    # Find all sorted.bam files in subdirectories
    for sample_dir in os.listdir(INPUT_DIR):
        sample_path = os.path.join(INPUT_DIR, sample_dir)
        if not os.path.isdir(sample_path):
            continue
        
        bam_file = os.path.join(sample_path, f"{sample_dir}.sorted.bam")
        if os.path.exists(bam_file):
            # Check if already processed
            output_sample_dir = os.path.join(OUTPUT_DIR, sample_dir)
            final_bam = os.path.join(output_sample_dir, f"{sample_dir}.final.filtered.sorted.bam")
            
            if os.path.exists(final_bam):
                print(f"Skip: Sample {sample_dir} already filtered")
            else:
                bam_files.append(bam_file)
    
    return sorted(bam_files)

def main():
    parser = argparse.ArgumentParser(
        description='Combined BAM filter tool: filter unmapped/secondary reads, filter mismatched starts/ends, and filter short reads (no PCR deduplication)'
    )
    
    # Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--all', action='store_true', help='Process all samples in batch mode')
    mode_group.add_argument('--sample', type=str, help='Process a single sample by ID')
    
    # Filter parameters
    parser.add_argument(
        '--mismatch_thresh', 
        type=int, 
        default=50,
        help='Maximum allowed soft clip or mismatch length at start/end (default: 50)'
    )
    parser.add_argument(
        '--filter_end', 
        action='store_true', 
        help='Also filter the end regions'
    )
    parser.add_argument(
        '--min_length', 
        type=int, 
        default=500,
        help='Minimum read length (default: 500)'
    )
    parser.add_argument(
        '--jobs', 
        '-j', 
        type=int, 
        default=1,
        help='Number of samples to process in parallel (default: 1)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("="*60)
    print("Combined BAM Filter Script v2.0 for variants_call_10")
    print("="*60)
    print(f"Input directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Filter parameters: mismatch_thresh={args.mismatch_thresh}, filter_end={args.filter_end}, min_length={args.min_length}bp")
    print("="*60)
    
    if args.sample:
        # Single sample mode
        print(f"\nSingle sample mode: {args.sample}")
        input_bam = os.path.join(INPUT_DIR, args.sample, f"{args.sample}.sorted.bam")
        
        if not os.path.exists(input_bam):
            print(f"Error: Input BAM file not found: {input_bam}")
            sys.exit(1)
        
        success = process_single_bam(input_bam, args.mismatch_thresh, args.filter_end, args.min_length)
        sys.exit(0 if success else 1)
        
    elif args.all:
        # Batch mode
        print("\nBatch processing mode: all samples")
        bam_files = find_all_bam_files()
        
        if not bam_files:
            print("No BAM files found to process!")
            sys.exit(0)
        
        print(f"Found {len(bam_files)} samples to process")
        
        # Process with multiprocessing if jobs > 1
        if args.jobs > 1:
            print(f"Processing {args.jobs} samples in parallel...")
            
            # Use partial to create a picklable function with fixed arguments
            process_func = partial(process_single_bam, 
                                   mismatch_thresh=args.mismatch_thresh, 
                                   filter_end=args.filter_end, 
                                   min_length=args.min_length)
            
            with ProcessPoolExecutor(max_workers=args.jobs) as executor:
                results = list(executor.map(process_func, bam_files))
            
            success_count = sum(results)
            failed_count = len(results) - success_count
        else:
            # Serial processing
            success_count = 0
            failed_count = 0
            
            for bam_file in bam_files:
                if process_single_bam(bam_file, args.mismatch_thresh, args.filter_end, args.min_length):
                    success_count += 1
                else:
                    failed_count += 1
        
        print(f"\n{'='*60}")
        print("Batch Processing Summary")
        print(f"{'='*60}")
        print(f"Total samples: {len(bam_files)}")
        print(f"Successfully processed: {success_count}")
        print(f"Failed: {failed_count}")
        print(f"{'='*60}")

if __name__ == "__main__":
    main()
