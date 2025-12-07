#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch Individual Consensus Mapping Script for Nanopore Sequencing Data
=======================================================================

Purpose:
    Align nanopore reads to their individual consensus sequences using minimap2.
    Each sample is mapped to its own consensus sequence (not a universal reference).

Required Environment:
    - minimap2 (alignment tool)
    - samtools (BAM processing)
    
    Install via conda/mamba:
        conda create -n mapping minimap2 samtools
        conda activate mapping

Input:
    - FASTQ files: {FASTQ_DIR}/{sample_id}_subsampled.trimmed_filtered.porechop/{sample_id}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz
    - Consensus sequences: {CONSENSUS_DIR}/{sample_id}/consensus.fasta

Output:
    - {WORK_DIR}/mapping/{sample_id}/
      * {sample_id}.sorted.bam - Sorted alignment file
      * {sample_id}.sorted.bam.bai - BAM index
      * {sample_id}.sorted.bam.flagstat.txt - Alignment statistics

Adjustable Parameters:
    --threads: Number of threads per sample for minimap2 (default: 2)
    --processes: Number of samples to process in parallel (default: 8)

Usage:
    python 11.1_mapping_batch_individual.py [-t threads] [-p processes]

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import glob
import subprocess
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

# ========== PATH CONFIGURATION ==========
# These paths can be overridden by:
#   1. Environment variables (HBV_FASTQ_DIR, HBV_CONSENSUS_DIR, HBV_WORK_DIR)
#   2. Command line arguments
#
# Default: Use relative paths from script location (for standalone use)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))

# FASTQ files directory
FASTQ_DIR = os.environ.get(
    "HBV_FASTQ_DIR",
    os.path.join(_PROJECT_DIR, "host_deconv_out_5")
)

# Consensus sequences directory
CONSENSUS_DIR = os.environ.get(
    "HBV_CONSENSUS_DIR",
    os.path.join(_PROJECT_DIR, "Medaka_consensus_8", "r2")
)

# Working directory (output root)
WORK_DIR = os.environ.get(
    "HBV_WORK_DIR",
    os.path.join(_PROJECT_DIR, "variants_call_10")
)

# Output directory for mapping results
OUTPUT_DIR = os.path.join(WORK_DIR, "mapping")
# ========================================

def run_mapping(fq_path, ref_fa, out_dir, sample_id, threads=2):
    """
    Process a single fastq file:
      1. minimap2 alignment
      2. samtools convert to BAM, sort, and index
      3. Generate alignment statistics
    
    Args:
        fq_path: Path to input FASTQ file
        ref_fa: Path to reference FASTA (consensus sequence)
        out_dir: Output directory
        sample_id: Sample identifier
        threads: Number of threads for minimap2
    
    Returns:
        Path to sorted BAM file
    """
    # Create output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Output file naming
    sample_out_dir = os.path.join(out_dir, sample_id)
    if not os.path.exists(sample_out_dir):
        os.makedirs(sample_out_dir, exist_ok=True)
        
    sam_path = os.path.join(sample_out_dir, f"{sample_id}.sam")
    sorted_bam = os.path.join(sample_out_dir, f"{sample_id}.sorted.bam")
    sorted_bam_index = sorted_bam + ".bai"
    stats_txt = os.path.join(sample_out_dir, f"{sample_id}.sorted.bam.flagstat.txt")

    # minimap2 alignment (-ax map-ont -t THREADS)
    cmd_minimap2 = [
        "minimap2",
        "-ax", "map-ont",
        "-t", str(threads),
        ref_fa,
        fq_path
    ]
    
    print(f"[mapping] Running minimap2 alignment for {sample_id}...")
    with open(sam_path, "w") as sam_out:
        subprocess.run(cmd_minimap2, stdout=sam_out, check=True)

    # Convert SAM to sorted BAM
    print(f"[mapping] Converting SAM to sorted BAM for {sample_id}...")
    cmd_sam2bam = [
        "samtools", "view",
        "-bS", sam_path
    ]
    cmd_sort = [
        "samtools", "sort",
        "-o", sorted_bam
    ]
    # Use pipe
    p1 = subprocess.Popen(cmd_sam2bam, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd_sort, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p2.communicate()  # Wait for completion

    # Index BAM file
    print(f"[mapping] Creating index for sorted BAM {sample_id}...")
    cmd_index = ["samtools", "index", sorted_bam]
    subprocess.run(cmd_index, check=True)

    # Generate alignment statistics
    print(f"[mapping] Generating flagstat for {sample_id}...")
    with open(stats_txt, "w") as fstats:
        subprocess.run(["samtools", "flagstat", sorted_bam], stdout=fstats, check=True)

    # Remove SAM file to save space
    if os.path.exists(sam_path):
        os.remove(sam_path)

    print(f"[mapping] Completed: {sample_id}.sorted.bam")
    return sorted_bam

def find_samples():
    """
    Find all samples with consensus sequences and filtered reads.
    
    Returns:
        Dictionary: {sample_id: (consensus_path, reads_path)}
    
    Directory structure:
        FASTQ: {FASTQ_DIR}/{sample_id}_subsampled.trimmed_filtered.porechop/{sample_id}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz
        Consensus: {CONSENSUS_DIR}/{sample_id}/consensus.fasta
    """
    samples = {}
    
    # Check if FASTQ directory exists
    if not os.path.exists(FASTQ_DIR):
        print(f"Error: FASTQ directory not found: {FASTQ_DIR}")
        return samples
    
    # Check if Consensus directory exists
    if not os.path.exists(CONSENSUS_DIR):
        print(f"Error: Consensus directory not found: {CONSENSUS_DIR}")
        return samples
    
    # Find all consensus directories (each represents a sample)
    for item in os.listdir(CONSENSUS_DIR):
        consensus_dir = os.path.join(CONSENSUS_DIR, item)
        
        # Skip non-directories
        if not os.path.isdir(consensus_dir):
            continue
        
        sample_id = item
        
        # Check consensus file exists
        consensus_file = os.path.join(consensus_dir, "consensus.fasta")
        if not os.path.exists(consensus_file):
            print(f"Warning: Sample {sample_id} missing consensus file: {consensus_file}")
            continue
        
        # Find corresponding FASTQ file
        # Pattern: {sample_id}_subsampled.trimmed_filtered.porechop/{sample_id}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz
        fastq_subdir = f"{sample_id}_subsampled.trimmed_filtered.porechop"
        fastq_filename = f"{sample_id}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz"
        reads_file = os.path.join(FASTQ_DIR, fastq_subdir, fastq_filename)
        
        if not os.path.exists(reads_file):
            # Try alternative pattern
            alt_pattern = os.path.join(FASTQ_DIR, fastq_subdir, "*.fastq.gz")
            alt_files = glob.glob(alt_pattern)
            if alt_files:
                reads_file = alt_files[0]
            else:
                print(f"Warning: Sample {sample_id} missing FASTQ file: {reads_file}")
                continue
        
        # Check if already processed
        bam_file = os.path.join(OUTPUT_DIR, sample_id, f"{sample_id}.sorted.bam")
        bai_file = bam_file + ".bai"
        flagstat_file = bam_file + ".flagstat.txt"
        
        if os.path.exists(bam_file) and os.path.exists(bai_file) and os.path.exists(flagstat_file):
            print(f"Skip: Sample {sample_id} already processed")
        else:
            samples[sample_id] = (consensus_file, reads_file)
    
    return samples

def process_single_sample(sample_data):
    """
    Wrapper function to process a single sample for multiprocessing.
    
    Args:
        sample_data: Tuple of (sample_id, consensus_path, reads_path, output_dir, threads)
    
    Returns:
        Dictionary with processing result
    """
    sample_id, consensus_path, reads_path, output_dir, threads = sample_data
    
    try:
        print(f"\n[Process {os.getpid()}] Starting sample {sample_id}")
        print(f"[Process {os.getpid()}] Consensus sequence: {consensus_path}")
        print(f"[Process {os.getpid()}] Filtered reads: {reads_path}")
        
        out_bam = run_mapping(reads_path, consensus_path, output_dir, sample_id, threads)
        
        print(f"[Process {os.getpid()}] Sample {sample_id} completed. BAM saved to: {out_bam}")
        return {"sample_id": sample_id, "status": "success", "output": out_bam}
    
    except Exception as e:
        print(f"[Process {os.getpid()}] Error processing sample {sample_id}: {str(e)}")
        return {"sample_id": sample_id, "status": "error", "error": str(e)}

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Batch process samples, aligning each to its own consensus sequence"
    )
    parser.add_argument(
        "-t", "--threads", 
        type=int, 
        default=2, 
        help="Number of threads for minimap2 per sample (default: 2)"
    )
    parser.add_argument(
        "-p", "--processes", 
        type=int, 
        default=8, 
        help="Number of samples to process in parallel (default: 8)"
    )
    args = parser.parse_args()
    
    threads = args.threads
    max_processes = args.processes
    
    print("=" * 60)
    print("Mapping Script v2.0 for variants_call_10")
    print("=" * 60)
    print(f"FASTQ directory: {FASTQ_DIR}")
    print(f"Consensus directory: {CONSENSUS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Each sample uses {threads} threads for minimap2")
    print(f"Processing {max_processes} samples in parallel")
    print(f"Total CPU core usage: {threads * max_processes} cores")
    print("=" * 60)
    
    # Create output directory
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find all samples
    samples = find_samples()
    
    if not samples:
        print("No matching samples found! Please check the data directory.")
        return
    
    print(f"Found {len(samples)} samples to process.\n")
    
    # Prepare sample data for multiprocessing
    sample_tasks = []
    for sample_id, (consensus_path, reads_path) in samples.items():
        sample_tasks.append((sample_id, consensus_path, reads_path, OUTPUT_DIR, threads))
    
    # Process samples using multiprocessing
    successful = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=max_processes) as executor:
        print(f"Starting {max_processes} parallel processes to handle samples...")
        results = executor.map(process_single_sample, sample_tasks)
        
        for result in results:
            if result["status"] == "success":
                successful += 1
                print(f"[OK] Sample {result['sample_id']} processed successfully")
            else:
                failed += 1
                print(f"[FAIL] Sample {result['sample_id']} failed: {result['error']}")
    
    print(f"\n=== Processing Complete Summary ===")
    print(f"Successfully processed: {successful} samples")
    print(f"Failed: {failed} samples")
    print(f"Total: {successful + failed} samples")

if __name__ == "__main__":
    # Set multiprocessing start method to avoid issues on some systems
    mp.set_start_method('spawn', force=True)
    main()
