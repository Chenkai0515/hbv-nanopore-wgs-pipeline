#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
iVar Variant Calling Script - Individual Consensus Reference
=============================================================

Purpose:
    Perform variant calling using iVar on filtered BAM files.
    Each sample is called against its own consensus sequence (not a universal reference).

Required Environment:
    - samtools (for mpileup generation)
    - ivar (for variant calling)
    
    Install via conda/mamba:
        conda create -n ivar samtools ivar
        conda activate ivar

Input:
    - Filtered BAM files: {WORK_DIR}/filtering/{sample_id}/{sample_id}.final.filtered.sorted.bam
    - Consensus sequences: {CONSENSUS_DIR}/{sample_id}/consensus.fasta

Output:
    - Variant files: {WORK_DIR}/variants/ivar/{sample_id}/
      * variants.tsv - Main variant calls (SNPs and indels)
      * variants.tsv.indel - Indel-specific calls (if any)

Adjustable Parameters:
    MIN_DEPTH: Minimum sequencing depth for variant calling (default: 100)
    MIN_FREQ: Minimum variant allele frequency threshold (default: 0.01)
    --threads: Reserved thread parameter (mpileup is single-threaded)
    --jobs: Number of samples to process in parallel (default: 1)
    --max-depth: Maximum mpileup depth limit (default: 0 for unlimited)
    --baseq: Minimum base quality for mpileup (default: 15)
    --mapq: Minimum mapping quality for mpileup (default: 10)

Usage:
    # Test single sample
    python 11.3a_variant_call_ivar_batch.py --test 10090
    
    # List all available samples
    python 11.3a_variant_call_ivar_batch.py --list
    
    # Process all samples (with confirmation)
    python 11.3a_variant_call_ivar_batch.py --all
    
    # Process all samples without confirmation
    python 11.3a_variant_call_ivar_batch.py --all -y
    
    # Parallel processing with 4 jobs
    python 11.3a_variant_call_ivar_batch.py --all -j 4

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import subprocess
import multiprocessing as mp
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

# ========== PATH CONFIGURATION ==========
# Working directory
WORK_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/variants_call_10"

# Consensus sequences directory
CONSENSUS_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r2"

# BAM files directory (from filtering step)
BAM_DIR = os.path.join(WORK_DIR, "filtering")

# Output directory
OUTPUT_DIR = os.path.join(WORK_DIR, "variants", "ivar")
# ========================================

# Variant calling parameters
MIN_DEPTH = 100   # ivar variants -m parameter (minimum sequencing depth)
MIN_FREQ  = 0.01  # ivar variants -t parameter (variant frequency threshold)

# Performance-related default parameters (can be overridden via command line)
MP_THREADS_DEFAULT = min(8, mp.cpu_count())
MP_MAX_DEPTH_DEFAULT = 0      # For high-coverage viral samples, unlimited depth (0 = unlimited)
MP_BASE_QUAL_DEFAULT = 15     # Minimum base quality (-Q 15)
MP_MAP_QUAL_DEFAULT  = 10     # Minimum mapping quality (default 10)

def run_variant_call(
    sample_id: str,
    mp_threads: int = MP_THREADS_DEFAULT,
    mp_max_depth: int = MP_MAX_DEPTH_DEFAULT,
    mp_base_qual: int = MP_BASE_QUAL_DEFAULT,
    mp_map_qual: int  = MP_MAP_QUAL_DEFAULT,
    force: bool = False,
    skip_existing: bool = False,
):
    """
    Perform variant calling on a filtered BAM file using its consensus sequence as reference.
    Supports configurable mpileup parameters for improved speed.
    
    Args:
        sample_id: Sample identifier
        mp_threads: Thread parameter (reserved, mpileup doesn't support multithreading)
        mp_max_depth: Maximum mpileup depth
        mp_base_qual: Minimum base quality
        mp_map_qual: Minimum mapping quality
        force: Force rerun even if output exists
        skip_existing: Skip if output already exists
    
    Returns:
        True if successful, False otherwise
    """
    # Sample-related paths
    bam_path = os.path.join(BAM_DIR, sample_id, f"{sample_id}.final.filtered.sorted.bam")
    consensus_fasta = os.path.join(CONSENSUS_DIR, sample_id, "consensus.fasta")
    out_dir = os.path.join(OUTPUT_DIR, sample_id)
    os.makedirs(out_dir, exist_ok=True)

    # Define output file names
    variants_tsv = os.path.join(out_dir, "variants.tsv")
    indel_tsv = variants_tsv.replace('.tsv', '.tsv.indel')

    # Skip if output already exists
    if skip_existing and os.path.exists(variants_tsv) and os.path.getsize(variants_tsv) > 0:
        print(f"[{sample_id}] Output already exists {variants_tsv}, skipping per --skip-existing")
        return True
    if not force and os.path.exists(variants_tsv) and os.path.getsize(variants_tsv) > 0:
        print(f"[{sample_id}] Detected existing output. Use --force to rerun or delete existing results")
        # Still return True, consider as completed
        return True

    # Check if files exist
    if not os.path.exists(bam_path):
        print(f"[Warning] BAM file does not exist: {bam_path}, skipping sample {sample_id}")
        return False

    if not os.path.exists(consensus_fasta):
        print(f"[Warning] Consensus sequence file does not exist: {consensus_fasta}, skipping sample {sample_id}")
        return False

    print(f"=== [Sample: {sample_id}] Starting variant calling ===")
    print(f"[{sample_id}] BAM file: {bam_path}")
    print(f"[{sample_id}] Reference sequence: {consensus_fasta}")

    # 1) Check if BAM index exists, create if not
    bam_index = bam_path + ".bai"
    if not os.path.exists(bam_index) or (os.path.getsize(bam_index) == 0):
        print(f"[{sample_id}] Creating BAM index...")
        try:
            subprocess.run(["samtools", "index", bam_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"[Error] Failed to create BAM index for sample {sample_id}: {e}")
            return False
    else:
        print(f"[{sample_id}] BAM index exists, skipping index creation")

    # 2) Check if consensus sequence has index, create if not
    fai_path = consensus_fasta + ".fai"
    if not os.path.exists(fai_path) or (os.path.getsize(fai_path) == 0):
        print(f"[{sample_id}] Creating consensus sequence index...")
        try:
            subprocess.run(["samtools", "faidx", consensus_fasta], check=True)
        except subprocess.CalledProcessError as e:
            print(f"[Error] Failed to create FASTA index for sample {sample_id}: {e}")
            return False

    # 3) samtools mpileup + ivar variants
    # Note: samtools mpileup doesn't support --threads parameter, only single-threaded
    cmd_mpileup = [
        "samtools", "mpileup",
        "-A",
        "-d", str(int(mp_max_depth)) if int(mp_max_depth) > 0 else "0",
        "-Q", str(int(mp_base_qual)),  # Changed to -Q 15 (default), filters low base quality to speed up
        "-B",
        "-f", consensus_fasta,
    ]
    # Optional: minimum mapping quality
    if int(mp_map_qual) > 0:
        cmd_mpileup.extend(["-q", str(int(mp_map_qual))])
    # Target BAM
    cmd_mpileup.append(bam_path)

    cmd_ivar = [
        "ivar", "variants",
        "-r", consensus_fasta,
        "-m", str(MIN_DEPTH),
        "-t", str(MIN_FREQ),
        "-p", variants_tsv.replace(".tsv","")  # ivar will automatically add .tsv / .tsv.indel
    ]

    print(f"[{sample_id}] Running variant calling via iVar...")
    print(f"[{sample_id}] mpileup command: {' '.join(cmd_mpileup)}")
    print(f"[{sample_id}] ivar command: {' '.join(cmd_ivar)}")

    # Use pipe to pass mpileup output to ivar
    try:
        # Note: Don't capture p1's stderr to avoid blocking due to many warnings; let it output directly to parent process stderr
        p1 = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd_ivar, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        stdout, stderr = p2.communicate()  # Wait for ivar to finish
        p1_ret = p1.wait()

        if p1_ret != 0:
            print(f"[Error] samtools mpileup failed for sample {sample_id}")
            return False
        if p2.returncode != 0:
            print(f"[Error] ivar variants failed for sample {sample_id}. Error message: {stderr.decode('utf-8', errors='ignore')}")
            return False

        # Check output results
        if os.path.exists(variants_tsv) and os.path.getsize(variants_tsv) > 0:
            # Count variants
            variant_count = 0
            with open(variants_tsv, 'r') as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    if i == 0:
                        pass
                    variant_count += 1
            # Subtract header (1 line)
            variant_count = max(0, variant_count - 1)
            print(f"[{sample_id}] Complete! Detected {variant_count} variants")
            print(f"[{sample_id}] Variant results saved to: {variants_tsv}")

            # If indel file exists, report it too
            if os.path.exists(indel_tsv) and os.path.getsize(indel_tsv) > 0:
                indel_count = 0
                with open(indel_tsv, 'r') as f:
                    for i, line in enumerate(f):
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        indel_count += 1
                indel_count = max(0, indel_count - 1)
                print(f"[{sample_id}] Indel file: {indel_count} insertion/deletion variants")
        else:
            print(f"[{sample_id}] Warning: No variant results file generated")
            return False

        return True

    except Exception as e:
        print(f"[Error] Exception occurred during variant calling for sample {sample_id}: {e}")
        return False


def get_sample_list():
    """
    Get list of all sample IDs (only keep samples with both BAM and consensus sequence).
    
    Returns:
        List of valid sample IDs
    """
    sample_dirs = []
    if os.path.exists(BAM_DIR):
        sample_dirs = [d for d in os.listdir(BAM_DIR) if os.path.isdir(os.path.join(BAM_DIR, d))]

    valid_samples = []
    for sample_id in sample_dirs:
        bam_path = os.path.join(BAM_DIR, sample_id, f"{sample_id}.final.filtered.sorted.bam")
        consensus_path = os.path.join(CONSENSUS_DIR, sample_id, "consensus.fasta")
        if os.path.exists(bam_path) and os.path.exists(consensus_path):
            valid_samples.append(sample_id)
        else:
            print(f"[Warning] Sample {sample_id} missing files, skipping")
            if not os.path.exists(bam_path):
                print(f"  - Missing BAM: {bam_path}")
            if not os.path.exists(consensus_path):
                print(f"  - Missing consensus: {consensus_path}")

    return sorted(valid_samples)


def test_single_sample(sample_id, args):
    """
    Test a single sample.
    
    Args:
        sample_id: Sample ID to test
        args: Command line arguments
    
    Returns:
        True if successful, False otherwise
    """
    print(f"\n=== Single sample test mode: Sample {sample_id} ===")

    # Check if sample is in available list
    valid_samples = get_sample_list()
    if sample_id not in valid_samples:
        print(f"[Error] Sample {sample_id} not in available sample list")
        print(f"Available samples: {', '.join(valid_samples[:10])}{'...' if len(valid_samples) > 10 else ''}")
        return False

    # Run variant calling
    success = run_variant_call(
        sample_id,
        mp_threads=args.threads,
        mp_max_depth=args.max_depth,
        mp_base_qual=args.baseq,
        mp_map_qual=args.mapq,
        force=args.force,
        skip_existing=args.skip_existing
    )

    if success:
        print(f"\n[OK] Sample {sample_id} test successful!")

        # Display result summary
        variants_tsv = os.path.join(OUTPUT_DIR, sample_id, "variants.tsv")
        if os.path.exists(variants_tsv):
            print(f"\nResult summary:")
            print(f"Variant file: {variants_tsv}")

            # Simply display first few lines
            print("\nFirst 5 variant lines:")
            with open(variants_tsv, 'r') as f:
                for i, line in enumerate(f):
                    if i < 6:  # header + 5 data lines
                        print(f"  {line.rstrip()}")
                    else:
                        break
    else:
        print(f"\n[FAIL] Sample {sample_id} test failed!")

    return success


def main():
    """
    Main function: Support batch processing and single sample testing.
    """
    parser = argparse.ArgumentParser(
        description="ivar variant calling script - using each sample's own consensus sequence as reference"
    )
    parser.add_argument("--test", type=str, help="Single sample test mode, specify sample ID")
    parser.add_argument("--list", action="store_true", help="List all available samples")
    parser.add_argument("--all", action="store_true", help="Process all samples in batch mode")

    # Performance and parameters
    parser.add_argument("--threads", type=int, default=MP_THREADS_DEFAULT, help="Reserved thread parameter (samtools mpileup doesn't support multithreading)")
    parser.add_argument("--jobs", "-j", type=int, default=1, help="Number of samples to process in parallel (default: 1)")
    parser.add_argument("--max-depth", "-D", type=int, default=MP_MAX_DEPTH_DEFAULT, help="mpileup maximum depth limit (-d, default 0 = unlimited)")
    parser.add_argument("--baseq", "-Q", type=int, default=MP_BASE_QUAL_DEFAULT, help="mpileup minimum base quality (-Q, default 15)")
    parser.add_argument("--mapq", "-q", type=int, default=MP_MAP_QUAL_DEFAULT, help="mpileup minimum mapping quality (-q, default 10)")
    parser.add_argument("-y", "--yes", action="store_true", help="Skip confirmation in batch mode")
    parser.add_argument("--force", action="store_true", help="Force rerun (even if output exists)")
    parser.add_argument("--skip-existing", action="store_true", help="Skip sample if output already exists")

    args = parser.parse_args()

    # Create output directory
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("iVar Variant Calling Script v2.0 for variants_call_10")
    print("=" * 60)
    print(f"BAM directory: {BAM_DIR}")
    print(f"Consensus sequence directory: {CONSENSUS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Minimum depth: {MIN_DEPTH}, Minimum frequency: {MIN_FREQ}")
    print(f"mpileup parameters: -Q {args.baseq}, -q {args.mapq}, -d {args.max_depth if args.max_depth>0 else 'unlimited'}")
    print(f"Parallel: --jobs {args.jobs}")
    print("=" * 60)

    # Get sample list
    sample_ids = get_sample_list()
    if not sample_ids:
        print(f"No valid samples found!")
        return

    print(f"\nFound {len(sample_ids)} valid samples")

    # List samples mode
    if args.list:
        print("\nAvailable sample list:")
        for i, sample_id in enumerate(sample_ids, 1):
            print(f"{i:3d}. {sample_id}")
        return

    # Single sample test mode
    if args.test:
        success = test_single_sample(args.test, args)
        return success
    
    # Batch processing mode
    if args.all:
        print(f"\nPreparing to process samples: {', '.join(sample_ids[:5])}{'...' if len(sample_ids) > 5 else ''}")

        # Confirm whether to continue
        if not args.yes:
            response = input(f"\nContinue processing all {len(sample_ids)} samples? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Batch processing canceled")
                return

        # Execute: Support multi-sample parallel processing (mpileup within each sample also has multithreading)
        total = len(sample_ids)
        success_count = 0
        fail_count = 0
        finished = 0

        if args.jobs <= 1:
            # Serial
            print(f"\nDetected CPU core count: {mp.cpu_count()}. Using serial processing (--jobs=1).")
            for sample_id in sample_ids:
                print(f"\n{'='*50}")
                success = run_variant_call(
                    sample_id,
                    mp_threads=args.threads,
                    mp_max_depth=args.max_depth,
                    mp_base_qual=args.baseq,
                    mp_map_qual=args.mapq,
                    force=args.force,
                    skip_existing=args.skip_existing
                )
                if success:
                    success_count += 1
                else:
                    fail_count += 1
                finished += 1
                print(f"Progress: {finished}/{total} (success: {success_count}, failed: {fail_count})")
        else:
            # Parallel
            print(f"\nDetected CPU core count: {mp.cpu_count()}. Parallel processing (--jobs={args.jobs}).")
            # Hint: mpileup is single-threaded, mainly parallelizing multiple samples
            print(f"[Info] Will process {args.jobs} samples simultaneously, each sample's mpileup is single-threaded.")

            with ThreadPoolExecutor(max_workers=args.jobs) as ex:
                future_to_sample = {
                    ex.submit(
                        run_variant_call,
                        sid,
                        args.threads,
                        args.max_depth,
                        args.baseq,
                        args.mapq,
                        args.force,
                        args.skip_existing
                    ): sid for sid in sample_ids
                }
                for fut in as_completed(future_to_sample):
                    sid = future_to_sample[fut]
                    ok = False
                    try:
                        ok = fut.result()
                    except Exception as e:
                        print(f"[Error] Sample {sid} raised exception during parallel execution: {e}")
                        ok = False

                    if ok:
                        success_count += 1
                    else:
                        fail_count += 1
                    finished += 1
                    print(f"Progress: {finished}/{total} (success: {success_count}, failed: {fail_count})")

        print(f"\nAll variant calling tasks completed!")
        print(f"Processing statistics:")
        print(f"  - Total samples: {total}")
        print(f"  - Success: {success_count}")
        print(f"  - Failed: {fail_count}")
        print(f"  - Results saved to: {OUTPUT_DIR}")
    else:
        print("\n[Error] Please specify operation mode: --test, --list, or --all")
        parser.print_help()

if __name__ == "__main__":
    main()
