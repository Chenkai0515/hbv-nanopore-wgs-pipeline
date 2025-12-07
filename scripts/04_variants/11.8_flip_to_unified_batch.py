#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Flip View B Variants to Unified Reference Perspective
======================================================

Purpose:
    Transform variants from "View B" (reads vs consensus) to unified reference perspective.
    This final step enables cohort-level comparison by placing all samples on the same
    reference coordinate system.
    
    Critical Algorithm:
    - For positions where consensus (C) equals unified reference (R): Keep variants as-is
    - For positions where C != R: Flip ref/alt to unified perspective
    - INDELs: Pass-through unchanged (alignment-dependent)

Required Environment:
    - pandas (data manipulation)
    - numpy (numerical operations)
    
    Install via conda/mamba:
        conda activate ivar  # or any environment with pandas and numpy

Input:
    - View A (Consensus vs Unified Ref): {VIEW_A_DIR}/*.csv
      * Columns: pos, mutat, ref (=R), alt (=C), length
      * Contains positions where consensus differs from unified reference
    
    - View B (Reads vs Consensus): {WORK_DIR}/variants/combined/{sample}_combined.tsv
      * Columns: sample_id, chrom, pos, ref (=C), alt, AF, DP, tiers, etc.
      * Contains detected variants relative to consensus

Output:
    - {WORK_DIR}/variants/unified/{sample}_combined.unified.tsv
      * SNPs: Flipped to unified reference (ref=R)
      * INDELs: Pass-through from View B
      * Additional rows: Consensus-ref differences (AF=1.0)

Adjustable Parameters:
    --test SAMPLE: Test mode for single sample
    --list: List available samples
    --all: Process all samples
    -y, --yes: Skip confirmation
    --view-a-dir PATH: Override View A directory path

Flipping Logic:
    Scenario 1 (R == C): Output original variants (ref=R=C)
    Scenario 2 (R != C, variant at pos): Generate alt=C with AF=1-sum(AFs), keep alt!=R,C
    Scenario 3 (R != C, no variant): Add row ref=R, alt=C, AF=1.0

Usage:
    # Test single sample
    python 11.8_flip_to_unified_batch.py --test 10090
    
    # List available samples
    python 11.8_flip_to_unified_batch.py --list
    
    # Process all samples
    python 11.8_flip_to_unified_batch.py --all -y

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import re
import glob
import sys
import argparse
import pandas as pd
import numpy as np

# ========== PATH CONFIGURATION ==========
# Working directory
WORK_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/variants_call_10"

# View A directory (consensus vs unified reference)
DEFAULT_A_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/consensus_ref_viewa_9"

# View B directory (combined variants)
B_DIR = os.path.join(WORK_DIR, "variants", "combined")

# Output directory
OUT_DIR = os.path.join(WORK_DIR, "variants", "unified")
# ========================================

def _is_snp_row(df: pd.DataFrame) -> pd.Series:
    """Robustly identify SNP rows in View A (prefer mutat=='SNP', else length==1)"""
    if 'mutat' in df.columns:
        return df['mutat'].astype(str).str.upper().eq('SNP')
    if 'length' in df.columns:
        return df['length'].astype(str).astype(float).eq(1.0)
    # Fallback: assume all are SNPs
    return pd.Series(True, index=df.index)

def flip_viewB_to_unified(A_df: pd.DataFrame, B_df: pd.DataFrame) -> pd.DataFrame:
    """
    Flip View B (consensus as REF) to unified reference perspective.
    
    Args:
        A_df: View A with positions where R!=C (A.ref=R, A.alt=C)
        B_df: View B with variants detected in reads vs consensus
    
    Returns:
        Flipped dataframe with unified reference as REF
    """
    
    # (1) Build mapping from View A SNPs only: pos -> (R, C)
    if not A_df.empty:
        A_snp = A_df[_is_snp_row(A_df)].copy()
    else:
        A_snp = A_df.copy()
    
    a_map = {int(r.pos): (str(r.ref), str(r.alt)) for r in A_snp.itertuples()}

    # (2) Split B into SNPs and INDELs
    snps = B_df[B_df['var_type'] == 'SNP'].copy()
    indels = B_df[B_df['var_type'] != 'SNP'].copy()  # Pass-through

    out_rows = []

    # (3) Process SNPs position by position
    for pos, grp in snps.groupby('pos', sort=True):
        grp = grp.copy()
        # B.ref is consensus base C (should be consistent within position)
        C = str(grp['ref'].iloc[0])

        # Unified reference base R: from View A if available, else R=C
        if pos in a_map:
            R, _C_from_A = a_map[pos]
        else:
            R = C

        # Calculate sum(AF) and f_C = 1 - sum(AF)
        def safe_sum(series):
            x = series.dropna()
            return float(x.sum()) if len(x) else np.nan

        sum_AF_primary = safe_sum(grp['AF_primary']) if 'AF_primary' in grp else np.nan
        sum_AF_clair3 = safe_sum(grp['AF_clair3']) if 'AF_clair3' in grp else np.nan
        sum_AF_ivar = safe_sum(grp['AF_ivar']) if 'AF_ivar' in grp else np.nan

        def one_minus(x):
            if np.isnan(x):
                return np.nan
            return max(0.0, min(1.0, 1.0 - float(x)))

        fC_primary = one_minus(sum_AF_primary)
        fC_clair3 = one_minus(sum_AF_clair3)
        fC_ivar = one_minus(sum_AF_ivar)

        # Representative DP (take max for safety)
        DP_primary = float(grp['DP_primary'].max()) if 'DP_primary' in grp else np.nan
        DP_clair3 = float(grp['DP_clair3'].max()) if 'DP_clair3' in grp else np.nan
        DP_ivar = float(grp['DP_ivar'].max()) if 'DP_ivar' in grp else np.nan

        if R == C:
            # S1/S2: R==C → Output as-is, explicitly set ref=R
            for r in grp.itertuples(index=False):
                row = r._asdict()
                row['ref'] = R
                out_rows.append(row)
        else:
            # S4-S7: R!=C → Generate alt=C, drop alt=R, keep alt!=R,C
            rep = grp.iloc[0]
            base_row = {
                'sample_id': rep.sample_id,
                'chrom': rep.chrom,
                'pos': pos,
                'var_type': 'SNP',
                'struct_len': 1,
                'DP_primary': DP_primary,
                'DP_clair3': DP_clair3,
                'DP_ivar': DP_ivar,
                'file_source': 'flipped_from_consensus'
            }
            
            # (a) ALT=C row (major allele in consensus)
            rowC = base_row.copy()
            rowC.update({
                'ref': R,
                'alt': C,
                'AF_primary': fC_primary,
                'AF_clair3': fC_clair3,
                'AF_ivar': fC_ivar,
                'AD_clair3': np.nan,
                'AD_ivar': np.nan,
                'source_label': 'FLIPPED',
                'tier': 'SNP-FLIPPED',
                'decision': 'KEEP',
                'reason': 'consensus_major_as_ALT_vs_unified_ref'
            })
            out_rows.append(rowC)

            # (b) Keep alt!=R and alt!=C variants, change ref to R
            for r in grp.itertuples(index=False):
                if r.alt in (R, C):
                    continue
                row = r._asdict()
                row['ref'] = R
                out_rows.append(row)

    # (4) Handle S3: R!=C but position not in View B (add ALT=C with AF=1.0)
    b_snp_positions = set(snps['pos'].astype(int).tolist())
    a_snp_positions = set(int(p) for p in a_map.keys())

    # Get sample_id and chrom from B for filling
    sample_id = None
    chrom = None
    if not B_df.empty:
        if 'sample_id' in B_df.columns and pd.notna(B_df['sample_id']).any():
            sample_id = str(B_df['sample_id'].dropna().iloc[0])
        if 'chrom' in B_df.columns and pd.notna(B_df['chrom']).any():
            chrom = str(B_df['chrom'].dropna().iloc[0])

    missing_positions = sorted(a_snp_positions - b_snp_positions)
    for pos in missing_positions:
        R, C = a_map[pos]
        out_rows.append({
            'sample_id': sample_id if sample_id is not None else '',
            'chrom': chrom if chrom is not None else '',
            'pos': pos,
            'ref': R,
            'alt': C,
            'var_type': 'SNP',
            'struct_len': 1,
            'AF_primary': 1.0,
            'DP_primary': np.nan,
            'AF_clair3': 1.0,
            'DP_clair3': np.nan,
            'AD_clair3': np.nan,
            'AF_ivar': 1.0,
            'DP_ivar': np.nan,
            'AD_ivar': np.nan,
            'source_label': 'FLIPPED',
            'tier': 'SNP-FLIPPED',
            'decision': 'KEEP',
            'reason': 'R!=C_and_viewB_missing__AF_major_set_to_1.0',
            'file_source': 'flipped_from_consensus'
        })

    # (5) Combine SNPs + INDELs and sort
    out_snps = pd.DataFrame(out_rows, columns=B_df.columns)
    out = pd.concat([out_snps, indels], ignore_index=True)
    out = out.sort_values(['pos', 'var_type', 'alt']).reset_index(drop=True)
    return out

def main():
    parser = argparse.ArgumentParser(
        description='Flip View B variants to unified reference perspective'
    )
    parser.add_argument('--test', type=str, help='Test mode: process single sample')
    parser.add_argument('--list', action='store_true', help='List available samples')
    parser.add_argument('--all', action='store_true', help='Process all samples')
    parser.add_argument('-y', '--yes', action='store_true', help='Skip confirmation')
    parser.add_argument('--view-a-dir', type=str, default=DEFAULT_A_DIR,
                        help=f'View A directory (default: {DEFAULT_A_DIR})')
    args = parser.parse_args()

    A_DIR = args.view_a_dir
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 60)
    print("Flip to Unified Script v2.0 for variants_call_10")
    print("=" * 60)
    print(f"View A directory: {A_DIR}")
    print(f"View B directory: {B_DIR}")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 60)

    # Index View A (CSV files)
    a_files = glob.glob(os.path.join(A_DIR, "*.csv"))
    if not a_files:
        print(f"[ERROR] No CSV files found in View A directory: {A_DIR}")
        sys.exit(1)
    
    a_index = {}
    for ap in a_files:
        base = os.path.basename(ap)
        for sid in set(re.findall(r'\d+', base)):
            a_index.setdefault(sid, []).append(ap)

    # Scan View B (TSV files)
    b_files = glob.glob(os.path.join(B_DIR, "*_combined.tsv"))
    if not b_files:
        print(f"[ERROR] No *_combined.tsv files found in {B_DIR}")
        sys.exit(2)

    # List mode
    if args.list:
        print(f"Found {len(b_files)} samples with combined.tsv:")
        for i, bp in enumerate(sorted(b_files), 1):
            base = os.path.basename(bp)
            sample_id = re.findall(r'\d+', base)[0] if re.findall(r'\d+', base) else 'unknown'
            has_view_a = "✓" if sample_id in a_index else "✗"
            print(f"{i:3d}. {sample_id} (View A: {has_view_a})")
        return

    # Test mode
    if args.test:
        b_files = [f for f in b_files if args.test in os.path.basename(f)]
        if not b_files:
            print(f"[ERROR] Sample {args.test} not found in {B_DIR}")
            sys.exit(1)
        print(f"[INFO] Test mode: processing sample {args.test}")
    
    # All mode
    elif args.all:
        if not args.yes:
            response = input(f"Process {len(b_files)} samples? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Cancelled")
                return
    else:
        print("[ERROR] Please specify --test, --list, or --all")
        parser.print_help()
        return

    print(f"[INFO] Processing {len(b_files)} samples")

    summary = []
    for bp in b_files:
        try:
            B = pd.read_csv(bp, sep='\t')
        except Exception as e:
            print(f"[WARN] Skipping {bp}: read failed ({e})")
            continue

        # Extract sample_id
        if 'sample_id' not in B.columns:
            m = re.findall(r'\d+', os.path.basename(bp))
            if m:
                B['sample_id'] = m[0]
            else:
                B['sample_id'] = ''

        sample_id = str(B['sample_id'].iloc[0]) if not B.empty else (
            re.findall(r'\d+', os.path.basename(bp))[0] if re.findall(r'\d+', os.path.basename(bp)) else ''
        )

        # Find corresponding View A file
        cand = a_index.get(sample_id, [])
        if not cand:
            cand = [ap for ap in a_files if sample_id and sample_id in os.path.basename(ap)]
        
        if not cand:
            print(f"[ERROR] No View A CSV found for sample {sample_id}; file {bp}")
            continue
        
        if len(cand) > 1:
            print(f"[WARN] Multiple View A CSVs for sample {sample_id}, using first: {os.path.basename(cand[0])}")

        A = pd.read_csv(cand[0])
        out_df = flip_viewB_to_unified(A, B)

        b_base = os.path.splitext(os.path.basename(bp))[0]
        out_path = os.path.join(OUT_DIR, f"{b_base}.unified.tsv")
        out_df.to_csv(out_path, sep='\t', index=False)

        summary.append((sample_id, bp, out_path, len(B), len(out_df)))
        print(f"[OK] {sample_id}: {len(B)} rows → {len(out_df)} rows (unified)")

    if summary:
        S = pd.DataFrame(summary, columns=['sample_id', 'input_b', 'output_unified', 'n_in', 'n_out'])
        print("\n" + "=" * 60)
        print("Processing Summary:")
        print("=" * 60)
        print(S.to_string(index=False))
        
        # Save summary
        S.to_csv(os.path.join(OUT_DIR, "flip_summary.csv"), index=False)
        print(f"\nSummary saved to: {os.path.join(OUT_DIR, 'flip_summary.csv')}")
    else:
        print("[INFO] No output generated. Please check input directories and file naming.")

if __name__ == "__main__":
    main()
