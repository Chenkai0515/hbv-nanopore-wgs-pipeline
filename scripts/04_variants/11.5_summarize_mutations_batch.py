#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mutation Spectrum Summarization Script
=======================================

Purpose:
    Generate per-sample and cohort-level mutation summaries from filtered iVar and Clair3 results.
    Integrates results from both variant callers to create tiered variant lists:
    - Tier-1 (High confidence): DUAL (Clair3 ∩ iVar) ∪ Clair3 PASS ∪ iVar-only KEEP
    - Tier-2 (Exploratory): iVar-only CANDIDATE
    
    Note: Each sample uses its own consensus reference. Cohort aggregation is performed
    at the event level (chrom, pos, ref, alt) for event counting and frequency distribution.

Required Environment:
    - pandas (data manipulation)
    - pysam (optional, for faster VCF parsing)
    
    Install via conda/mamba:
        conda activate ivar  # or Clair3 (both have pandas)

Input:
    - From Step 4: {WORK_DIR}/variants/filtered/{sample_id}/
      * ivar_only_final.tsv - Filtered iVar-only variants
      * {sample_id}.ivar_mid.tsv - All iVar variants with source labels
    - From Step 3b: {WORK_DIR}/variants/clair3/{sample_id}_clair3/merge_output.vcf.gz

Output:
    - Per sample: {WORK_DIR}/variants/filtered/{sample_id}/
      * accepted.tsv - High-confidence variants (Tier-1)
      * candidates.tsv - Candidate variants (Tier-2)
    - Cohort level: {WORK_DIR}/variants/filtered/cohort_summary/
      * cohort_accepted.tsv - All accepted variants
      * cohort_candidates.tsv - All candidate variants
      * cohort_events_tier1.tsv - Event-level aggregation
      * bcp_prec_summary_tier1.tsv - BCP/PreC mutations (Tier-1)
      * bcp_prec_summary_tier2.tsv - BCP/PreC mutations (Tier-2)

Adjustable Parameters:
    --test: Test mode for single sample
    --list: List available samples
    --all: Process all samples (with confirmation)
    -y, --yes: Skip confirmation

Target Sites (BCP/PreC):
    - Position 1762: BCP A1762T mutation
    - Position 1764: BCP G1764A mutation
    - Position 1896: PreC G1896A mutation

Usage:
    # Test single sample
    python 11.5_summarize_mutations_batch.py --test 10090
    
    # List available samples
    python 11.5_summarize_mutations_batch.py --list
    
    # Process all samples
    python 11.5_summarize_mutations_batch.py --all -y

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import sys
import gzip
import glob
import argparse
from typing import List, Tuple, Optional, Dict

import pandas as pd

# ========== PATH CONFIGURATION ==========
# These paths can be overridden by:
#   1. Environment variables (HBV_WORK_DIR)
#   2. Command line arguments
#
# Default: Use relative paths from script location (for standalone use)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))

# Working directory
WORK_DIR = os.environ.get(
    "HBV_WORK_DIR",
    os.path.join(_PROJECT_DIR, "variants_call_10")
)

# Input paths (derived from WORK_DIR)
IVAR_FINAL_ROOT = os.path.join(WORK_DIR, "variants", "filtered")
CLAIR3_DIR = os.path.join(WORK_DIR, "variants", "clair3")

# Output path
OUT_DIR = os.path.join(IVAR_FINAL_ROOT, "cohort_summary")
# ========================================

GENOME_LENGTH = 3221
OFFSET = 1823  # Must match 11.6_coordinate_transform_batch.py

def ref_to_rot(pos_ref: int) -> int:
    """Convert X02763 coordinates to rotated coordinates"""
    return ((pos_ref - 1) - OFFSET) % GENOME_LENGTH + 1




# Target sites (X02763 reference coordinates)
TARGET_SITES = {
    1762: {"label": "BCP_A1762T", "expected_alt": "T"},
    1764: {"label": "BCP_G1764A", "expected_alt": "A"},
    1896: {"label": "PreC_G1896A", "expected_alt": "A"},
}

# Check for pysam
try:
    import pysam
    _HAS_PYSAM = True
except Exception:
    _HAS_PYSAM = False

def ensure_dir(p: str):
    """Create directory if it doesn't exist"""
    os.makedirs(p, exist_ok=True)

def list_samples(root: str) -> List[str]:
    """Find all samples in filtered root directory"""
    sids = []
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d):
            continue
        sid = os.path.basename(d)
        if sid == "cohort_summary":  # Skip summary directory
            continue
        if (os.path.exists(os.path.join(d, "ivar_only_final.tsv")) or
            os.path.exists(os.path.join(d, f"{sid}.ivar_mid.tsv"))):
            sids.append(sid)
    return sids

def var_type_and_len(ref: str, alt: str) -> Tuple[str, int]:
    """Determine variant type and length from VCF representation"""
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) == 1 and len(alt) == 1:
        return "SNP", 1
    if alt.startswith(ref):
        return "INS", len(alt) - len(ref)
    if ref.startswith(alt):
        return "DEL", len(ref) - len(alt)
    return "MNV", max(len(ref), len(alt))

def parse_clair3_vcf(vcf_path: str) -> pd.DataFrame:
    """Parse Clair3 VCF file and extract PASS variants"""
    rows = []
    if not os.path.exists(vcf_path):
        return pd.DataFrame(columns=[
            "chrom", "pos", "ref", "alt", "filter", "qual",
            "AF_clair3", "DP_clair3", "AD_clair3", "var_type", "struct_len"
        ])

    if _HAS_PYSAM:
        vf = pysam.VariantFile(vcf_path)
        sample_names = list(vf.header.samples)
        sm = sample_names[0] if sample_names else None
        
        for rec in vf:
            filt = "PASS" if (len(rec.filter.keys()) == 0 or "PASS" in rec.filter.keys()) else ",".join(list(rec.filter.keys()))
            if filt != "PASS":
                continue
            
            chrom = str(rec.chrom)
            pos = int(rec.pos)
            ref = rec.ref.upper()
            alts = [a.upper() for a in (rec.alts or [])]
            
            for alt in alts:
                AF = None
                DP = None
                AD = None
                
                # Try INFO fields
                if "AF" in rec.info:
                    try:
                        if isinstance(rec.info["AF"], tuple) and len(rec.info["AF"]) == len(alts):
                            AF = float(rec.info["AF"][alts.index(alt)])
                        else:
                            AF = float(rec.info["AF"])
                    except Exception:
                        pass
                
                if "DP" in rec.info:
                    try:
                        DP = int(rec.info["DP"])
                    except Exception:
                        pass
                
                # Try FORMAT fields
                if sm is not None and sm in rec.samples:
                    s = rec.samples[sm]
                    if AF is None and "AF" in s:
                        try:
                            AF = float(s["AF"])
                        except Exception:
                            pass
                    if DP is None and "DP" in s and s["DP"] is not None:
                        try:
                            DP = int(s["DP"])
                        except Exception:
                            pass
                    if "AD" in s and s["AD"] is not None:
                        try:
                            ad = list(s["AD"])
                            if len(ad) >= 2:
                                AD = int(ad[alts.index(alt) + 1]) if len(ad) > alts.index(alt) + 1 else None
                                if AF is None and DP:
                                    AF = AD / DP if DP > 0 and AD is not None else None
                        except Exception:
                            pass
                
                if AF is None and DP is not None and AD is not None and DP > 0:
                    AF = AD / DP
                
                vt, ln = var_type_and_len(ref, alt)
                rows.append({
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "filter": filt, "qual": float(rec.qual) if rec.qual is not None else None,
                    "AF_clair3": AF, "DP_clair3": DP, "AD_clair3": AD,
                    "var_type": vt, "struct_len": ln
                })
        return pd.DataFrame(rows)
    
    # Fallback: text parsing
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if len(toks) < 8:
                continue
            
            chrom, pos, _id, ref, alt, qual, filt, info = toks[:8]
            if filt != "PASS":
                continue
            
            pos = int(pos)
            ref = ref.upper()
            alts = [a.upper() for a in alt.split(",")]
            
            # Parse INFO
            info_map = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_map[k] = v
            
            # Parse FORMAT
            fmt_keys = toks[8].split(":") if len(toks) > 8 else []
            sample_fields = toks[9].split(":") if (len(toks) > 9 and toks[8]) else []
            fmt_map = {k: (sample_fields[i] if i < len(sample_fields) else None) 
                      for i, k in enumerate(fmt_keys)} if fmt_keys else {}
            
            for a_i, a in enumerate(alts):
                AF = None
                DP = None
                AD = None
                
                # Try to extract AF, DP, AD
                if "AF" in info_map:
                    try:
                        if "," in info_map["AF"]:
                            AF = float(info_map["AF"].split(",")[a_i])
                        else:
                            AF = float(info_map["AF"])
                    except Exception:
                        pass
                
                if "DP" in info_map:
                    try:
                        DP = int(info_map["DP"])
                    except Exception:
                        pass
                
                if "AF" in fmt_map and fmt_map["AF"] not in (None, "."):
                    try:
                        AF = float(fmt_map["AF"])
                    except Exception:
                        pass
                
                if "DP" in fmt_map and fmt_map["DP"] not in (None, "."):
                    try:
                        DP = int(fmt_map["DP"])
                    except Exception:
                        pass
                
                if "AD" in fmt_map and fmt_map["AD"] not in (None, "."):
                    try:
                        ad = [int(x) for x in fmt_map["AD"].split(",")]
                        if len(ad) >= 2:
                            AD = ad[a_i + 1] if len(ad) > a_i + 1 else None
                            if AF is None and DP:
                                AF = AD / DP if DP > 0 and AD is not None else None
                    except Exception:
                        pass
                
                if AF is None and DP is not None and AD is not None and DP > 0:
                    AF = AD / DP
                
                vt, ln = var_type_and_len(ref, a)
                rows.append({
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": a,
                    "filter": "PASS", "qual": float(qual) if qual not in (None, ".") else None,
                    "AF_clair3": AF, "DP_clair3": DP, "AD_clair3": AD,
                    "var_type": vt, "struct_len": ln
                })
    
    return pd.DataFrame(rows)

def load_ivar_only_final(path: str) -> pd.DataFrame:
    """Load iVar-only filtered results from Step 4"""
    if not os.path.exists(path):
        return pd.DataFrame(columns=[
            "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
            "ALT_FREQ", "TOTAL_DP", "ALT_DP", "ALT_RV", "ALT_FWD", "ALT_QUAL",
            "HmerLen", "NearPrimer", "source", "tier", "decision", "reason"
        ])
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "pos": int, "ref": str, "alt": str})
    return df

def load_ivar_mid(path: str) -> pd.DataFrame:
    """Load intermediate iVar results (contains DUAL/IV_ONLY source labels)"""
    if not os.path.exists(path):
        return pd.DataFrame(columns=[
            "chrom", "pos", "ref", "alt", "var_type", "struct_len", "ALT_FREQ", "TOTAL_DP",
            "ALT_DP", "ALT_RV", "ALT_FWD", "ALT_QUAL", "PASS", "HmerLen", "NearPrimer",
            "source", "event_key"
        ])
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "pos": int, "ref": str, "alt": str})
    
    # Ensure event_key exists
    if "event_key" not in df.columns:
        df["event_key"] = df.apply(
            lambda r: (r["chrom"], int(r["pos"]), str(r["ref"]).upper(), str(r["alt"]).upper()),
            axis=1
        )
    return df

def make_event_key(chrom: str, pos: int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    """Create standardized event key"""
    return (str(chrom), int(pos), str(ref).upper(), str(alt).upper())

def main():
    parser = argparse.ArgumentParser(
        description="Generate per-sample and cohort mutation summaries"
    )
    parser.add_argument("--test", type=str, help="Test mode: process single sample")
    parser.add_argument("--list", action="store_true", help="List all available samples")
    parser.add_argument("--all", action="store_true", help="Process all samples")
    parser.add_argument("-y", "--yes", action="store_true", help="Skip confirmation")
    args = parser.parse_args()

    ensure_dir(OUT_DIR)
    sids = list_samples(IVAR_FINAL_ROOT)
    
    if not sids:
        print(f"[ERROR] No samples found in {IVAR_FINAL_ROOT}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("Summarize Mutations Script v2.0 for variants_call_10")
    print("=" * 60)
    print(f"Input directory: {IVAR_FINAL_ROOT}")
    print(f"Clair3 directory: {CLAIR3_DIR}")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 60)

    # List mode
    if args.list:
        print(f"Found {len(sids)} samples:")
        for i, sid in enumerate(sids, 1):
            print(f"{i:3d}. {sid}")
        return

    # Test mode
    if args.test:
        if args.test not in sids:
            print(f"[ERROR] Sample {args.test} not found", file=sys.stderr)
            sys.exit(1)
        sids = [args.test]
        print(f"[INFO] Test mode: processing sample {args.test}")
    
    # All mode
    elif args.all:
        if not args.yes:
            response = input(f"Process all {len(sids)} samples? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Cancelled")
                return
    else:
        print("[ERROR] Please specify --test, --list, or --all")
        parser.print_help()
        return

    print(f"[INFO] Samples to process: {len(sids)}")

    all_accepted = []
    all_candidates = []

    for sid in sids:
        subdir = os.path.join(IVAR_FINAL_ROOT, sid)
        ivar_final = os.path.join(subdir, "ivar_only_final.tsv")
        ivar_mid = os.path.join(subdir, f"{sid}.ivar_mid.tsv")
        clair3_vcf = os.path.join(CLAIR3_DIR, f"{sid}_clair3", "merge_output.vcf.gz")

        # Load data
        df_iv = load_ivar_only_final(ivar_final)
        df_mid = load_ivar_mid(ivar_mid)
        df_c3 = parse_clair3_vcf(clair3_vcf)

        # iVar-only KEEP / CANDIDATE
        iv_keep = df_iv[df_iv["decision"] == "KEEP"].copy() if not df_iv.empty else pd.DataFrame()
        iv_cand = df_iv[df_iv["decision"] == "CANDIDATE"].copy() if not df_iv.empty else pd.DataFrame()

        # Clair3 PASS events
        c3 = df_c3.copy()
        if not c3.empty:
            c3["event_key"] = c3.apply(lambda r: make_event_key(r["chrom"], r["pos"], r["ref"], r["alt"]), axis=1)
            c3_keys = set(c3["event_key"].tolist())
        else:
            c3_keys = set()

        # All iVar events
        iv_all = df_mid.copy()
        if not iv_all.empty:
            iv_all["event_key"] = iv_all.apply(lambda r: make_event_key(r["chrom"], r["pos"], r["ref"], r["alt"]), axis=1)
            iv_all_keys = set(iv_all["event_key"].tolist())
            iv_all_idx: Dict[Tuple[str, int, str, str], Tuple[float, int, int]] = {}
            for _, r in iv_all.iterrows():
                iv_all_idx[r["event_key"]] = (
                    float(r["ALT_FREQ"]) if pd.notna(r["ALT_FREQ"]) else None,
                    int(r["TOTAL_DP"]) if pd.notna(r["TOTAL_DP"]) else None,
                    int(r["ALT_DP"]) if pd.notna(r["ALT_DP"]) else None
                )
        else:
            iv_all_keys = set()
            iv_all_idx = {}

        # iVar-KEEP event keys
        iv_keep_keys = set()
        if not iv_keep.empty:
            if "event_key" not in iv_keep.columns:
                iv_keep["event_key"] = iv_keep.apply(
                    lambda r: make_event_key(r["chrom"], r["pos"], r["ref"], r["alt"]),
                    axis=1
                )
            iv_keep_keys = set(iv_keep["event_key"].tolist())

        # Build accepted list
        # 1) Clair3 PASS
        if not c3.empty:
            c3["sample_id"] = sid
            c3["source_label"] = c3["event_key"].apply(lambda k: "DUAL" if k in iv_all_keys else "CLAIR3_ONLY")
            c3["tier"] = "CLAIR3-PASS"
            c3["decision"] = "KEEP"
            c3["reason"] = "clair3_pass"
            c3["AF_primary"] = c3["AF_clair3"]
            c3["DP_primary"] = c3["DP_clair3"]
            c3["AF_ivar"] = c3["event_key"].apply(lambda k: iv_all_idx.get(k, (None, None, None))[0])
            c3["DP_ivar"] = c3["event_key"].apply(lambda k: iv_all_idx.get(k, (None, None, None))[1])
            c3["AD_ivar"] = c3["event_key"].apply(lambda k: iv_all_idx.get(k, (None, None, None))[2])
            
            c3_accept = c3[[
                "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
                "AF_primary", "DP_primary", "AF_clair3", "DP_clair3", "AD_clair3",
                "AF_ivar", "DP_ivar", "AD_ivar", "source_label", "tier", "decision", "reason"
            ]].copy()
        else:
            c3_accept = pd.DataFrame(columns=[
                "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
                "AF_primary", "DP_primary", "AF_clair3", "DP_clair3", "AD_clair3",
                "AF_ivar", "DP_ivar", "AD_ivar", "source_label", "tier", "decision", "reason"
            ])

        # 2) iVar-only KEEP (not in Clair3 PASS)
        if not iv_keep.empty:
            iv_keep_only = iv_keep[~iv_keep["event_key"].isin(c3_keys)].copy()
            if not iv_keep_only.empty:
                iv_keep_only["sample_id"] = sid
                iv_keep_only["AF_primary"] = iv_keep_only["ALT_FREQ"]
                iv_keep_only["DP_primary"] = iv_keep_only["TOTAL_DP"]
                iv_keep_only["AF_clair3"] = None
                iv_keep_only["DP_clair3"] = None
                iv_keep_only["AD_clair3"] = None
                iv_keep_only["AF_ivar"] = iv_keep_only["ALT_FREQ"]
                iv_keep_only["DP_ivar"] = iv_keep_only["TOTAL_DP"]
                iv_keep_only["AD_ivar"] = iv_keep_only["ALT_DP"]
                iv_keep_only["source_label"] = "IVAR_ONLY_KEEP"
                
                iv_keep_accept = iv_keep_only[[
                    "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
                    "AF_primary", "DP_primary", "AF_clair3", "DP_clair3", "AD_clair3",
                    "AF_ivar", "DP_ivar", "AD_ivar", "source_label", "tier", "decision", "reason"
                ]].copy()
            else:
                iv_keep_accept = pd.DataFrame(columns=c3_accept.columns)
        else:
            iv_keep_accept = pd.DataFrame(columns=c3_accept.columns)

        accepted = pd.concat([c3_accept, iv_keep_accept], ignore_index=True)
        
        # Save per-sample accepted
        acc_path = os.path.join(IVAR_FINAL_ROOT, sid, "accepted.tsv")
        accepted.to_csv(acc_path, sep="\t", index=False)
        all_accepted.append(accepted)
        print(f"[OK] {sid}: {len(accepted)} accepted variants")

        # Build candidates list
        cand = pd.DataFrame(columns=[
            "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
            "AF", "DP", "AD", "tier", "decision", "reason", "source_label"
        ])
        if not iv_cand.empty:
            tmp = iv_cand.copy()
            tmp["sample_id"] = sid
            tmp["AF"] = tmp["ALT_FREQ"]
            tmp["DP"] = tmp["TOTAL_DP"]
            tmp["AD"] = tmp["ALT_DP"]
            tmp["source_label"] = "IVAR_ONLY_CANDIDATE"
            cand = tmp[[
                "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
                "AF", "DP", "AD", "tier", "decision", "reason", "source_label"
            ]]
        
        cand_path = os.path.join(IVAR_FINAL_ROOT, sid, "candidates.tsv")
        cand.to_csv(cand_path, sep="\t", index=False)
        all_candidates.append(cand)

    # Cohort aggregation
    cohort_acc = pd.concat(all_accepted, ignore_index=True) if all_accepted else pd.DataFrame()
    cohort_cand = pd.concat(all_candidates, ignore_index=True) if all_candidates else pd.DataFrame()

    ensure_dir(OUT_DIR)
    if not cohort_acc.empty:
        cohort_acc.to_csv(os.path.join(OUT_DIR, "cohort_accepted.tsv"), sep="\t", index=False)
    if not cohort_cand.empty:
        cohort_cand.to_csv(os.path.join(OUT_DIR, "cohort_candidates.tsv"), sep="\t", index=False)

    # Event-level aggregation (Tier-1)
    if not cohort_acc.empty:
        grp_cols = ["chrom", "pos", "ref", "alt", "var_type", "struct_len"]
        
        def q25(x):
            return x.quantile(0.25)
        
        def q75(x):
            return x.quantile(0.75)
        
        af_series = cohort_acc["AF_primary"].copy()
        af_series = af_series.fillna(cohort_acc.get("AF_ivar"))
        cohort_acc["_af_for_stat"] = af_series

        agg = (cohort_acc
               .groupby(grp_cols, as_index=False)
               .agg(n_samples=("sample_id", "nunique"),
                    n_events=("sample_id", "size"),
                    n_dual=("source_label", lambda s: (s == "DUAL").sum()),
                    n_clair3_only=("source_label", lambda s: (s == "CLAIR3_ONLY").sum()),
                    n_ivar_only_keep=("source_label", lambda s: (s == "IVAR_ONLY_KEEP").sum()),
                    af_median=("_af_for_stat", "median"),
                    af_mean=("_af_for_stat", "mean"),
                    af_q25=("_af_for_stat", q25),
                    af_q75=("_af_for_stat", q75)
                    ))
        agg.to_csv(os.path.join(OUT_DIR, "cohort_events_tier1.tsv"), sep="\t", index=False)

    # Target site summaries
    
    def summarize_sites(df: pd.DataFrame, label: str) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame(columns=[
                "site_label", "pos", "sample_id", "chrom", "ref", "alt",
                "AF", "DP", "source_label", "tier", "decision"
            ])
        
        out_rows = []
        for pos_ref, meta in TARGET_SITES.items():
            lbl = meta["label"]
            # Convert X02763 coordinates to rotated coordinates
            pos_rot = ref_to_rot(pos_ref)
            # Ensure pos column type matches
            sub = df[df["pos"] == pos_rot].copy()
            if sub.empty:
                continue
            
            # Get AF value
            if "AF_primary" in sub.columns:
                af_col = sub["AF"] if "AF" in sub.columns else pd.Series([None] * len(sub))
                af_ivar_col = sub["AF_ivar"] if "AF_ivar" in sub.columns else pd.Series([None] * len(sub))
                sub["AF_use"] = sub["AF_primary"].fillna(af_col).fillna(af_ivar_col)
            elif "AF" in sub.columns:
                sub["AF_use"] = sub["AF"]
            else:
                sub["AF_use"] = None
            
            for _, r in sub.iterrows():
                out_rows.append({
                    "site_label": lbl,
                    "pos": pos_ref,
                    "sample_id": r["sample_id"],
                    "chrom": r["chrom"],
                    "ref": r["ref"],
                    "alt": r["alt"],
                    "AF": r["AF_use"],
                    "DP": r.get("DP_primary", r.get("DP")),
                    "source_label": r.get("source_label"),
                    "tier": r.get("tier"),
                    "decision": r.get("decision"),
                })
        return pd.DataFrame(out_rows)

    if not cohort_acc.empty:
        bcp_prec_acc = summarize_sites(cohort_acc, "Tier1")
        bcp_prec_acc.to_csv(os.path.join(OUT_DIR, "bcp_prec_summary_tier1.tsv"), sep="\t", index=False)
    
    if not cohort_cand.empty:
        bcp_prec_cand = summarize_sites(cohort_cand, "Tier2")
        bcp_prec_cand.to_csv(os.path.join(OUT_DIR, "bcp_prec_summary_tier2.tsv"), sep="\t", index=False)

    print("[DONE] Summarization complete. Output directory:", OUT_DIR)
    print("  - cohort_accepted.tsv / cohort_candidates.tsv")
    print("  - cohort_events_tier1.tsv")
    print("  - bcp_prec_summary_tier1.tsv / bcp_prec_summary_tier2.tsv")
    print("  - Per-sample: accepted.tsv / candidates.tsv")

if __name__ == "__main__":
    main()
