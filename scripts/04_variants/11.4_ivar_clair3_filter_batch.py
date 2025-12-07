#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
iVar and Clair3 Joint Variant Filtering Script
===============================================

Purpose:
    Comprehensive filtering and quality assessment of variants detected by iVar and Clair3.
    Combines results from both tools to generate high-quality variant calls with tiered
    confidence levels and cohort-level blacklist filtering.

Required Environment:
    - pandas (data manipulation)
    - pysam (optional, for faster VCF parsing)
    
    Install via conda/mamba:
        conda activate ivar  # or Clair3 (both have pandas)
        pip install pysam  # optional for better performance

Input:
    - iVar results: {WORK_DIR}/variants/ivar/{sample_id}/variants.tsv
    - Clair3 results: {WORK_DIR}/variants/clair3/{sample_id}_clair3/merge_output.vcf.gz
    - Consensus sequences: {CONSENSUS_DIR}/{sample_id}/consensus.fasta

Output:
    - Per sample: {WORK_DIR}/variants/filtered/{sample_id}/
      * ivar_only_final.tsv - Filtered iVar-only variants with tier/decision
      * {sample_id}.ivar_mid.tsv - Intermediate results (all iVar variants)
      * accepted.tsv - High-confidence variants (from Step 5)
      * candidates.tsv - Candidate variants for further validation (from Step 5)
    - Cohort level: {WORK_DIR}/variants/filtered/
      * blacklist.tsv - Recurrent low-frequency artifacts
      * summary_counts.tsv - Per-sample variant counts
      * reason_breakdown.tsv - Rejection reason statistics

Adjustable Parameters:
    --profile: Filter profile (balanced/strict/lenient, default: lenient)
    --test: Test mode for single sample
    --list: List available samples
    --all: Process all samples (with confirmation)
    -y, --yes: Skip confirmation
    
Filter Profiles:
    - strict: High stringency, maximum specificity
    - balanced: Recommended, good sensitivity/specificity balance
    - lenient: Lower thresholds, maximum sensitivity (default)

Quality Tiers:
    - HC (High Confidence): AF ≥20%, high coverage, strict QC
    - MC (Medium Confidence): 5% ≤ AF <20%, good coverage, standard QC
    - LF (Low Frequency): 1% ≤ AF <5%, high coverage, candidate status
    
Variant Sources:
    - DUAL: Detected by both iVar and Clair3 (highest confidence)
    - IV_ONLY: Detected by iVar only (subject to filtering)
    - CLAIR3_ONLY: Detected by Clair3 only

Usage:
    # Test single sample
    python 11.4_ivar_clair3_filter_batch.py --test 10090
    
    # List available samples
    python 11.4_ivar_clair3_filter_batch.py --list
    
    # Process all samples
    python 11.4_ivar_clair3_filter_batch.py --all
    
    # Use strict profile
    python 11.4_ivar_clair3_filter_batch.py --all --profile strict -y

Version: 2.0 (Adapted for variants_call_10 pipeline)
"""

import os
import sys
import gzip
import glob
import math
import re
import argparse
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Set

import pandas as pd

# ========== PATH CONFIGURATION ==========
# These paths can be overridden by:
#   1. Environment variables (HBV_WORK_DIR, HBV_CONSENSUS_DIR)
#   2. Command line arguments (--work-dir, --consensus-dir)
#
# Default: Use relative paths from script location (for standalone use)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))

# Working directory (variants_call_10)
WORK_DIR = os.environ.get(
    "HBV_WORK_DIR",
    os.path.join(_PROJECT_DIR, "variants_call_10")
)

# Consensus sequences directory (Medaka_consensus_8/r2)
CONSENSUS_DIR = os.environ.get(
    "HBV_CONSENSUS_DIR",
    os.path.join(_PROJECT_DIR, "Medaka_consensus_8", "r2")
)

# Input paths (derived from WORK_DIR)
IVAR_DIR = os.path.join(WORK_DIR, "variants", "ivar")
CLAIR3_DIR = os.path.join(WORK_DIR, "variants", "clair3")

# Output path (derived from WORK_DIR)
OUT_ROOT = os.path.join(WORK_DIR, "variants", "filtered")
# ========================================

# Filter profile (can be overridden by command line)
DEFAULT_FILTER_PROFILE = "lenient"

# Primer region (unified coordinates) - set to 0 to disable primer filtering
PRIMER_START, PRIMER_END, PRIMER_FLANK = 0, 0, 0

# Blacklist (PoN) construction thresholds
BLACKLIST_AF_MAX = 0.02          # Max AF for "low-frequency recurrent" (<=2%)
BLACKLIST_MIN_FRAC = 0.10        # Appears in ≥10% of samples
BLACKLIST_MIN_COUNT_FLOOR = 10   # Or at least 10 samples (whichever is larger)

# Blacklist application rules
BLACKLIST_FORCE_AF = 0.05   # Force reject if AF < 5%
BLACKLIST_CAND_AF = 0.10    # Mark as candidate if 5% ≤ AF < 10%

# Filter profile definitions
PROFILES = {
    "strict": {
        "AF_HIGH": 0.20, "AF_MID_LOW": 0.05, "AF_LOW": 0.01,
        # SNP
        "SNP_HC_DP_MIN": 200, "SNP_HC_AD_MIN": 50,
        "SNP_MC_DP_MIN": 1000, "SNP_MC_AD_MIN": 20,
        "SNP_LF_DP_MIN": 3000, "SNP_LF_AD_MIN": 30,
        # INDEL
        "INDEL_HC_AF_MIN": 0.10, "INDEL_HC_DP_MIN": 500, "INDEL_HC_AD_MIN": 50,
        "INDEL_MC_DP_MIN": 1000, "INDEL_MC_AD_MIN": 25,
        "INDEL_LF_DP_MIN": 3000, "INDEL_LF_AD_MIN": 40,
        # Quality/strand/homopolymer
        "ALT_QUAL_MIN": 25,
        "STRAND_MIN_READS": 5, "STRAND_MIN_FRAC": 0.20,
        "HMER_FAIL_LEN_SNP": 6, "HMER_FAIL_LEN_INDEL": 6,
        # Structure length
        "INS_MAX_HC": 10, "DEL_MAX_HC": 33, "INS_MAX_LF": 5, "DEL_MAX_LF": 15,
        # Window clustering
        "CLUSTER_WIN": 10, "CLUSTER_COUNT": 5,
    },
    "balanced": {
        "AF_HIGH": 0.20, "AF_MID_LOW": 0.05, "AF_LOW": 0.01,
        # SNP (lowered MC/LF coverage thresholds; ALT_QUAL=20 matches iVar testing)
        "SNP_HC_DP_MIN": 200, "SNP_HC_AD_MIN": 30,
        "SNP_MC_DP_MIN": 500, "SNP_MC_AD_MIN": 15,
        "SNP_LF_DP_MIN": 2000, "SNP_LF_AD_MIN": 25,
        # INDEL (still strict, but slightly relaxed from strict)
        "INDEL_HC_AF_MIN": 0.12, "INDEL_HC_DP_MIN": 500, "INDEL_HC_AD_MIN": 40,
        "INDEL_MC_DP_MIN": 800, "INDEL_MC_AD_MIN": 25,
        "INDEL_LF_DP_MIN": 2500, "INDEL_LF_AD_MIN": 35,
        "ALT_QUAL_MIN": 20,
        # Strand strategy: only enforce when "evaluable" (ALT_FWD+ALT_REV ≥ 10)
        "STRAND_MIN_READS": 5, "STRAND_MIN_FRAC": 0.20,
        "HMER_FAIL_LEN_SNP": 6, "HMER_FAIL_LEN_INDEL": 6,
        "INS_MAX_HC": 10, "DEL_MAX_HC": 33, "INS_MAX_LF": 5, "DEL_MAX_LF": 15,
        "CLUSTER_WIN": 10, "CLUSTER_COUNT": 5,
    },
    "lenient": {
        "AF_HIGH": 0.20, "AF_MID_LOW": 0.05, "AF_LOW": 0.01,
        "SNP_HC_DP_MIN": 150, "SNP_HC_AD_MIN": 20,
        "SNP_MC_DP_MIN": 400, "SNP_MC_AD_MIN": 10,
        "SNP_LF_DP_MIN": 1500, "SNP_LF_AD_MIN": 20,
        "INDEL_HC_AF_MIN": 0.10, "INDEL_HC_DP_MIN": 400, "INDEL_HC_AD_MIN": 30,
        "INDEL_MC_DP_MIN": 600, "INDEL_MC_AD_MIN": 20,
        "INDEL_LF_DP_MIN": 2000, "INDEL_LF_AD_MIN": 30,
        "ALT_QUAL_MIN": 20,
        "STRAND_MIN_READS": 4, "STRAND_MIN_FRAC": 0.15,
        "HMER_FAIL_LEN_SNP": 6, "HMER_FAIL_LEN_INDEL": 6,
        "INS_MAX_HC": 12, "DEL_MAX_HC": 40, "INS_MAX_LF": 6, "DEL_MAX_LF": 20,
        "CLUSTER_WIN": 12, "CLUSTER_COUNT": 6,
    }
}

# Will be set by command line or default
P = None

# Check for pysam
try:
    import pysam
    _HAS_PYSAM = True
except Exception:
    _HAS_PYSAM = False

def ensure_dir(p: str):
    """Create directory if it doesn't exist"""
    os.makedirs(p, exist_ok=True)

def list_samples_from_ivar(ivar_dir: str) -> List[str]:
    """Find all samples with iVar variants.tsv"""
    out = []
    for d in sorted(glob.glob(os.path.join(ivar_dir, "*"))):
        if os.path.isdir(d) and os.path.exists(os.path.join(d, "variants.tsv")):
            out.append(os.path.basename(d))
    return out

def load_consensus_seq(fa: str) -> Optional[str]:
    """Load consensus sequence from FASTA file"""
    if not os.path.exists(fa):
        return None
    seq, seen = [], False
    with open(fa, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if seen:
                    break
                seen = True
                continue
            seq.append(line.strip())
    return "".join(seq).upper() if seq else None

def hmer_len_at(seq: str, pos1: int) -> int:
    """Calculate homopolymer length at given position (1-based)"""
    if not seq:
        return 0
    i = pos1 - 1
    if i < 0 or i >= len(seq):
        return 0
    b = seq[i]
    if b not in "ACGT":
        return 0
    L = i
    while L - 1 >= 0 and seq[L - 1] == b:
        L -= 1
    R = i
    while R + 1 < len(seq) and seq[R + 1] == b:
        R += 1
    return R - L + 1

def is_near_primer(pos1: int) -> bool:
    """Check if position is near primer region"""
    return (PRIMER_START - PRIMER_FLANK) <= pos1 <= (PRIMER_END + PRIMER_FLANK)

def parse_clair3_events(vcf_path: str) -> Set[Tuple[str, int, str, str]]:
    """Parse Clair3 VCF to extract variant events"""
    ev = set()
    if not os.path.exists(vcf_path):
        return ev
    
    if _HAS_PYSAM:
        vf = pysam.VariantFile(vcf_path)
        for rec in vf:
            chrom = str(rec.chrom)
            pos = int(rec.pos)
            ref = rec.ref.upper()
            for alt in rec.alts or []:
                ev.add((chrom, pos, ref, alt.upper()))
        return ev
    
    # Fallback: text parsing
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = toks[0], int(toks[1]), toks[3], toks[4]
            for a in alt.split(","):
                ev.add((chrom, pos, ref.upper(), a.upper()))
    return ev

def ivar_to_vcf_rep(row: pd.Series) -> Tuple[str, int, str, str, str, int]:
    """Convert iVar representation to VCF-style representation"""
    chrom = str(row["REGION"])
    pos = int(row["POS"])
    ref = str(row["REF"]).upper()
    alt = str(row["ALT"]).upper()
    
    if re.fullmatch(r"[ACGT]", alt):
        return chrom, pos, ref, alt, "SNP", 1
    if alt.startswith("+"):
        ins = alt[1:]
        return chrom, pos, ref, ref + ins, "INS", len(ins)
    if alt.startswith("-"):
        dele = alt[1:]
        return chrom, pos, ref + dele, ref, "DEL", len(dele)
    return chrom, pos, ref, alt, "UNK", 0

def load_ivar_tsv(tsv: str) -> pd.DataFrame:
    """Load and parse iVar variants.tsv file"""
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    need = ["REGION", "POS", "REF", "ALT", "REF_DP", "REF_RV", "REF_QUAL",
            "ALT_DP", "ALT_RV", "ALT_QUAL", "ALT_FREQ", "TOTAL_DP", "PVAL", "PASS"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise RuntimeError(f"Missing columns {miss} in {tsv}")
    
    # Convert numeric columns
    num = ["POS", "REF_DP", "REF_RV", "REF_QUAL", "ALT_DP", "ALT_RV", "ALT_QUAL", "ALT_FREQ", "TOTAL_DP", "PVAL"]
    for c in num:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    
    df["PASS"] = df["PASS"].astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    df["ALT_FWD"] = df["ALT_DP"] - df["ALT_RV"]
    df["REF_FWD"] = df["REF_DP"] - df["REF_RV"]
    return df

def strand_eval_status(alt_fwd: float, alt_rev: float) -> str:
    """Evaluate strand bias status"""
    total = (alt_fwd or 0) + (alt_rev or 0)
    if total < 2 * P["STRAND_MIN_READS"]:
        return "unevaluable"
    if (alt_fwd >= P["STRAND_MIN_READS"] and alt_rev >= P["STRAND_MIN_READS"] and
        min(alt_fwd, alt_rev) / max(1, (alt_fwd + alt_rev)) >= P["STRAND_MIN_FRAC"]):
        return "ok"
    return "fail"

@dataclass
class TierResult:
    tier: str
    decision: str
    reasons: List[str]

def classify_row(r: pd.Series) -> TierResult:
    """Classify a variant row into tier and decision"""
    vt = r["var_type"]
    AF = float(r["ALT_FREQ"])
    DP = int(r["TOTAL_DP"])
    AD = int(r["ALT_DP"])
    ALTQ = float(r["ALT_QUAL"])
    HMER = int(r.get("HmerLen", 0))
    NEAR = bool(r.get("NearPrimer", False))
    ALT_FWD = float(r.get("ALT_FWD", 0))
    ALT_REV = float(r.get("ALT_RV", 0))
    LEN = int(r.get("struct_len", 0))
    reasons = []

    # Initial filters
    if not r["PASS"]:
        return TierResult("REJECT", "REJECT", ["ivar_PASS_FALSE"])
    if ALTQ < P["ALT_QUAL_MIN"]:
        return TierResult("REJECT", "REJECT", ["low_alt_qual"])
    if NEAR:
        return TierResult("REJECT", "REJECT", ["near_primer"])

    strand_status = strand_eval_status(ALT_FWD, ALT_REV)
    
    def strand_guard(ad_needed: int = 0, af_needed: float = 0.0) -> Optional[str]:
        """Check strand bias with fallback to AD/AF coverage"""
        if strand_status == "ok":
            return None
        if strand_status == "fail":
            if AD >= ad_needed and AF >= af_needed:
                return "strand_fail_but_covered"
            return "strand_fail"
        return "strand_unevaluable"

    # SNP classification
    if vt == "SNP":
        # High confidence
        sg = strand_guard(ad_needed=100, af_needed=0.25)
        if (AF >= P["AF_HIGH"] and DP >= P["SNP_HC_DP_MIN"] and 
            AD >= max(P["SNP_HC_AD_MIN"], int(0.15 * DP)) and HMER < P["HMER_FAIL_LEN_SNP"] and
            (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("SNP-HC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Medium confidence
        sg = strand_guard(ad_needed=50, af_needed=0.08)
        if (P["AF_MID_LOW"] <= AF < P["AF_HIGH"] and DP >= P["SNP_MC_DP_MIN"] and
            AD >= max(P["SNP_MC_AD_MIN"], int(0.03 * DP)) and HMER < P["HMER_FAIL_LEN_SNP"] and
            (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("SNP-MC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Low frequency
        sg = strand_guard(ad_needed=P["SNP_LF_AD_MIN"], af_needed=0.015)
        if (P["AF_LOW"] <= AF < P["AF_MID_LOW"] and DP >= P["SNP_LF_DP_MIN"] and
            AD >= P["SNP_LF_AD_MIN"] and HMER <= 4 and
            (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("SNP-LF", "CANDIDATE", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Rejected
        why = []
        if AF < P["AF_LOW"]: why.append("af_low")
        if DP < min(P["SNP_LF_DP_MIN"], P["SNP_MC_DP_MIN"]): why.append("dp_low")
        if AD < min(P["SNP_LF_AD_MIN"], P["SNP_MC_AD_MIN"]): why.append("ad_low")
        if HMER >= P["HMER_FAIL_LEN_SNP"]: why.append("homopolymer")
        if strand_status == "fail": why.append("strand_fail")
        return TierResult("REJECT", "REJECT", why or ["snp_thresholds"])

    # INDEL classification
    if vt == "INS":
        if LEN > 50:
            return TierResult("REJECT", "REJECT", ["long_insertion_gt50"])
        sg = strand_guard(ad_needed=150, af_needed=0.20)
        
        # High confidence
        if (AF >= P["INDEL_HC_AF_MIN"] and DP >= P["INDEL_HC_DP_MIN"] and
            AD >= max(P["INDEL_HC_AD_MIN"], int(0.08 * DP)) and
            (HMER < P["HMER_FAIL_LEN_INDEL"] or (AF >= 0.30 and DP >= 1000)) and
            LEN <= P["INS_MAX_HC"] and (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-HC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Medium confidence
        sg = strand_guard(ad_needed=60, af_needed=0.08)
        if (P["AF_MID_LOW"] <= AF < P["INDEL_HC_AF_MIN"] and DP >= P["INDEL_MC_DP_MIN"] and
            AD >= P["INDEL_MC_AD_MIN"] and HMER < P["HMER_FAIL_LEN_INDEL"] and
            LEN <= P["INS_MAX_HC"] and (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-MC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Low frequency
        sg = strand_guard(ad_needed=P["INDEL_LF_AD_MIN"], af_needed=0.02)
        if (P["AF_LOW"] <= AF < P["AF_MID_LOW"] and DP >= P["INDEL_LF_DP_MIN"] and
            AD >= P["INDEL_LF_AD_MIN"] and HMER <= 4 and LEN <= P["INS_MAX_LF"] and
            (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-LF", "CANDIDATE", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Rejected
        why = []
        if AF < P["AF_LOW"]: why.append("af_low")
        if DP < min(P["INDEL_LF_DP_MIN"], P["INDEL_MC_DP_MIN"]): why.append("dp_low")
        if AD < min(P["INDEL_LF_AD_MIN"], P["INDEL_MC_AD_MIN"]): why.append("ad_low")
        if HMER >= P["HMER_FAIL_LEN_INDEL"]: why.append("homopolymer")
        if LEN > P["INS_MAX_HC"]: why.append("len_too_long")
        if strand_status == "fail": why.append("strand_fail")
        return TierResult("REJECT", "REJECT", why or ["ins_thresholds"])

    if vt == "DEL":
        if LEN > 100:
            return TierResult("REJECT", "REJECT", ["long_deletion_gt100"])
        sg = strand_guard(ad_needed=150, af_needed=0.20)
        
        # High confidence
        if (AF >= P["INDEL_HC_AF_MIN"] and DP >= P["INDEL_HC_DP_MIN"] and
            AD >= max(P["INDEL_HC_AD_MIN"], int(0.08 * DP)) and
            (HMER < P["HMER_FAIL_LEN_INDEL"] or (AF >= 0.30 and DP >= 1000)) and
            LEN <= P["DEL_MAX_HC"] and (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-HC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Medium confidence
        sg = strand_guard(ad_needed=60, af_needed=0.08)
        if (P["AF_MID_LOW"] <= AF < P["INDEL_HC_AF_MIN"] and DP >= P["INDEL_MC_DP_MIN"] and
            AD >= P["INDEL_MC_AD_MIN"] and HMER < P["HMER_FAIL_LEN_INDEL"] and
            LEN <= P["DEL_MAX_HC"] and (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-MC", "KEEP", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Low frequency
        sg = strand_guard(ad_needed=P["INDEL_LF_AD_MIN"], af_needed=0.02)
        if (P["AF_LOW"] <= AF < P["AF_MID_LOW"] and DP >= P["INDEL_LF_DP_MIN"] and
            AD >= P["INDEL_LF_AD_MIN"] and HMER <= 4 and LEN <= P["DEL_MAX_LF"] and
            (sg in (None, "strand_unevaluable", "strand_fail_but_covered"))):
            return TierResult("INDEL-LF", "CANDIDATE", [] if sg in (None, "strand_unevaluable") else [sg])
        
        # Rejected
        why = []
        if AF < P["AF_LOW"]: why.append("af_low")
        if DP < min(P["INDEL_LF_DP_MIN"], P["INDEL_MC_DP_MIN"]): why.append("dp_low")
        if AD < min(P["INDEL_LF_AD_MIN"], P["INDEL_MC_AD_MIN"]): why.append("ad_low")
        if HMER >= P["HMER_FAIL_LEN_INDEL"]: why.append("homopolymer")
        if LEN > P["DEL_MAX_HC"]: why.append("len_too_long")
        if strand_status == "fail": why.append("strand_fail")
        return TierResult("REJECT", "REJECT", why or ["del_thresholds"])

    return TierResult("REJECT", "REJECT", ["unknown_type"])

def mirror_indel_resolution(df: pd.DataFrame) -> pd.DataFrame:
    """Resolve mirror indels in homopolymer regions"""
    if df.empty:
        return df
    sub = df[(df["var_type"].isin(["INS", "DEL"]))]
    if sub.empty:
        return df
    
    to_reject = set()
    for (chrom, pos), g in sub.groupby(["chrom", "pos"]):
        if g["HmerLen"].max() < 6:
            continue
        ins = g[(g["var_type"] == "INS") & (g["struct_len"] == 1)]
        dele = g[(g["var_type"] == "DEL") & (g["struct_len"] == 1)]
        if ins.empty or dele.empty:
            continue
        
        ins_bases = {idx: str(r["alt"])[-1] for idx, r in ins.iterrows()}
        del_bases = {idx: str(r["ref"])[-1] for idx, r in dele.iterrows()}
        
        for bi in set(ins_bases.values()):
            del_idxs = [idx for idx, b in del_bases.items() if b == bi]
            ins_idxs = [idx for idx, b in ins_bases.items() if b == bi]
            if not del_idxs or not ins_idxs:
                continue
            
            # Take highest AF as primary
            best_ins = max(ins_idxs, key=lambda i: float(df.loc[i, "ALT_FREQ"]))
            best_del = max(del_idxs, key=lambda i: float(df.loc[i, "ALT_FREQ"]))
            
            if float(df.loc[best_ins, "ALT_FREQ"]) >= float(df.loc[best_del, "ALT_FREQ"]):
                primary = best_ins
                secondary = del_idxs
            else:
                primary = best_del
                secondary = ins_idxs
            
            if df.loc[primary, "decision"] != "KEEP":
                continue
            
            for j in secondary:
                AF2 = float(df.loc[j, "ALT_FREQ"])
                alt_fwd = float(df.loc[j, "ALT_FWD"])
                alt_rev = float(df.loc[j, "ALT_RV"])
                strand_status = strand_eval_status(alt_fwd, alt_rev)
                
                if AF2 <= 0.03 or strand_status == "fail" or float(df.loc[j, "ALT_QUAL"]) < P["ALT_QUAL_MIN"]:
                    to_reject.add(j)
    
    if to_reject:
        df.loc[list(to_reject), "tier"] = df.loc[list(to_reject), "tier"].astype(str) + "|MIRROR"
        df.loc[list(to_reject), "decision"] = "REJECT"
        df.loc[list(to_reject), "reason"] = df.loc[list(to_reject), "reason"].astype(str) + ";mirror_indel_in_homopolymer"
    
    return df

def cluster_suspect_resolution(df: pd.DataFrame) -> pd.DataFrame:
    """Identify and reject clustered low-frequency indels"""
    if df.empty:
        return df
    sub = df[(df["var_type"].isin(["INS", "DEL"])) & (df["ALT_FREQ"] < P["AF_MID_LOW"])].copy()
    if sub.empty:
        return df
    
    for chrom, g in sub.groupby("chrom"):
        pos = sorted(g["pos"].tolist())
        i = j = 0
        while i < len(pos):
            while j < len(pos) and pos[j] - pos[i] <= P["CLUSTER_WIN"]:
                j += 1
            if j - i >= P["CLUSTER_COUNT"]:
                s, e = pos[i], pos[j - 1]
                m = (df["chrom"].eq(chrom) & df["pos"].between(s, e) & 
                     df["var_type"].isin(["INS", "DEL"]) & (df["ALT_FREQ"] < P["AF_MID_LOW"]))
                df.loc[m, "tier"] = df.loc[m, "tier"].astype(str) + "|CLUSTER"
                df.loc[m, "decision"] = "REJECT"
                df.loc[m, "reason"] = df.loc[m, "reason"].astype(str) + ";cluster_suspect"
            i += 1
    
    return df

def process_sample(sample_id: str,
                   ivar_path: str,
                   c3_vcf_path: str,
                   consensus_fa: str,
                   out_dir: str,
                   black_counter: Dict[Tuple[str, int, str, str], Set[str]],
                   dual_counts: Dict[str, int]) -> pd.DataFrame:
    """Process a single sample"""
    os.makedirs(out_dir, exist_ok=True)
    
    df = load_ivar_tsv(ivar_path)
    rep = df.apply(ivar_to_vcf_rep, axis=1, result_type="expand")
    rep.columns = ["chrom", "pos", "vcf_ref", "vcf_alt", "var_type", "struct_len"]
    df = pd.concat([df, rep], axis=1)
    df = df[df["var_type"].isin(["SNP", "INS", "DEL"])].copy()
    df["NearPrimer"] = df["POS"].apply(lambda p: is_near_primer(int(p)))
    
    # Calculate homopolymer length
    seq = load_consensus_seq(consensus_fa) if os.path.exists(consensus_fa) else None
    df["HmerLen"] = df.apply(lambda r: hmer_len_at(seq, int(r["POS"])) if seq else 0, axis=1)
    
    # Standard keys
    df["ref"] = df["vcf_ref"].astype(str).str.upper()
    df["alt"] = df["vcf_alt"].astype(str).str.upper()
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["POS"].astype(int)
    df["event_key"] = df.apply(lambda r: (r["chrom"], int(r["pos"]), r["ref"], r["alt"]), axis=1)
    
    # Clair3 intersection
    c3 = parse_clair3_events(c3_vcf_path)
    df["source"] = df["event_key"].apply(lambda k: "DUAL" if k in c3 else "IV_ONLY")
    dual_counts[sample_id] = int((df["source"] == "DUAL").sum())

    # Collect blacklist: only iVar-only & AF<=2%
    iv_low = df[(df["source"] == "IV_ONLY") & (df["ALT_FREQ"] <= BLACKLIST_AF_MAX)]
    for k in iv_low["event_key"].unique():
        black_counter[k].add(sample_id)

    # Intermediate table
    mid_cols = ["chrom", "pos", "ref", "alt", "var_type", "struct_len", "ALT_FREQ", "TOTAL_DP",
                "ALT_DP", "ALT_RV", "ALT_FWD", "ALT_QUAL", "PASS", "HmerLen", "NearPrimer", "source", "event_key"]
    df_mid = df[mid_cols].copy()
    df_mid.to_csv(os.path.join(out_dir, f"{sample_id}.ivar_mid.tsv"), sep="\t", index=False)

    return df_mid[df_mid["source"] == "IV_ONLY"].copy()

def classify_df(df_iv: pd.DataFrame) -> pd.DataFrame:
    """Classify all variants in dataframe"""
    if df_iv.empty:
        for c in ["tier", "decision", "reason"]:
            df_iv[c] = []
        return df_iv
    
    tiers = []
    decs = []
    reasons = []
    for _, r in df_iv.iterrows():
        tr = classify_row(r)
        tiers.append(tr.tier)
        decs.append(tr.decision)
        reasons.append(";".join(tr.reasons) if tr.reasons else "")
    
    df_iv["tier"] = tiers
    df_iv["decision"] = decs
    df_iv["reason"] = reasons
    return df_iv

def apply_blacklist(df: pd.DataFrame, bl: set) -> pd.DataFrame:
    """Apply blacklist with flexible rules based on AF"""
    if df.empty or not bl:
        return df
    
    hit = df["event_key"].isin(bl)
    if not hit.any():
        return df
    
    # AF < 5% → force reject; 5% ≤ AF < 10% → downgrade to candidate; AF ≥ 10% → override
    mask_low = hit & (df["ALT_FREQ"] < BLACKLIST_FORCE_AF)
    mask_mid = hit & (df["ALT_FREQ"].between(BLACKLIST_FORCE_AF, BLACKLIST_CAND_AF, inclusive="left"))
    
    df.loc[mask_low, "decision"] = "REJECT"
    df.loc[mask_low, "tier"] = df.loc[mask_low, "tier"].astype(str) + "|BL"
    df.loc[mask_low, "reason"] = df.loc[mask_low, "reason"].astype(str) + ";blacklisted"

    df.loc[mask_mid & (df["decision"] == "KEEP"), "decision"] = "CANDIDATE"
    df.loc[mask_mid, "tier"] = df.loc[mask_mid, "tier"].astype(str) + "|BL_SOFT"
    df.loc[mask_mid, "reason"] = df.loc[mask_mid, "reason"].astype(str) + ";blacklist_soft"
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description="iVar and Clair3 joint variant filtering with tiered classification"
    )
    parser.add_argument("--test", type=str, help="Test mode: process single sample")
    parser.add_argument("--list", action="store_true", help="List all available samples")
    parser.add_argument("--all", action="store_true", help="Process all samples")
    parser.add_argument("--profile", type=str, default=DEFAULT_FILTER_PROFILE,
                        choices=["balanced", "strict", "lenient"],
                        help="Filter profile (default: lenient)")
    parser.add_argument("--use-blacklist", action="store_true", 
                        help="Apply blacklist (PoN) filtering (default: disabled)")
    parser.add_argument("-y", "--yes", action="store_true", help="Skip confirmation")
    args = parser.parse_args()

    # Set profile
    global P
    P = PROFILES[args.profile]

    os.makedirs(OUT_ROOT, exist_ok=True)
    samples = list_samples_from_ivar(IVAR_DIR)
    
    if not samples:
        print(f"[ERROR] No samples with variants.tsv found in {IVAR_DIR}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("iVar/Clair3 Filter Script v2.0 for variants_call_10")
    print("=" * 60)
    print(f"iVar directory: {IVAR_DIR}")
    print(f"Clair3 directory: {CLAIR3_DIR}")
    print(f"Consensus directory: {CONSENSUS_DIR}")
    print(f"Output directory: {OUT_ROOT}")
    print("=" * 60)

    # List mode
    if args.list:
        print(f"Found {len(samples)} samples with iVar results:")
        for i, sid in enumerate(samples, 1):
            print(f"{i:3d}. {sid}")
        return

    # Test mode
    if args.test:
        if args.test not in samples:
            print(f"[ERROR] Sample {args.test} not found in available samples", file=sys.stderr)
            sys.exit(1)
        samples = [args.test]
        print(f"[INFO] Test mode: processing sample {args.test}")
    
    # All mode
    elif args.all:
        if not args.yes:
            response = input(f"Process all {len(samples)} samples with profile '{args.profile}'? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Cancelled")
                return
    else:
        print("[ERROR] Please specify --test, --list, or --all")
        parser.print_help()
        return

    print(f"[INFO] Filter profile: {args.profile}")
    print(f"[INFO] Samples to process: {len(samples)}")

    # First pass: preprocessing + PoN statistics
    black_counter: Dict[Tuple[str, int, str, str], Set[str]] = defaultdict(set)
    dual_counts: Dict[str, int] = {}
    
    for sid in samples:
        ivar_path = os.path.join(IVAR_DIR, sid, "variants.tsv")
        c3_vcf = os.path.join(CLAIR3_DIR, f"{sid}_clair3", "merge_output.vcf.gz")
        cons_fa = os.path.join(CONSENSUS_DIR, sid, "consensus.fasta")
        out_dir = os.path.join(OUT_ROOT, sid)
        
        try:
            _ = process_sample(sid, ivar_path, c3_vcf, cons_fa, out_dir, black_counter, dual_counts)
            print(f"[OK] {sid}: Preprocessed")
        except Exception as e:
            print(f"[WARN] {sid} processing failed: {e}", file=sys.stderr)
            continue

    # Build blacklist
    n = len(samples)
    thr = max(int(math.ceil(n * BLACKLIST_MIN_FRAC)), BLACKLIST_MIN_COUNT_FLOOR)
    bl_rows = []
    for k, sset in black_counter.items():
        if len(sset) >= thr:
            chrom, pos, ref, alt = k
            bl_rows.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "n_samples_lowAF": len(sset), "frac": round(len(sset) / n, 4)
            })
    
    bl_rows = sorted(bl_rows, key=lambda x: x["n_samples_lowAF"], reverse=True)
    bl_set = set([(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in bl_rows])
    pd.DataFrame(bl_rows).to_csv(os.path.join(OUT_ROOT, "blacklist.tsv"), sep="\t", index=False)
    print(f"[INFO] Blacklist events: {len(bl_rows)} (AF≤{BLACKLIST_AF_MAX}, ≥{BLACKLIST_MIN_FRAC*100:.0f}% samples or ≥{BLACKLIST_MIN_COUNT_FLOOR} samples)")
    print(f"[INFO] Blacklist application: {'ENABLED' if args.use_blacklist else 'DISABLED'}")

    # Second pass: tier classification + conflict resolution + clustering + blacklist application
    summary = []
    reason_counter = Counter()
    
    for sid in samples:
        mid = os.path.join(OUT_ROOT, sid, f"{sid}.ivar_mid.tsv")
        if not os.path.exists(mid):
            summary.append({"sample_id": sid, "ivar_only_total": 0, "keep": 0, "candidate": 0, "reject": 0, "dual": dual_counts.get(sid, 0)})
            continue
        
        df_mid = pd.read_csv(mid, sep="\t", dtype={"chrom": str, "pos": int})
        df_iv = df_mid[df_mid["source"] == "IV_ONLY"].copy()
        
        if df_iv.empty:
            out = os.path.join(OUT_ROOT, sid, "ivar_only_final.tsv")
            pd.DataFrame(columns=[
                "sample_id", "chrom", "pos", "ref", "alt", "var_type", "struct_len",
                "ALT_FREQ", "TOTAL_DP", "ALT_DP", "ALT_RV", "ALT_FWD", "ALT_QUAL",
                "HmerLen", "NearPrimer", "source", "tier", "decision", "reason"
            ]).to_csv(out, sep="\t", index=False)
            summary.append({"sample_id": sid, "ivar_only_total": 0, "keep": 0, "candidate": 0, "reject": 0, "dual": dual_counts.get(sid, 0)})
            continue

        # Classify
        df_iv = classify_df(df_iv)
        # Mirror indel resolution
        df_iv = mirror_indel_resolution(df_iv)
        # Cluster resolution
        df_iv = cluster_suspect_resolution(df_iv)
        # Apply blacklist (only if --use-blacklist is enabled)
        if args.use_blacklist:
            df_iv = apply_blacklist(df_iv, bl_set)

        df_iv.insert(0, "sample_id", sid)
        out = os.path.join(OUT_ROOT, sid, "ivar_only_final.tsv")
        df_iv.to_csv(out, sep="\t", index=False)

        keep = int((df_iv["decision"] == "KEEP").sum())
        cand = int((df_iv["decision"] == "CANDIDATE").sum())
        rej = int((df_iv["decision"] == "REJECT").sum())
        summary.append({
            "sample_id": sid,
            "ivar_only_total": int(df_iv.shape[0]),
            "keep": keep,
            "candidate": cand,
            "reject": rej,
            "dual": dual_counts.get(sid, 0)
        })

        # Rejection reason statistics
        for rs in df_iv.loc[df_iv["decision"] == "REJECT", "reason"].fillna("NA"):
            for token in str(rs).split(";"):
                token = token.strip()
                if not token:
                    continue
                reason_counter[token] += 1

    pd.DataFrame(summary).to_csv(os.path.join(OUT_ROOT, "summary_counts.tsv"), sep="\t", index=False)
    
    # Reason statistics
    rc = pd.DataFrame(sorted(reason_counter.items(), key=lambda x: x[1], reverse=True),
                      columns=["reason", "count"])
    rc.to_csv(os.path.join(OUT_ROOT, "reason_breakdown.tsv"), sep="\t", index=False)

    print("[INFO] Complete. Main outputs:")
    print(" - blacklist.tsv")
    print(" - summary_counts.tsv")
    print(" - reason_breakdown.tsv")
    print(" - <sample>/ivar_only_final.tsv")
    print(" - <sample>/<sample>.ivar_mid.tsv")

if __name__ == "__main__":
    main()
