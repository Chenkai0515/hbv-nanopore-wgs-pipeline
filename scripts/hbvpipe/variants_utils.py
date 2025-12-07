"""
Variant Quality Assessment and Tiering Utilities

This module provides functions for assessing variant quality and assigning
confidence tiers based on the dual-caller approach (iVar + Clair3).

Tiering System:
--------------
Variants are classified into confidence tiers based on:
- Allele frequency (AF)
- Caller agreement
- Coverage depth
- Strand bias

Tiers:
- HC (High Confidence): AF >= 20%, good evidence
- MC (Medium Confidence): AF >= 5%, moderate evidence  
- LF (Low Frequency): AF >= 1%, weak evidence
- FLIPPED: Consensus differs from reference (AF effectively 1.0)

Blacklist (Panel of Normals):
----------------------------
Systematic artifacts are identified by:
- High frequency across samples (>30% of cohort)
- Low AF in individual samples (<5%)
These positions are flagged as potential artifacts.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class Variant:
    """Represents a single variant call."""
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float
    depth: int
    caller: str  # 'ivar', 'clair3', 'both', 'consensus'
    
    @property
    def is_snp(self) -> bool:
        """Check if variant is a SNP."""
        return len(self.ref) == 1 and len(self.alt) == 1
    
    @property
    def is_indel(self) -> bool:
        """Check if variant is an indel."""
        return len(self.ref) != len(self.alt)
    
    @property
    def var_type(self) -> str:
        """Get variant type string."""
        if self.is_snp:
            return "SNP"
        elif len(self.ref) > len(self.alt):
            return "DEL"
        else:
            return "INS"


@dataclass
class TierResult:
    """Result of tier assignment."""
    tier: str  # 'HC', 'MC', 'LF', 'FLIPPED'
    decision: str  # 'KEEP', 'CANDIDATE', 'REJECT'
    reason: str


def assign_tier(
    variant: Variant,
    hc_min_af: float = 0.20,
    mc_min_af: float = 0.05,
    lf_min_af: float = 0.01,
    min_depth: int = 10,
    require_both_callers: bool = False
) -> TierResult:
    """
    Assign confidence tier to a variant.
    
    Args:
        variant: Variant object to tier
        hc_min_af: Minimum AF for High Confidence (default 20%)
        mc_min_af: Minimum AF for Medium Confidence (default 5%)
        lf_min_af: Minimum AF for Low Frequency (default 1%)
        min_depth: Minimum read depth required
        require_both_callers: Require agreement from both callers for HC
    
    Returns:
        TierResult with tier, decision, and reason
    
    Examples:
        >>> v = Variant("HBV", 1762, "A", "T", 0.25, 500, "both")
        >>> result = assign_tier(v)
        >>> result.tier
        'HC'
        >>> result.decision
        'KEEP'
    """
    af = variant.af
    depth = variant.depth
    caller = variant.caller
    var_type = variant.var_type
    
    # Check if this is a consensus-level difference (flipped)
    if caller == 'consensus' or af >= 0.95:
        return TierResult(
            tier=f"{var_type}-FLIPPED",
            decision="KEEP",
            reason="Consensus differs from reference"
        )
    
    # Check minimum depth
    if depth < min_depth:
        return TierResult(
            tier=f"{var_type}-LF",
            decision="CANDIDATE",
            reason=f"Low depth ({depth} < {min_depth})"
        )
    
    # Assign tier based on AF
    if af >= hc_min_af:
        tier = f"{var_type}-HC"
        
        # For HC, optionally require both callers
        if require_both_callers and caller != 'both':
            return TierResult(
                tier=tier,
                decision="CANDIDATE",
                reason=f"Single caller only ({caller})"
            )
        
        return TierResult(
            tier=tier,
            decision="KEEP",
            reason=f"High confidence (AF={af:.2%})"
        )
    
    elif af >= mc_min_af:
        return TierResult(
            tier=f"{var_type}-MC",
            decision="KEEP",
            reason=f"Medium confidence (AF={af:.2%})"
        )
    
    elif af >= lf_min_af:
        return TierResult(
            tier=f"{var_type}-LF",
            decision="CANDIDATE",
            reason=f"Low frequency (AF={af:.2%})"
        )
    
    else:
        return TierResult(
            tier=f"{var_type}-LF",
            decision="REJECT",
            reason=f"Below minimum AF ({af:.2%} < {lf_min_af:.0%})"
        )


def build_blacklist(
    variants_by_sample: Dict[str, List[Variant]],
    min_cohort_freq: float = 0.3,
    max_af_threshold: float = 0.05
) -> Set[Tuple[str, int, str, str]]:
    """
    Build a blacklist of likely artifact positions (Panel of Normals approach).
    
    Identifies positions that appear frequently across samples at low AF,
    suggesting systematic artifacts rather than true variants.
    
    Args:
        variants_by_sample: Dict mapping sample_id to list of variants
        min_cohort_freq: Minimum fraction of samples with variant (default 30%)
        max_af_threshold: Maximum AF to consider as artifact (default 5%)
    
    Returns:
        Set of (chrom, pos, ref, alt) tuples to blacklist
    
    Examples:
        >>> samples = {
        ...     "S1": [Variant("HBV", 100, "A", "G", 0.02, 500, "ivar")],
        ...     "S2": [Variant("HBV", 100, "A", "G", 0.03, 600, "ivar")],
        ...     "S3": [Variant("HBV", 100, "A", "G", 0.01, 400, "ivar")],
        ... }
        >>> blacklist = build_blacklist(samples, min_cohort_freq=0.5)
        >>> ("HBV", 100, "A", "G") in blacklist
        True
    """
    n_samples = len(variants_by_sample)
    if n_samples == 0:
        return set()
    
    # Count occurrences of each variant position
    variant_counts: Dict[Tuple[str, int, str, str], List[float]] = {}
    
    for sample_id, variants in variants_by_sample.items():
        for v in variants:
            key = (v.chrom, v.pos, v.ref, v.alt)
            if key not in variant_counts:
                variant_counts[key] = []
            variant_counts[key].append(v.af)
    
    # Identify blacklist positions
    blacklist: Set[Tuple[str, int, str, str]] = set()
    
    for key, af_values in variant_counts.items():
        sample_freq = len(af_values) / n_samples
        
        # Check if appears in enough samples at low AF
        if sample_freq >= min_cohort_freq:
            # Check if most occurrences are at low AF
            low_af_count = sum(1 for af in af_values if af <= max_af_threshold)
            if low_af_count / len(af_values) >= 0.5:
                blacklist.add(key)
    
    return blacklist


def merge_variant_calls(
    ivar_variants: List[Variant],
    clair3_variants: List[Variant],
    position_tolerance: int = 0
) -> List[Variant]:
    """
    Merge variant calls from iVar and Clair3.
    
    Variants at the same position are merged, with preference for:
    1. Both callers agree (highest confidence)
    2. Single caller with higher AF
    
    Args:
        ivar_variants: Variants from iVar
        clair3_variants: Variants from Clair3
        position_tolerance: Allow this many bp difference for matching
    
    Returns:
        List of merged Variant objects with updated caller field
    
    Examples:
        >>> ivar = [Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")]
        >>> clair3 = [Variant("HBV", 100, "A", "G", 0.23, 480, "clair3")]
        >>> merged = merge_variant_calls(ivar, clair3)
        >>> merged[0].caller
        'both'
    """
    merged: Dict[Tuple[str, int, str, str], Variant] = {}
    
    # Add iVar variants
    for v in ivar_variants:
        key = (v.chrom, v.pos, v.ref, v.alt)
        merged[key] = Variant(
            chrom=v.chrom,
            pos=v.pos,
            ref=v.ref,
            alt=v.alt,
            af=v.af,
            depth=v.depth,
            caller='ivar'
        )
    
    # Merge with Clair3 variants
    for v in clair3_variants:
        key = (v.chrom, v.pos, v.ref, v.alt)
        
        if key in merged:
            # Both callers found this variant
            existing = merged[key]
            merged[key] = Variant(
                chrom=v.chrom,
                pos=v.pos,
                ref=v.ref,
                alt=v.alt,
                af=max(existing.af, v.af),  # Take higher AF
                depth=max(existing.depth, v.depth),
                caller='both'
            )
        else:
            # Clair3 only
            merged[key] = Variant(
                chrom=v.chrom,
                pos=v.pos,
                ref=v.ref,
                alt=v.alt,
                af=v.af,
                depth=v.depth,
                caller='clair3'
            )
    
    # Sort by position
    return sorted(merged.values(), key=lambda v: (v.chrom, v.pos))


def filter_by_blacklist(
    variants: List[Variant],
    blacklist: Set[Tuple[str, int, str, str]],
    mark_only: bool = True
) -> List[Tuple[Variant, bool]]:
    """
    Filter or mark variants that appear in blacklist.
    
    Args:
        variants: List of variants to filter
        blacklist: Set of (chrom, pos, ref, alt) to filter
        mark_only: If True, mark but don't remove; if False, remove
    
    Returns:
        List of (variant, is_blacklisted) tuples
    """
    result = []
    for v in variants:
        key = (v.chrom, v.pos, v.ref, v.alt)
        is_blacklisted = key in blacklist
        
        if mark_only or not is_blacklisted:
            result.append((v, is_blacklisted))
    
    return result

