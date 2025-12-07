"""
Tests for variant quality assessment and tiering utilities.
"""

import sys
from pathlib import Path

import pytest

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

from hbvpipe.variants_utils import (
    Variant,
    TierResult,
    assign_tier,
    build_blacklist,
    merge_variant_calls,
)


# ============================================================================
# Tests: Variant Class
# ============================================================================


class TestVariant:
    """Tests for Variant dataclass."""

    def test_snp_detection(self):
        """SNP is correctly identified."""
        v = Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")
        assert v.is_snp is True
        assert v.is_indel is False
        assert v.var_type == "SNP"

    def test_deletion_detection(self):
        """Deletion is correctly identified."""
        v = Variant("HBV", 100, "AT", "A", 0.15, 400, "clair3")
        assert v.is_snp is False
        assert v.is_indel is True
        assert v.var_type == "DEL"

    def test_insertion_detection(self):
        """Insertion is correctly identified."""
        v = Variant("HBV", 100, "A", "ATG", 0.10, 300, "both")
        assert v.is_snp is False
        assert v.is_indel is True
        assert v.var_type == "INS"


# ============================================================================
# Tests: Tier Assignment
# ============================================================================


class TestAssignTier:
    """Tests for assign_tier function."""

    def test_high_confidence_snp(self):
        """AF >= 20% and good depth is HC."""
        v = Variant("HBV", 100, "A", "G", 0.25, 500, "both")
        result = assign_tier(v)
        assert result.tier == "SNP-HC"
        assert result.decision == "KEEP"

    def test_medium_confidence_snp(self):
        """AF 5-20% is MC."""
        v = Variant("HBV", 100, "A", "G", 0.10, 500, "ivar")
        result = assign_tier(v)
        assert result.tier == "SNP-MC"
        assert result.decision == "KEEP"

    def test_low_frequency_snp(self):
        """AF 1-5% is LF (CANDIDATE)."""
        v = Variant("HBV", 100, "A", "G", 0.02, 500, "ivar")
        result = assign_tier(v)
        assert result.tier == "SNP-LF"
        assert result.decision == "CANDIDATE"

    def test_below_minimum_af(self):
        """AF < 1% is rejected."""
        v = Variant("HBV", 100, "A", "G", 0.005, 500, "ivar")
        result = assign_tier(v)
        assert result.tier == "SNP-LF"
        assert result.decision == "REJECT"

    def test_consensus_flipped(self):
        """Consensus-level variant (AF ~1.0) is FLIPPED."""
        v = Variant("HBV", 100, "A", "G", 0.98, 500, "consensus")
        result = assign_tier(v)
        assert "FLIPPED" in result.tier
        assert result.decision == "KEEP"

    def test_low_depth(self):
        """Low depth results in CANDIDATE."""
        v = Variant("HBV", 100, "A", "G", 0.25, 5, "both")
        result = assign_tier(v, min_depth=10)
        assert result.decision == "CANDIDATE"
        assert "depth" in result.reason.lower()

    def test_indel_tiering(self):
        """Indels get correct tier prefix."""
        v = Variant("HBV", 100, "AT", "A", 0.15, 500, "ivar")
        result = assign_tier(v)
        assert result.tier == "DEL-MC"

    def test_custom_thresholds(self):
        """Custom AF thresholds work correctly."""
        v = Variant("HBV", 100, "A", "G", 0.15, 500, "both")

        # With default thresholds, 15% is MC
        result_default = assign_tier(v)
        assert result_default.tier == "SNP-MC"

        # With lower HC threshold, 15% becomes HC
        result_custom = assign_tier(v, hc_min_af=0.10)
        assert result_custom.tier == "SNP-HC"

    def test_require_both_callers(self):
        """When requiring both callers, single caller is CANDIDATE."""
        v = Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")
        result = assign_tier(v, require_both_callers=True)
        assert result.tier == "SNP-HC"
        assert result.decision == "CANDIDATE"

        v_both = Variant("HBV", 100, "A", "G", 0.25, 500, "both")
        result_both = assign_tier(v_both, require_both_callers=True)
        assert result_both.decision == "KEEP"


# ============================================================================
# Tests: Blacklist Building
# ============================================================================


class TestBuildBlacklist:
    """Tests for build_blacklist function (Panel of Normals)."""

    def test_frequent_low_af_variant(self):
        """Variant in >30% of samples at low AF is blacklisted."""
        samples = {
            "S1": [Variant("HBV", 100, "A", "G", 0.02, 500, "ivar")],
            "S2": [Variant("HBV", 100, "A", "G", 0.03, 600, "ivar")],
            "S3": [Variant("HBV", 100, "A", "G", 0.01, 400, "ivar")],
            "S4": [Variant("HBV", 100, "A", "G", 0.02, 500, "ivar")],
        }
        blacklist = build_blacklist(samples, min_cohort_freq=0.5)
        assert ("HBV", 100, "A", "G") in blacklist

    def test_infrequent_variant_not_blacklisted(self):
        """Variant in <30% of samples is not blacklisted."""
        samples = {
            "S1": [Variant("HBV", 100, "A", "G", 0.02, 500, "ivar")],
            "S2": [],
            "S3": [],
            "S4": [],
        }
        blacklist = build_blacklist(samples, min_cohort_freq=0.3)
        assert ("HBV", 100, "A", "G") not in blacklist

    def test_high_af_variant_not_blacklisted(self):
        """Variant at high AF is not blacklisted (likely real)."""
        samples = {
            "S1": [Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")],
            "S2": [Variant("HBV", 100, "A", "G", 0.30, 600, "ivar")],
            "S3": [Variant("HBV", 100, "A", "G", 0.20, 400, "ivar")],
        }
        blacklist = build_blacklist(samples, min_cohort_freq=0.5, max_af_threshold=0.05)
        assert ("HBV", 100, "A", "G") not in blacklist

    def test_empty_samples(self):
        """Empty samples returns empty blacklist."""
        blacklist = build_blacklist({})
        assert blacklist == set()


# ============================================================================
# Tests: Variant Merging
# ============================================================================


class TestMergeVariantCalls:
    """Tests for merge_variant_calls function."""

    def test_same_variant_both_callers(self):
        """Same variant from both callers is merged with caller='both'."""
        ivar = [Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")]
        clair3 = [Variant("HBV", 100, "A", "G", 0.23, 480, "clair3")]

        merged = merge_variant_calls(ivar, clair3)

        assert len(merged) == 1
        assert merged[0].caller == "both"
        assert merged[0].af == 0.25  # Takes higher AF

    def test_different_variants(self):
        """Different variants remain separate."""
        ivar = [Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")]
        clair3 = [Variant("HBV", 200, "T", "C", 0.30, 600, "clair3")]

        merged = merge_variant_calls(ivar, clair3)

        assert len(merged) == 2
        assert merged[0].pos == 100
        assert merged[0].caller == "ivar"
        assert merged[1].pos == 200
        assert merged[1].caller == "clair3"

    def test_ivar_only_variant(self):
        """Variant from iVar only has caller='ivar'."""
        ivar = [Variant("HBV", 100, "A", "G", 0.25, 500, "ivar")]
        clair3 = []

        merged = merge_variant_calls(ivar, clair3)

        assert len(merged) == 1
        assert merged[0].caller == "ivar"

    def test_clair3_only_variant(self):
        """Variant from Clair3 only has caller='clair3'."""
        ivar = []
        clair3 = [Variant("HBV", 100, "A", "G", 0.25, 500, "clair3")]

        merged = merge_variant_calls(ivar, clair3)

        assert len(merged) == 1
        assert merged[0].caller == "clair3"

    def test_sorted_output(self):
        """Output is sorted by position."""
        ivar = [Variant("HBV", 300, "A", "G", 0.25, 500, "ivar")]
        clair3 = [
            Variant("HBV", 100, "T", "C", 0.30, 600, "clair3"),
            Variant("HBV", 200, "G", "A", 0.20, 400, "clair3"),
        ]

        merged = merge_variant_calls(ivar, clair3)

        positions = [v.pos for v in merged]
        assert positions == [100, 200, 300]
