"""
Tests for coordinate transformation functions.

The HBV pipeline uses a rotated reference (starting at DR1, offset 1823).
These tests verify the coordinate transformation between:
- Rotated reference coordinates (used in pipeline)
- Original/standard coordinates (used in databases)
"""

import sys
from pathlib import Path

import pytest

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

from hbvpipe.coordinates import rotated_to_original, original_to_rotated


# ============================================================================
# Tests: Rotated to Original Conversion
# ============================================================================

class TestRotatedToOriginal:
    """Tests for rotated_to_original conversion."""
    
    def test_position_one(self):
        """Position 1 in rotated = position 1824 in original."""
        assert rotated_to_original(1) == 1824
    
    def test_position_100(self):
        """Position 100 in rotated = position 1923 in original."""
        assert rotated_to_original(100) == 1923
    
    def test_wraparound(self):
        """Test position that wraps around the genome."""
        # Position 1399 in rotated = (1398 + 1823) % 3221 + 1 = 1
        assert rotated_to_original(1399) == 1
    
    def test_last_position(self):
        """Last position in rotated = position 1823 in original."""
        assert rotated_to_original(3221) == 1823
    
    def test_bcp_position_1762(self):
        """Test BCP mutation position 1762 (standard)."""
        # Original 1762 should map from rotated 3160
        # Verify: (3160 - 1 + 1823) % 3221 + 1 = 1762
        assert rotated_to_original(3160) == 1762
    
    def test_bcp_position_1764(self):
        """Test BCP mutation position 1764 (standard)."""
        assert rotated_to_original(3162) == 1764
    
    def test_precore_position_1896(self):
        """Test PreCore mutation position 1896 (standard)."""
        # Original 1896 should map from rotated 73
        # Verify: (73 - 1 + 1823) % 3221 + 1 = 1896
        assert rotated_to_original(73) == 1896
    
    def test_invalid_position_zero(self):
        """Position 0 should raise ValueError."""
        with pytest.raises(ValueError):
            rotated_to_original(0)
    
    def test_invalid_position_negative(self):
        """Negative position should raise ValueError."""
        with pytest.raises(ValueError):
            rotated_to_original(-1)
    
    def test_invalid_position_too_large(self):
        """Position > genome_length should raise ValueError."""
        with pytest.raises(ValueError):
            rotated_to_original(3222)


# ============================================================================
# Tests: Original to Rotated Conversion
# ============================================================================

class TestOriginalToRotated:
    """Tests for original_to_rotated conversion."""
    
    def test_position_1824(self):
        """Original position 1824 = rotated position 1."""
        assert original_to_rotated(1824) == 1
    
    def test_position_1(self):
        """Original position 1 = rotated position 1399."""
        assert original_to_rotated(1) == 1399
    
    def test_position_1823(self):
        """Original position 1823 = rotated position 3221."""
        assert original_to_rotated(1823) == 3221
    
    def test_bcp_position_1762(self):
        """BCP position 1762 should map to rotated 3160."""
        assert original_to_rotated(1762) == 3160
    
    def test_precore_position_1896(self):
        """PreCore position 1896 should map to rotated 73."""
        assert original_to_rotated(1896) == 73


# ============================================================================
# Tests: Round-Trip Conversion
# ============================================================================

class TestRoundTrip:
    """Test that conversions are reversible."""
    
    def test_roundtrip_all_positions(self):
        """Every position should survive a round-trip conversion."""
        genome_length = 3221
        
        for pos in range(1, genome_length + 1):
            # Rotated -> Original -> Rotated
            original = rotated_to_original(pos)
            back_to_rotated = original_to_rotated(original)
            assert back_to_rotated == pos, f"Failed at rotated position {pos}"
    
    def test_roundtrip_original_all_positions(self):
        """Every original position should survive a round-trip."""
        genome_length = 3221
        
        for pos in range(1, genome_length + 1):
            # Original -> Rotated -> Original
            rotated = original_to_rotated(pos)
            back_to_original = rotated_to_original(rotated)
            assert back_to_original == pos, f"Failed at original position {pos}"


# ============================================================================
# Tests: Edge Cases
# ============================================================================

class TestEdgeCases:
    """Test edge cases and boundary conditions."""
    
    def test_custom_offset(self):
        """Test with a custom offset value."""
        # With offset 0, positions should be unchanged
        assert rotated_to_original(100, offset=0) == 100
    
    def test_custom_genome_length(self):
        """Test with a custom genome length."""
        # Smaller genome
        assert rotated_to_original(50, offset=10, genome_length=100) == 60
    
    def test_offset_equals_genome_length(self):
        """When offset equals genome length, positions unchanged."""
        genome_length = 100
        for pos in range(1, genome_length + 1):
            result = rotated_to_original(pos, offset=genome_length, 
                                         genome_length=genome_length)
            assert result == pos


# ============================================================================
# Tests: Known Mutation Sites
# ============================================================================

class TestKnownMutationSites:
    """Test conversion of clinically relevant mutation sites."""
    
    @pytest.mark.parametrize("original,rotated", [
        (1762, 3160),  # BCP A1762T
        (1764, 3162),  # BCP G1764A
        (1896, 73),    # PreCore G1896A
        (1899, 76),    # PreCore G1899A
        (1824, 1),     # DR1 (rotation start)
    ])
    def test_clinical_positions(self, original, rotated):
        """Test clinically relevant position conversions."""
        assert original_to_rotated(original) == rotated
        assert rotated_to_original(rotated) == original

