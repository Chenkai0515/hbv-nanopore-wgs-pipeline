"""
Tests for PAF file parsing and interval utilities.
"""

import sys
from pathlib import Path

import pytest

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

from hbvpipe.paf_utils import (
    parse_paf_line,
    merge_intervals,
    classify_read_by_alignment,
    PAFRecord,
)


# ============================================================================
# Tests: PAF Line Parsing
# ============================================================================

class TestParsePafLine:
    """Tests for parse_paf_line function."""
    
    def test_valid_line(self):
        """Parse a valid PAF line."""
        line = "read1\t1000\t10\t950\t+\tHBV\t3221\t100\t1040\t890\t940\t60"
        record = parse_paf_line(line)
        
        assert record is not None
        assert record.query_name == "read1"
        assert record.query_length == 1000
        assert record.query_start == 10
        assert record.query_end == 950
        assert record.strand == "+"
        assert record.target_name == "HBV"
        assert record.mapq == 60
    
    def test_short_line(self):
        """Line with fewer than 12 columns should return None."""
        line = "read1\t1000\t10\t950\t+"
        record = parse_paf_line(line)
        assert record is None
    
    def test_empty_line(self):
        """Empty line should return None."""
        record = parse_paf_line("")
        assert record is None
    
    def test_invalid_numbers(self):
        """Line with non-numeric values should return None."""
        line = "read1\tXXX\t10\t950\t+\tHBV\t3221\t100\t1040\t890\t940\t60"
        record = parse_paf_line(line)
        assert record is None


class TestPAFRecordProperties:
    """Tests for PAFRecord computed properties."""
    
    @pytest.fixture
    def sample_record(self):
        return PAFRecord(
            query_name="read1",
            query_length=1000,
            query_start=100,
            query_end=900,
            strand="+",
            target_name="HBV",
            target_length=3221,
            target_start=50,
            target_end=850,
            num_matches=750,
            block_length=800,
            mapq=60
        )
    
    def test_query_aligned_length(self, sample_record):
        """Test query aligned length calculation."""
        assert sample_record.query_aligned_length == 800  # 900 - 100
    
    def test_target_aligned_length(self, sample_record):
        """Test target aligned length calculation."""
        assert sample_record.target_aligned_length == 800  # 850 - 50
    
    def test_identity(self, sample_record):
        """Test identity calculation."""
        # 750 / 800 = 0.9375
        assert abs(sample_record.identity - 0.9375) < 0.001
    
    def test_query_coverage(self, sample_record):
        """Test query coverage calculation."""
        # 800 / 1000 = 0.8
        assert abs(sample_record.query_coverage - 0.8) < 0.001


# ============================================================================
# Tests: Interval Merging
# ============================================================================

class TestMergeIntervals:
    """Tests for merge_intervals function."""
    
    def test_no_overlap(self):
        """Non-overlapping intervals stay separate."""
        intervals = [(0, 100), (200, 300), (400, 500)]
        result = merge_intervals(intervals)
        assert result == [(0, 100), (200, 300), (400, 500)]
    
    def test_simple_overlap(self):
        """Overlapping intervals are merged."""
        intervals = [(0, 100), (50, 150)]
        result = merge_intervals(intervals)
        assert result == [(0, 150)]
    
    def test_multiple_overlaps(self):
        """Multiple overlapping intervals merge into one."""
        intervals = [(0, 100), (50, 150), (140, 200)]
        result = merge_intervals(intervals)
        assert result == [(0, 200)]
    
    def test_adjacent_intervals(self):
        """Adjacent intervals (touching) are merged."""
        intervals = [(0, 100), (100, 200)]
        result = merge_intervals(intervals)
        assert result == [(0, 200)]
    
    def test_merge_distance(self):
        """Intervals within merge_distance are merged."""
        intervals = [(0, 100), (110, 200)]
        result = merge_intervals(intervals, merge_distance=10)
        assert result == [(0, 200)]
    
    def test_unsorted_input(self):
        """Unsorted input is handled correctly."""
        intervals = [(200, 300), (0, 100), (150, 250)]
        result = merge_intervals(intervals)
        assert result == [(0, 100), (150, 300)]
    
    def test_empty_input(self):
        """Empty input returns empty list."""
        result = merge_intervals([])
        assert result == []
    
    def test_single_interval(self):
        """Single interval is returned as-is."""
        intervals = [(0, 100)]
        result = merge_intervals(intervals)
        assert result == [(0, 100)]
    
    def test_contained_interval(self):
        """Interval contained within another is absorbed."""
        intervals = [(0, 300), (50, 100)]
        result = merge_intervals(intervals)
        assert result == [(0, 300)]


# ============================================================================
# Tests: Read Classification
# ============================================================================

class TestClassifyReadByAlignment:
    """Tests for classify_read_by_alignment function."""
    
    def test_viral_read(self):
        """Read with >50% viral and <10% host is classified as viral."""
        result = classify_read_by_alignment(
            read_id="read1",
            read_length=1000,
            viral_intervals=[(0, 600)],
            host_intervals=[(700, 750)],
            min_viral_frac=0.5,
            min_host_frac=0.1
        )
        assert result.classification == "viral"
        assert result.viral_fraction == 0.6
        assert result.host_fraction == 0.05
    
    def test_host_read(self):
        """Read with >10% host and <50% viral is classified as host."""
        result = classify_read_by_alignment(
            read_id="read1",
            read_length=1000,
            viral_intervals=[(0, 100)],
            host_intervals=[(200, 900)],
            min_viral_frac=0.5,
            min_host_frac=0.1
        )
        assert result.classification == "host"
    
    def test_chimeric_read(self):
        """Read with significant viral AND host is classified as chimeric."""
        result = classify_read_by_alignment(
            read_id="read1",
            read_length=1000,
            viral_intervals=[(0, 600)],
            host_intervals=[(500, 800)],  # Overlaps but both significant
            min_viral_frac=0.5,
            min_host_frac=0.1
        )
        assert result.classification == "chimeric"
    
    def test_unaligned_read(self):
        """Read with minimal alignments is classified as unaligned."""
        result = classify_read_by_alignment(
            read_id="read1",
            read_length=1000,
            viral_intervals=[(0, 50)],
            host_intervals=[(100, 150)],
            min_viral_frac=0.5,
            min_host_frac=0.1
        )
        assert result.classification == "unaligned"
    
    def test_overlapping_intervals_merged(self):
        """Overlapping intervals are merged before calculation."""
        result = classify_read_by_alignment(
            read_id="read1",
            read_length=1000,
            viral_intervals=[(0, 300), (200, 500), (400, 600)],
            host_intervals=[],
            min_viral_frac=0.5,
            min_host_frac=0.1
        )
        # Merged: (0, 600) = 600 bases = 60%
        assert result.classification == "viral"
        assert result.viral_fraction == 0.6

