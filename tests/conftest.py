"""
Pytest configuration and fixtures for HBV Pipeline tests.
"""

import pytest
import tempfile
import os
from pathlib import Path


# ============================================================================
# Path Fixtures
# ============================================================================

@pytest.fixture
def project_root():
    """Return the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def scripts_dir(project_root):
    """Return the scripts directory."""
    return project_root / "scripts"


@pytest.fixture
def temp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


# ============================================================================
# Reference Data Fixtures
# ============================================================================

@pytest.fixture
def hbv_genome_length():
    """HBV genome length."""
    return 3221


@pytest.fixture
def rotation_offset():
    """Rotation offset for coordinate transformation."""
    return 1823


@pytest.fixture
def sample_reference_sequence():
    """Provide a short reference sequence for testing."""
    # First 100 bases of a typical HBV sequence
    return (
        "TTCCACTGCCTTCCACCAAACTCTGCAAGATCCCAGAGTGAGAGGCCTGTATTTCCCTGCTGGTGGCTCC"
        "AGTTCAGGAACAGTAAACCCTGTTCTGACTAC"
    )


# ============================================================================
# FASTQ Fixtures
# ============================================================================

@pytest.fixture
def sample_fastq_content():
    """Provide sample FASTQ content for testing."""
    return """\
@read_001 runid=test
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_002 runid=test
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""


@pytest.fixture
def sample_fastq_file(temp_dir, sample_fastq_content):
    """Create a temporary FASTQ file."""
    fastq_path = temp_dir / "test_sample.fastq"
    fastq_path.write_text(sample_fastq_content)
    return fastq_path


# ============================================================================
# Variant Fixtures
# ============================================================================

@pytest.fixture
def sample_variant_dict():
    """Provide a sample variant dictionary."""
    return {
        "position": 1762,
        "ref": "A",
        "alt": "T",
        "af": 0.25,
        "depth": 500,
        "source": "DUAL",
        "tier": "HC"
    }


@pytest.fixture
def sample_ivar_output():
    """Provide sample iVar output content."""
    return """\
REGION	POS	REF	ALT	REF_DP	REF_RV	REF_QUAL	ALT_DP	ALT_RV	ALT_QUAL	ALT_FREQ	TOTAL_DP	PVAL	PASS
consensus	100	A	T	400	200	35	100	50	33	0.2	500	0.001	TRUE
consensus	200	G	C	450	225	36	50	25	32	0.1	500	0.01	TRUE
"""


# ============================================================================
# Coordinate Transformation Fixtures
# ============================================================================

@pytest.fixture
def coordinate_test_cases(hbv_genome_length, rotation_offset):
    """Provide test cases for coordinate transformation."""
    return [
        # (rotated_pos, expected_original_pos)
        (1, 1824),
        (100, 1923),
        (1398, 0),  # Wraps around (should be 3221 or handled specially)
        (1399, 1),
        (3221, 1823),
    ]


# ============================================================================
# Helper Functions
# ============================================================================

@pytest.fixture
def create_test_bam(temp_dir):
    """Factory fixture to create test BAM files."""
    def _create_bam(name="test.bam"):
        # This would create a minimal BAM file for testing
        # For now, just create an empty file
        bam_path = temp_dir / name
        bam_path.touch()
        return bam_path
    return _create_bam


# ============================================================================
# Configuration Fixtures
# ============================================================================

@pytest.fixture
def sample_config():
    """Provide a sample configuration dictionary."""
    return {
        "project": {
            "work_dir": "/tmp/test_project",
            "raw_fastq_dir": "/tmp/test_fastq",
            "results_dir": "/tmp/test_results"
        },
        "reference": {
            "hbv_unified": "/path/to/reference.fasta",
            "human_ref": "/path/to/human.fa.gz",
            "hbv_genome_length": 3221,
            "hbv_rotation_offset": 1823
        },
        "resources": {
            "threads": 4,
            "java_memory_gb": 4
        },
        "variants": {
            "filter_profile": "balanced",
            "tiers": {
                "hc_min_af": 0.20,
                "mc_min_af": 0.05,
                "lf_min_af": 0.01
            }
        }
    }

