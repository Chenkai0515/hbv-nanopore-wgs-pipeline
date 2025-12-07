"""
HBV Pipeline - Core Library

A collection of reusable functions for HBV genome analysis:
- Coordinate transformation between rotated and original reference
- PAF file parsing and interval operations
- Variant tiering and quality assessment
"""

from .coordinates import (
    rotated_to_original,
    original_to_rotated,
    transform_positions,
    HBV_GENOME_LENGTH,
    HBV_ROTATION_OFFSET,
)

from .paf_utils import (
    parse_paf_line,
    parse_paf_file,
    merge_intervals,
    classify_read_by_alignment,
)

from .variants_utils import (
    assign_tier,
    build_blacklist,
    merge_variant_calls,
)

__version__ = "1.0.0"
__author__ = "Chenkai Jiang"

__all__ = [
    # Coordinates
    "rotated_to_original",
    "original_to_rotated",
    "transform_positions",
    "HBV_GENOME_LENGTH",
    "HBV_ROTATION_OFFSET",
    # PAF utils
    "parse_paf_line",
    "parse_paf_file",
    "merge_intervals",
    "classify_read_by_alignment",
    # Variants
    "assign_tier",
    "build_blacklist",
    "merge_variant_calls",
]
