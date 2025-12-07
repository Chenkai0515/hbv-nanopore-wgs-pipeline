"""
Coordinate Transformation for HBV Genome Analysis

The HBV pipeline uses a rotated reference sequence that starts at DR1 (position 1824
in the standard X02763 reference). This module provides functions to convert between:

- Rotated coordinates: Used internally by the pipeline (position 1 = DR1)
- Original coordinates: Standard database coordinates (X02763, GenBank)

Mathematical Background:
------------------------
The HBV genome is circular (~3.2kb). For linear analysis, we "cut" the circle
at DR1 (position 1824). This creates a rotated reference where:

    rotated_pos = 1  corresponds to  original_pos = 1824
    rotated_pos = 1399  corresponds to  original_pos = 1 (wraparound)
    rotated_pos = 3221  corresponds to  original_pos = 1823

Conversion Formulas:
    rotated -> original: ((rotated_pos - 1) + offset) % genome_length + 1
    original -> rotated: ((original_pos - 1) - offset) % genome_length + 1

where offset = 1823 (the position before DR1 in original coordinates)
"""

from typing import List, Tuple, Union

# HBV genome constants
HBV_GENOME_LENGTH: int = 3221  # Standard HBV genome length
HBV_ROTATION_OFFSET: int = 1823  # DR1 position - 1 in original coordinates


def rotated_to_original(
    rotated_pos: int,
    offset: int = HBV_ROTATION_OFFSET,
    genome_length: int = HBV_GENOME_LENGTH
) -> int:
    """
    Convert rotated reference position to original (standard) position.
    
    This transforms pipeline-internal coordinates back to standard database
    coordinates (e.g., X02763, GenBank).
    
    Args:
        rotated_pos: Position in rotated reference (1-based)
        offset: Rotation offset, default 1823 for HBV
        genome_length: Total genome length, default 3221
    
    Returns:
        Position in original reference (1-based)
    
    Raises:
        ValueError: If position is out of valid range [1, genome_length]
    
    Examples:
        >>> rotated_to_original(1)  # DR1 start
        1824
        >>> rotated_to_original(73)  # PreCore G1896A site
        1896
        >>> rotated_to_original(3160)  # BCP A1762T site
        1762
    """
    if rotated_pos < 1 or rotated_pos > genome_length:
        raise ValueError(
            f"Position {rotated_pos} out of valid range [1, {genome_length}]"
        )
    
    original = ((rotated_pos - 1) + offset) % genome_length + 1
    return original


def original_to_rotated(
    original_pos: int,
    offset: int = HBV_ROTATION_OFFSET,
    genome_length: int = HBV_GENOME_LENGTH
) -> int:
    """
    Convert original (standard) position to rotated reference position.
    
    This transforms standard database coordinates to pipeline-internal
    coordinates.
    
    Args:
        original_pos: Position in original reference (1-based)
        offset: Rotation offset, default 1823 for HBV
        genome_length: Total genome length, default 3221
    
    Returns:
        Position in rotated reference (1-based)
    
    Raises:
        ValueError: If position is out of valid range [1, genome_length]
    
    Examples:
        >>> original_to_rotated(1824)  # DR1 start
        1
        >>> original_to_rotated(1896)  # PreCore G1896A site
        73
        >>> original_to_rotated(1762)  # BCP A1762T site
        3160
    """
    if original_pos < 1 or original_pos > genome_length:
        raise ValueError(
            f"Position {original_pos} out of valid range [1, {genome_length}]"
        )
    
    rotated = ((original_pos - 1) - offset) % genome_length + 1
    return rotated


def transform_positions(
    positions: List[int],
    to_original: bool = True,
    offset: int = HBV_ROTATION_OFFSET,
    genome_length: int = HBV_GENOME_LENGTH
) -> List[int]:
    """
    Transform a list of positions between coordinate systems.
    
    Args:
        positions: List of positions to transform
        to_original: If True, convert rotated->original; else original->rotated
        offset: Rotation offset
        genome_length: Genome length
    
    Returns:
        List of transformed positions
    
    Examples:
        >>> transform_positions([1, 73, 3160], to_original=True)
        [1824, 1896, 1762]
    """
    if to_original:
        return [rotated_to_original(p, offset, genome_length) for p in positions]
    else:
        return [original_to_rotated(p, offset, genome_length) for p in positions]


def transform_interval(
    start: int,
    end: int,
    to_original: bool = True,
    offset: int = HBV_ROTATION_OFFSET,
    genome_length: int = HBV_GENOME_LENGTH
) -> Tuple[int, int]:
    """
    Transform a coordinate interval between coordinate systems.
    
    Note: For intervals that span the wraparound point, the result may need
    special handling as the transformed start could be > end.
    
    Args:
        start: Interval start position (1-based, inclusive)
        end: Interval end position (1-based, inclusive)
        to_original: Direction of transformation
        offset: Rotation offset
        genome_length: Genome length
    
    Returns:
        Tuple of (transformed_start, transformed_end)
    """
    if to_original:
        return (
            rotated_to_original(start, offset, genome_length),
            rotated_to_original(end, offset, genome_length)
        )
    else:
        return (
            original_to_rotated(start, offset, genome_length),
            original_to_rotated(end, offset, genome_length)
        )


# Clinical mutation sites mapping (original -> rotated)
CLINICAL_SITES = {
    # BCP mutations
    1762: {"rotated": 3160, "name": "BCP_A1762T"},
    1764: {"rotated": 3162, "name": "BCP_G1764A"},
    # PreCore mutations
    1896: {"rotated": 73, "name": "PreC_G1896A"},
    1899: {"rotated": 76, "name": "PreC_G1899A"},
    # DR1 (rotation point)
    1824: {"rotated": 1, "name": "DR1_start"},
}

