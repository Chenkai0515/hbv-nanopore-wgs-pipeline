"""
PAF File Utilities for HBV Pipeline

This module provides functions for parsing PAF (Pairwise mApping Format) files
and classifying reads based on alignment characteristics.

PAF Format (minimap2 output):
    Col 1:  Query name
    Col 2:  Query length
    Col 3:  Query start (0-based)
    Col 4:  Query end
    Col 5:  Strand (+/-)
    Col 6:  Target name
    Col 7:  Target length
    Col 8:  Target start
    Col 9:  Target end
    Col 10: Number of matching bases
    Col 11: Alignment block length
    Col 12: Mapping quality

Common use cases:
- Classify reads as viral/host based on alignment patterns
- Extract HBV-aligned regions for masking
- Identify chimeric reads and host tails
"""

import gzip
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Iterator, TextIO


@dataclass
class PAFRecord:
    """Parsed PAF alignment record."""
    query_name: str
    query_length: int
    query_start: int  # 0-based
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    block_length: int
    mapq: int
    
    @property
    def query_aligned_length(self) -> int:
        """Length of query region in alignment."""
        return self.query_end - self.query_start
    
    @property
    def target_aligned_length(self) -> int:
        """Length of target region in alignment."""
        return self.target_end - self.target_start
    
    @property
    def identity(self) -> float:
        """Percent identity (matches / block_length)."""
        if self.block_length == 0:
            return 0.0
        return self.num_matches / self.block_length
    
    @property
    def query_coverage(self) -> float:
        """Fraction of query covered by alignment."""
        if self.query_length == 0:
            return 0.0
        return self.query_aligned_length / self.query_length


def parse_paf_line(line: str) -> Optional[PAFRecord]:
    """
    Parse a single PAF line into a PAFRecord.
    
    Args:
        line: Tab-separated PAF line
    
    Returns:
        PAFRecord if parsing successful, None otherwise
    
    Examples:
        >>> record = parse_paf_line("read1\\t1000\\t0\\t950\\t+\\tHBV\\t3221\\t100\\t1050\\t900\\t950\\t60")
        >>> record.query_name
        'read1'
        >>> record.mapq
        60
    """
    parts = line.strip().split('\t')
    if len(parts) < 12:
        return None
    
    try:
        return PAFRecord(
            query_name=parts[0],
            query_length=int(parts[1]),
            query_start=int(parts[2]),
            query_end=int(parts[3]),
            strand=parts[4],
            target_name=parts[5],
            target_length=int(parts[6]),
            target_start=int(parts[7]),
            target_end=int(parts[8]),
            num_matches=int(parts[9]),
            block_length=int(parts[10]),
            mapq=int(parts[11])
        )
    except (ValueError, IndexError):
        return None


def parse_paf_file(
    paf_path: str,
    min_mapq: int = 0,
    min_match_len: int = 0
) -> Iterator[PAFRecord]:
    """
    Parse a PAF file and yield records that pass filters.
    
    Args:
        paf_path: Path to PAF file (supports .gz)
        min_mapq: Minimum mapping quality filter
        min_match_len: Minimum alignment length filter
    
    Yields:
        PAFRecord objects passing filters
    
    Examples:
        >>> for record in parse_paf_file("alignments.paf", min_mapq=20):
        ...     print(record.query_name, record.mapq)
    """
    opener = gzip.open if paf_path.endswith('.gz') else open
    
    with opener(paf_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            record = parse_paf_line(line)
            if record is None:
                continue
            
            if record.mapq < min_mapq:
                continue
            
            if record.query_aligned_length < min_match_len:
                continue
            
            yield record


def merge_intervals(
    intervals: List[Tuple[int, int]],
    merge_distance: int = 0
) -> List[Tuple[int, int]]:
    """
    Merge overlapping or nearby intervals.
    
    Args:
        intervals: List of (start, end) tuples (0-based, half-open)
        merge_distance: Merge intervals within this distance
    
    Returns:
        List of merged (start, end) tuples
    
    Examples:
        >>> merge_intervals([(0, 100), (90, 200), (250, 300)])
        [(0, 200), (250, 300)]
        >>> merge_intervals([(0, 100), (110, 200)], merge_distance=10)
        [(0, 200)]
    """
    if not intervals:
        return []
    
    # Sort by start position
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [list(sorted_ivs[0])]
    
    for start, end in sorted_ivs[1:]:
        prev_start, prev_end = merged[-1]
        
        # Check if current interval overlaps or is within merge_distance
        if start <= prev_end + merge_distance:
            # Extend the previous interval
            merged[-1][1] = max(prev_end, end)
        else:
            # Start a new interval
            merged.append([start, end])
    
    return [(s, e) for s, e in merged]


def group_alignments_by_query(
    records: List[PAFRecord]
) -> Dict[str, List[PAFRecord]]:
    """
    Group PAF records by query name.
    
    Args:
        records: List of PAFRecord objects
    
    Returns:
        Dictionary mapping query name to list of alignments
    """
    by_query: Dict[str, List[PAFRecord]] = {}
    for record in records:
        if record.query_name not in by_query:
            by_query[record.query_name] = []
        by_query[record.query_name].append(record)
    return by_query


@dataclass
class ReadClassification:
    """Classification result for a read based on alignments."""
    read_id: str
    read_length: int
    classification: str  # 'viral', 'host', 'chimeric', 'unaligned'
    viral_intervals: List[Tuple[int, int]]
    host_intervals: List[Tuple[int, int]]
    viral_fraction: float
    host_fraction: float


def classify_read_by_alignment(
    read_id: str,
    read_length: int,
    viral_intervals: List[Tuple[int, int]],
    host_intervals: List[Tuple[int, int]],
    viral_targets: Optional[List[str]] = None,
    host_targets: Optional[List[str]] = None,
    min_viral_frac: float = 0.5,
    min_host_frac: float = 0.1,
) -> ReadClassification:
    """
    Classify a read based on viral and host alignments.
    
    Args:
        read_id: Read identifier
        read_length: Total read length
        viral_intervals: List of (start, end) aligned to viral reference
        host_intervals: List of (start, end) aligned to host reference
        viral_targets: Names of viral reference sequences (optional)
        host_targets: Names of host reference sequences (optional)
        min_viral_frac: Minimum viral coverage to classify as viral
        min_host_frac: Minimum host coverage to classify as host
    
    Returns:
        ReadClassification with classification result
    
    Classification Rules:
        - viral: >50% aligned to viral, <10% to host
        - host: >10% aligned to host, <50% to viral
        - chimeric: Significant alignments to both
        - unaligned: Neither meets threshold
    """
    # Merge overlapping intervals
    viral_merged = merge_intervals(viral_intervals)
    host_merged = merge_intervals(host_intervals)
    
    # Calculate coverage fractions
    viral_bases = sum(e - s for s, e in viral_merged)
    host_bases = sum(e - s for s, e in host_merged)
    
    viral_frac = viral_bases / read_length if read_length > 0 else 0
    host_frac = host_bases / read_length if read_length > 0 else 0
    
    # Determine classification
    if viral_frac >= min_viral_frac and host_frac < min_host_frac:
        classification = 'viral'
    elif host_frac >= min_host_frac and viral_frac < min_viral_frac:
        classification = 'host'
    elif viral_frac >= min_viral_frac and host_frac >= min_host_frac:
        classification = 'chimeric'
    else:
        classification = 'unaligned'
    
    return ReadClassification(
        read_id=read_id,
        read_length=read_length,
        classification=classification,
        viral_intervals=viral_merged,
        host_intervals=host_merged,
        viral_fraction=viral_frac,
        host_fraction=host_frac
    )


def extract_query_intervals(
    records: List[PAFRecord],
    min_mapq: int = 20,
    merge_distance: int = 10
) -> Dict[str, List[Tuple[int, int]]]:
    """
    Extract query-coordinate intervals from PAF records.
    
    Args:
        records: List of PAFRecord objects
        min_mapq: Minimum mapping quality
        merge_distance: Distance for merging nearby intervals
    
    Returns:
        Dictionary mapping read_id to list of (start, end) intervals
    """
    intervals_by_read: Dict[str, List[Tuple[int, int]]] = {}
    
    for record in records:
        if record.mapq < min_mapq:
            continue
        
        if record.query_name not in intervals_by_read:
            intervals_by_read[record.query_name] = []
        
        intervals_by_read[record.query_name].append(
            (record.query_start, record.query_end)
        )
    
    # Merge intervals for each read
    return {
        read_id: merge_intervals(ivs, merge_distance)
        for read_id, ivs in intervals_by_read.items()
    }

