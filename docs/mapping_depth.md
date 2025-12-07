# Mapping and Depth Analysis

## Overview

This module handles alignment of viral-enriched reads to the unified HBV reference and comprehensive coverage analysis.

## Reference Alignment

### Minimap2 Settings

```bash
minimap2 -ax map-ont \     # ONT preset
         -t 16 \           # Threads
         --secondary=yes \ # Include secondary alignments
         reference.fasta \
         reads.fastq.gz
```

**Preset Selection:**
- `map-ont`: Standard ONT reads (default)
- `lr:hq`: For high-quality (Q20+) reads

### Alignment Output

Two BAM files per sample:
1. `*.primary.sorted.bam` - Primary alignments only
2. `*.all.sorted.bam` - All alignments including secondary

**Why keep secondary alignments?**
- HBV has repetitive regions
- Multi-mapping helps identify problematic areas
- Useful for copy number analysis

## Coverage Analysis

### Per-Base Depth

Using mosdepth for efficient depth calculation:

```bash
mosdepth \
    --by 1 \              # Per-base resolution
    --threads 4 \
    --flag 2308 \         # Exclude unmapped/secondary/supplementary
    output_prefix \
    input.bam
```

### Sliding Window Analysis

Coverage is analyzed in sliding windows to identify:
- Low coverage regions (gaps, deletions)
- High variance regions (amplification bias)

**Default Parameters:**
- Window size: 200 bp
- Step size: 50 bp

### Coverage Thresholds

| Category | Formula | Default |
|----------|---------|---------|
| Low coverage | < max(LOW_ABS, median × LOW_FRAC) | <100× or <20% of median |
| High variance | \|cov - median\| > MAD × MAD_MULT | >4 MAD from median |

## Homopolymer Masking

ONT sequencing has systematic errors in homopolymer regions.

### Mask Levels

| Level | Run Length | Purpose |
|-------|------------|---------|
| `ge4` | ≥4 bp | Flag for review |
| `ge6` | ≥6 bp | High error rate |
| `ge6_pad5` | ≥6 bp ± 5bp | Extended buffer |

### HBV Homopolymer Regions

Notable homopolymer regions in HBV:

| Position (rotated) | Sequence | Gene |
|--------------------|----------|------|
| ~500-510 | AAAAA | PreS2 |
| ~1500-1510 | TTTTT | Polymerase |
| ~2800-2810 | GGGG | X gene |

## Scripts

### Step 7: Reference Mapping (`7_TA1_map_V2_unmaskreads.sh`)

```bash
export BAMS_DIR="/path/to/output"
export REF="/path/to/reference.fasta"
export FASTQ_DIR="/path/to/viral_reads"

bash 7_TA1_map_V2_unmaskreads.sh
```

### Step 8: Depth Analysis (`8_TA2_depth_v2.sh`)

```bash
export BAMS_DIR="/path/to/bams"
export REF="/path/to/reference.fasta"
export OUT_DIR="/path/to/depth_output"

bash 8_TA2_depth_v2.sh
```

## Output Structure

```
TA1_map_6_V2/
├── bam/
│   ├── {sample}.primary.sorted.bam
│   ├── {sample}.primary.sorted.bam.bai
│   ├── {sample}.all.sorted.bam
│   └── {sample}.all.sorted.bam.bai
├── stats/
│   ├── {sample}.flagstat.txt
│   └── {sample}.idxstats.txt
├── hotspots/
│   └── {sample}.multimapped.tsv
└── logs/

TA2_depth_7/
├── mosdepth/
│   ├── {sample}.mosdepth.summary.txt
│   └── {sample}.per-base.bed.gz
├── windows/
│   ├── {sample}.windows.bed.gz
│   └── {sample}.windows.stats.tsv
├── flags/
│   ├── {sample}.lowcov.bed
│   └── {sample}.highvar.bed
├── homopolymers/
│   ├── ref.hpoly_ge4.bed
│   ├── ref.hpoly_ge6.bed
│   └── ref.hpoly_ge6.pad5.bed
├── summary/
│   ├── per_sample_summary.tsv
│   └── all_samples_combined.tsv
└── logs/
```

## Coverage Interpretation

### Expected Coverage Profile

For HBV amplicon sequencing:
- Overall: 500-5000× typical
- Should be relatively uniform
- Slight dips at amplicon boundaries normal

### Low Coverage Causes

| Pattern | Likely Cause | Action |
|---------|--------------|--------|
| Focal gap | Deletion | Verify with reads |
| Amplicon boundary | PCR dropout | Check primer design |
| GC-rich region | Sequencing bias | Expected |
| Homopolymer | ONT limitation | Use caution |

### High Variance Causes

| Pattern | Likely Cause | Action |
|---------|--------------|--------|
| Amplification spike | PCR bias | Normalize |
| Multi-mapped region | Repeat | Check hotspots |
| Copy number change | Biological | Investigate |

## Summary Statistics

### Per-Sample Summary

| Metric | Description |
|--------|-------------|
| mean_depth | Average coverage |
| median_depth | Median coverage |
| min_depth | Minimum coverage |
| max_depth | Maximum coverage |
| pct_covered_10x | % of genome ≥10× |
| pct_covered_100x | % of genome ≥100× |
| n_lowcov_windows | Count of low coverage windows |
| n_highvar_windows | Count of high variance windows |

### Quality Thresholds

| Metric | Recommended | Minimum |
|--------|-------------|---------|
| Mean depth | >500× | >100× |
| % covered ≥10× | >99% | >95% |
| % covered ≥100× | >95% | >80% |
| Low coverage windows | <5 | <20 |

## Integration with Downstream

Coverage analysis informs:

1. **Consensus calling**: Low coverage → N bases
2. **Variant calling**: Depth affects confidence
3. **Homopolymer masking**: Flags unreliable regions

The BED files from this step are used in variant filtering to:
- Exclude variants in low-coverage regions
- Flag variants in homopolymer regions
- Identify multi-mapping artifacts

