# 03_mapping_consensus - Mapping & Consensus Module

This module handles reference alignment, coverage analysis, and consensus sequence generation using Medaka's two-round polishing approach.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    Mapping & Consensus Pipeline                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Viral Enriched FASTQ                                                      │
│       │                                                                     │
│       ▼                                                                     │
│   ┌───────────────────┐                                                     │
│   │ 7. Mapping        │  minimap2 to unified HBV reference                  │
│   │    (TA1)          │  Output: sorted BAM + stats                         │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐                                                     │
│   │ 8. Depth Analysis │  mosdepth coverage + homopolymer masks              │
│   │    (TA2)          │  Output: per-base depth, low-cov windows            │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────────────────────────────────────────────┐            │
│   │              Medaka Two-Round Consensus                    │            │
│   │  ┌─────────────────────────────────────────────────────┐  │            │
│   │  │ 9.1 Round 1: Initial alignment + QC                 │  │            │
│   │  │     Input: viral_enriched → unified_ref             │  │            │
│   │  │     Output: r1/minimap2_align/                      │  │            │
│   │  └───────────────────────┬─────────────────────────────┘  │            │
│   │                          │                                 │            │
│   │                          ▼                                 │            │
│   │  ┌─────────────────────────────────────────────────────┐  │            │
│   │  │ 9.2 Round 1: Medaka consensus from BAM              │  │            │
│   │  │     Input: BAM from 9.1                             │  │            │
│   │  │     Output: r1/{sample}/consensus.fasta             │  │            │
│   │  └───────────────────────┬─────────────────────────────┘  │            │
│   │                          │                                 │            │
│   │                          ▼                                 │            │
│   │  ┌─────────────────────────────────────────────────────┐  │            │
│   │  │ 9.3 Round 2: Re-polish against R1 consensus         │  │            │
│   │  │     Input: reads + r1/consensus.fasta               │  │            │
│   │  │     Output: r2/{sample}/consensus.fasta (FINAL)     │  │            │
│   │  └─────────────────────────────────────────────────────┘  │            │
│   └───────────────────────────────────────────────────────────┘            │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐                                                     │
│   │ 10. Consensus     │  Align R2 consensus to unified reference            │
│   │     Comparison    │  Output: View A (consensus vs reference)            │
│   └───────────────────┘                                                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Coordinate System

This pipeline uses a **rotated HBV reference** to handle the circular genome:

```
Original HBV Genome (circular, 3221 bp):
        ┌───────────────────────────────────────┐
        │                                       │
        │         PreS1    PreS2    S           │
     5' ├──┬──────┬────────┬───────┬──┤ 3'     │
        │  │      │        │       │  │         │
        │  Core   │        │   Pol │  │         │
        │  1-1838 │        │       │  │         │
        │         │        │       │  │         │
        └─────────┴────────┴───────┴──┴─────────┘
              ▲
              │ Original position 1

Rotated Reference (starts at DR1, position 1824):
     5' ├──────────────────────────────────────┤ 3'
        │ DR1 ─────────────────────────── Core │
        │ (new pos 1)                (new pos) │
        └──────────────────────────────────────┘

Offset: 1823 bp
Conversion: original_pos = ((current_pos - 1) + 1823) % 3221 + 1
```

**Why rotate?**
- DR1 (Direct Repeat 1) is a natural breakpoint for linearization
- Avoids splitting important genes across the reference boundary
- Facilitates integration site analysis

## Scripts

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `7_TA1_map_V2_unmaskreads.sh` | Viral FASTQ | Sorted BAM | minimap2 alignment with secondary analysis |
| `8_TA2_depth_v2.sh` | BAM | Depth BED/TSV | Coverage analysis, homopolymer masks |
| `9.1_run_r1_minimap2_and_qc_v2.sh` | Viral FASTQ | R1 BAM | Round 1 alignment + QC stats |
| `9.2_run_medaka_from_bam_v2.sh` | R1 BAM | R1 consensus | Medaka consensus from BAM |
| `9.3_run_medaka_round2_v2.sh` | Reads + R1 | R2 consensus | Round 2 polish (final) |
| `10_consensus_comparison_with_transform.py` | R2 consensus | CSV comparison | View A: consensus vs reference |

## Usage

### Step 7: Mapping to Unified Reference

```bash
conda activate hbv_base

# Set paths
export BAMS_DIR="/path/to/output/bam"
export REF="/path/to/reordered_1824_3221_reference.fasta"
export FASTQ_DIR="/path/to/host_deconv_out"

bash 7_TA1_map_V2_unmaskreads.sh
```

**Output:**
- `bam/*.primary.sorted.bam` - Primary alignments
- `bam/*.all.sorted.bam` - All alignments (including secondary)
- `stats/*.flagstat.txt` - Alignment statistics
- `hotspots/*.tsv` - Multi-mapping regions

### Step 8: Depth Analysis

```bash
conda activate hbv_variants

# Set paths
export BAMS_DIR="/path/to/TA1_map/bam"
export REF="/path/to/reordered_reference.fasta"
export OUT_DIR="/path/to/depth_output"

bash 8_TA2_depth_v2.sh
```

**Output:**
- `mosdepth/*.mosdepth.summary.txt` - Coverage summary
- `perbase/*.perbase.bed.gz` - Per-base depth
- `windows/*.windows.bed.gz` - Sliding window coverage
- `flags/*.lowcov.bed` - Low coverage regions
- `flags/*.highvar.bed` - High variance regions
- `homopolymers/*.hpoly_ge4.bed` - Homopolymer masks

### Step 9.1-9.3: Medaka Two-Round Consensus

```bash
conda activate hbv_medaka

# Round 1: Alignment + QC
export FASTQ_ROOT="/path/to/host_deconv_out"
export REF_FASTA="/path/to/reordered_reference.fasta"
export OUT_ROOT="/path/to/Medaka_consensus/r1"

bash 9.1_run_r1_minimap2_and_qc_v2.sh

# Round 1: Medaka consensus
bash 9.2_run_medaka_from_bam_v2.sh

# Round 2: Polish against R1 consensus
bash 9.3_run_medaka_round2_v2.sh
```

**Output:**
- `r1/{sample}/consensus.fasta` - Round 1 consensus
- `r2/{sample}/consensus.fasta` - **Final consensus** (use this!)
- `r2/{sample}/consensus.fasta.fai` - FASTA index

### Step 10: Consensus Comparison

```bash
conda activate hbv_base

python 10_consensus_comparison_with_transform.py \
    --consensus-dir /path/to/Medaka_consensus/r2 \
    --reference /path/to/reordered_reference.fasta \
    --output-dir /path/to/consensus_comparison \
    --offset 1823 \
    --genome-length 3221
```

**Output (View A):**
- `{sample}_comparison_transformed.csv` - Position-by-position differences
- `integrated_summary.csv` - All samples summary
- `variant_position_distribution.txt` - Hotspot analysis

## Key Parameters

### Mapping Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MM2_PRESET` | map-ont | minimap2 preset |
| `MIN_MAPQ` | 20 | Minimum mapping quality |
| `THREADS` | 16 | CPU threads |

### Depth Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WIN` | 200 | Sliding window size (bp) |
| `STEP` | 50 | Window step size (bp) |
| `LOW_FRAC` | 0.20 | Low coverage = < 20% of median |
| `LOW_ABS` | 100 | Absolute low coverage threshold |
| `MAD_MULT` | 4 | MAD multiplier for high variance |
| `HPOLY_MIN4/6` | 4/6 | Homopolymer run lengths to mask |
| `HPOLY_PAD` | 5 | Padding around homopolymers |

### Medaka Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Model | r1041_e82_400bps_sup_v500 | ONT R10.4.1 SUP model |
| Min depth | 10 | Minimum depth for consensus |

## Why Two-Round Consensus?

```
Round 1: Reads → Unified Reference → R1 Consensus
         ─────────────────────────────────────
         High-error reads aligned to reference.
         Reference bias may affect consensus.

Round 2: Reads → R1 Consensus → R2 Consensus (Final)
         ─────────────────────────────────────────────
         Reads now align to sample-specific consensus.
         Removes reference bias, captures true variants.
```

This approach:
1. **Reduces reference bias** - Sample-specific consensus captures true sequence
2. **Improves accuracy** - Re-alignment to closer reference improves base calls
3. **Handles genotype diversity** - Works across HBV genotypes A-H

## Output Directory Structure

```
Medaka_consensus_8/
├── r1/
│   ├── minimap2_align/
│   │   ├── indices/
│   │   ├── bam/
│   │   │   └── {sample}.primary.sorted.bam
│   │   ├── stats/
│   │   │   └── {sample}.stats.tsv
│   │   └── logs/
│   └── {sample}/
│       ├── consensus.fasta
│       └── medaka.log
└── r2/
    └── {sample}/
        ├── consensus.fasta      ← FINAL CONSENSUS
        ├── consensus.fasta.fai
        └── medaka.log

TA2_depth_7/
├── mosdepth/
├── perbase/
├── windows/
├── flags/
├── homopolymers/
├── summary/
│   ├── per_sample_summary.tsv
│   └── all_samples_combined.tsv
└── logs/

consensus_ref_viewa_9/
├── {sample}_comparison_transformed.csv
├── integrated_summary.csv
└── variant_position_distribution.txt
```

## Coverage QC Interpretation

### Low Coverage Windows

Regions with depth < max(LOW_ABS, LOW_FRAC × median) are flagged:

| Cause | Action |
|-------|--------|
| True deletion | Verify with reads; may be real variant |
| PCR dropout | Check primer design |
| GC bias | Common in extreme GC regions |
| Integration site | Expected at human-HBV junctions |

### High Variance Windows

Regions with |coverage - median| > MAD_MULT × MAD:

| Cause | Action |
|-------|--------|
| Amplification hotspot | May indicate PCR bias |
| Repeated region | Check for multi-mapping |
| Copy number variation | Biological signal |

### Homopolymer Regions

ONT has systematic errors in homopolymers:
- 4+ bp runs: Flag for review
- 6+ bp runs: High error rate, use caution
- Padded regions (±5bp): Buffer for alignment uncertainty

## Troubleshooting

### Low consensus quality
- Check input read coverage (need ≥100× for good consensus)
- Verify Medaka model matches sequencing chemistry
- Review R1 alignment stats before R2

### Missing positions in consensus
- Low coverage regions produce 'N' bases
- Check depth analysis output for gaps
- May need more sequencing depth

### Consensus differs significantly from reference
- Expected for different genotypes
- Verify sample genotype assignment
- Check for sample contamination

## Dependencies

| Tool | Version | Environment |
|------|---------|-------------|
| minimap2 | ≥2.24 | hbv_base |
| samtools | ≥1.15 | hbv_base |
| medaka | ≥1.8.0 | hbv_medaka |
| mosdepth | ≥0.3.3 | hbv_variants |
| Python | ≥3.9 | hbv_base |
| BioPython | ≥1.79 | hbv_base |

