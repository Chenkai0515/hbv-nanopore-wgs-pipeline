# 04_variants - Variant Calling Module

This module performs variant detection using a dual-tool approach (iVar + Clair3), applies tiered quality filtering, and transforms coordinates to a unified reference system.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        Variant Calling Pipeline                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Viral Enriched FASTQ + R2 Consensus                                       │
│       │                                                                     │
│       ▼                                                                     │
│   ┌───────────────────────────────────────────────────────────────┐        │
│   │ 11.1 Mapping to Per-Sample Consensus                          │        │
│   │      Input: reads → individual consensus                      │        │
│   │      Output: mapping/{sample}.bam                             │        │
│   └─────────────────────────────┬─────────────────────────────────┘        │
│                                 │                                           │
│                                 ▼                                           │
│   ┌───────────────────────────────────────────────────────────────┐        │
│   │ 11.2 BAM Filtering                                            │        │
│   │      Remove: low-qual, soft-clipped, short alignments         │        │
│   │      Output: filtering/{sample}.filtered.bam                  │        │
│   └─────────────────────────────┬─────────────────────────────────┘        │
│                                 │                                           │
│                  ┌──────────────┴──────────────┐                           │
│                  │                             │                            │
│                  ▼                             ▼                            │
│   ┌─────────────────────────┐   ┌─────────────────────────┐               │
│   │ 11.3a iVar Variants     │   │ 11.3b Clair3 Variants   │               │
│   │  (Amplicon-aware)       │   │  (Deep Learning)        │               │
│   │  Output: ivar/*.tsv     │   │  Output: clair3/*.vcf   │               │
│   └───────────┬─────────────┘   └───────────┬─────────────┘               │
│               │                             │                              │
│               └──────────────┬──────────────┘                              │
│                              │                                             │
│                              ▼                                             │
│   ┌───────────────────────────────────────────────────────────────┐        │
│   │ 11.4 Joint Filtering & Tiering                                │        │
│   │      - Merge iVar + Clair3 calls                              │        │
│   │      - Apply profile filters (strict/balanced/lenient)        │        │
│   │      - Assign quality tiers (HC/MC/LF)                        │        │
│   │      - Build cohort blacklist (PoN)                           │        │
│   │      Output: filtered/{sample}/                               │        │
│   └─────────────────────────────┬─────────────────────────────────┘        │
│                                 │                                           │
│                                 ▼                                           │
│   ┌───────────────────────────────────────────────────────────────┐        │
│   │ 11.5 Summarize Mutations                                      │        │
│   │      - Generate accepted.tsv (final calls)                    │        │
│   │      - Generate candidates.tsv (for validation)               │        │
│   │      Output: filtered/{sample}/accepted.tsv                   │        │
│   └─────────────────────────────┬─────────────────────────────────┘        │
│                                 │                                           │
│                                 ▼                                           │
│   ┌───────────────────────────────────────────────────────────────┐        │
│   │ 11.6-11.8 Coordinate Transformation                           │        │
│   │      - Transform: consensus coords → rotated ref coords       │        │
│   │      - Merge with View A (consensus vs reference)             │        │
│   │      - Flip to unified reference coordinates                  │        │
│   │      Output: unified/{sample}_unified_variants.csv            │        │
│   └───────────────────────────────────────────────────────────────┘        │
│                                                                             │
│   FINAL OUTPUT: variants/unified/{sample}_unified_variants.csv             │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Dual-Tool Strategy

### Why iVar + Clair3?

| Tool | Strengths | Weaknesses |
|------|-----------|------------|
| **iVar** | Amplicon-aware, handles strand bias, good for low-freq variants | May miss complex variants |
| **Clair3** | Deep learning, handles indels well, less reference bias | May overcall in low complexity regions |

By combining both:
- **DUAL** calls (both tools agree) = Highest confidence
- **Single-tool** calls = Reviewed with tier-specific thresholds

### Variant Sources

```
┌─────────────────────────────────────────────────────────────┐
│                   Variant Source Categories                 │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│   DUAL         │ Called by BOTH iVar and Clair3            │
│   ─────────────┼────────────────────────────────────────── │
│                │ Highest confidence, preferred for         │
│                │ clinical reporting                         │
│                │                                            │
│   IV_ONLY      │ Called by iVar only                       │
│   ─────────────┼────────────────────────────────────────── │
│                │ Subject to additional filtering            │
│                │ Often true variants missed by Clair3       │
│                │                                            │
│   CLAIR3_ONLY  │ Called by Clair3 only                     │
│   ─────────────┼────────────────────────────────────────── │
│                │ May include complex variants               │
│                │ Review for homopolymer artifacts           │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Quality Tiers

| Tier | AF Range | Criteria | Interpretation |
|------|----------|----------|----------------|
| **HC** | ≥20% | High coverage, strict QC, DUAL preferred | Report with confidence |
| **MC** | 5-20% | Good coverage, standard QC | Standard analysis |
| **LF** | 1-5% | Adequate coverage | Candidates for validation |

## Scripts

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `11_run_variants_call_10.sh` | — | — | **Main orchestrator** (runs all steps) |
| `11.1_mapping_batch_individual.py` | FASTQ + Consensus | BAM | Map reads to per-sample consensus |
| `11.2combined_bam_filter_batch.py` | BAM | Filtered BAM | Quality filtering |
| `11.3a_variant_call_ivar_batch.py` | BAM | TSV | iVar variant calling |
| `11.3b_batch_clair3_variants.sh` | BAM + Consensus | VCF | Clair3 variant calling |
| `11.4_ivar_clair3_filter_batch.py` | TSV + VCF | Filtered TSV | Joint filtering & tiering |
| `11.5_summarize_mutations_batch.py` | Filtered TSV | accepted/candidates | Final variant summary |
| `11.6_coordinate_transform_batch.py` | TSV | Transformed TSV | Consensus → rotated ref coords |
| `11.7_merge_variants_batch.py` | TSV + View A | Merged TSV | Combine View A + View B |
| `11.8_flip_to_unified_batch.py` | Merged TSV | Unified TSV | Final coordinate transformation |

## Usage

### Run Complete Pipeline

```bash
conda activate hbv_variants

# Run all steps
./11_run_variants_call_10.sh

# Run with options
./11_run_variants_call_10.sh --test 10090      # Single sample test
./11_run_variants_call_10.sh --from 4          # Start from step 4
./11_run_variants_call_10.sh --step 5          # Run only step 5
./11_run_variants_call_10.sh -y                # Skip confirmations
```

### Run Individual Steps

#### 11.1 Mapping

```bash
python 11.1_mapping_batch_individual.py \
    --fastq-dir /path/to/host_deconv_out \
    --consensus-dir /path/to/Medaka_consensus/r2 \
    --output-dir /path/to/variants_call/mapping \
    --threads 16
```

#### 11.2 BAM Filtering

```bash
python 11.2combined_bam_filter_batch.py \
    --input-dir /path/to/variants_call/mapping \
    --output-dir /path/to/variants_call/filtering \
    --min-mapq 20 \
    --max-soft-clip-frac 0.3 \
    --min-aligned-len 500
```

#### 11.3a iVar Calling

```bash
python 11.3a_variant_call_ivar_batch.py \
    --bam-dir /path/to/variants_call/filtering \
    --output-dir /path/to/variants_call/variants/ivar \
    --min-quality 20 \
    --min-frequency 0.01 \
    --min-depth 10
```

#### 11.3b Clair3 Calling

```bash
# Activate Clair3 environment
conda activate hbv_clair3

bash 11.3b_batch_clair3_variants.sh \
    --all \
    --threads 16 \
    -y
```

**Clair3 Experimental Parameters:**
```bash
# For viral genomes, enable these experimental flags:
--var_pct_full=1.0          # Process full pileup
--ref_pct_full=1.0          # Full reference coverage
--haploid_sensitive         # Haploid-aware mode
--enable_variant_calling_at_sequence_head_and_tail
```

#### 11.4 Joint Filtering

```bash
python 11.4_ivar_clair3_filter_batch.py \
    --all \
    --profile lenient \
    -y
```

**Filter Profiles:**
- `strict`: Highest specificity (clinical reporting)
- `balanced`: Standard analysis (recommended)
- `lenient`: Maximum sensitivity (exploratory)

#### 11.5 Summarize Mutations

```bash
python 11.5_summarize_mutations_batch.py \
    --input-dir /path/to/variants_call/variants/filtered \
    --output-dir /path/to/variants_call/variants/filtered \
    --all
```

#### 11.6-11.8 Coordinate Transformation

```bash
# Transform coordinates
python 11.6_coordinate_transform_batch.py \
    --input-dir /path/to/variants_call/variants/filtered \
    --output-dir /path/to/variants_call/variants/transformed \
    --offset 1823

# Merge with View A
python 11.7_merge_variants_batch.py \
    --variants-dir /path/to/variants_call/variants/transformed \
    --viewa-dir /path/to/consensus_ref_viewa \
    --output-dir /path/to/variants_call/variants/combined

# Flip to unified reference
python 11.8_flip_to_unified_batch.py \
    --input-dir /path/to/variants_call/variants/combined \
    --output-dir /path/to/variants_call/variants/unified
```

## View A + View B System

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Coordinate Views Explained                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   VIEW A: Consensus vs Unified Reference                                    │
│   ──────────────────────────────────────                                    │
│   What: Differences between sample consensus and reference                  │
│   Source: Step 10 (10_consensus_comparison_with_transform.py)               │
│   Use: Identify fixed differences (genotype-specific, consensus-level)      │
│                                                                             │
│   VIEW B: Reads vs Consensus                                                │
│   ──────────────────────────────────────                                    │
│   What: Variants detected in reads aligned to consensus                     │
│   Source: Steps 11.3-11.5 (iVar + Clair3)                                  │
│   Use: Identify within-sample variation (quasispecies, mutations)          │
│                                                                             │
│   UNIFIED: Combined View                                                    │
│   ──────────────────────────────────────                                    │
│   What: All variants transformed to unified reference coordinates           │
│   Source: Steps 11.6-11.8                                                  │
│   Use: Cross-sample comparison, cohort analysis                            │
│                                                                             │
│   Transformation Chain:                                                     │
│   ────────────────────                                                     │
│   Reads → Consensus (View B coords)                                         │
│                ↓                                                            │
│   Consensus → Rotated Reference (transformed)                               │
│                ↓                                                            │
│   + View A differences → Unified Reference (final)                          │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Blacklist (Panel of Normals)

The pipeline builds a cohort-level blacklist to filter recurrent artifacts:

```
Blacklist Construction:
────────────────────────
1. Collect all variants across cohort
2. Count samples with each variant
3. Flag variants appearing in ≥30% of samples at AF <5%
4. These are likely:
   - Systematic sequencing errors
   - Reference bias artifacts
   - Homopolymer-associated errors

Application:
────────────────────────
- Blacklisted variants marked but NOT removed
- Allows manual review of flagged calls
- Can be overridden for known mutations
```

## Output Directory Structure

```
variants_call_10/
├── mapping/
│   └── {sample}.bam
├── filtering/
│   └── {sample}.filtered.bam
├── variants/
│   ├── ivar/
│   │   └── {sample}/
│   │       └── variants.tsv
│   ├── clair3/
│   │   └── {sample}_clair3/
│   │       └── merge_output.vcf.gz
│   ├── filtered/
│   │   ├── {sample}/
│   │   │   ├── ivar_only_final.tsv
│   │   │   ├── accepted.tsv        ← Accepted variants
│   │   │   └── candidates.tsv      ← Candidates for validation
│   │   ├── blacklist.tsv           ← Cohort-level artifacts
│   │   └── summary_counts.tsv
│   ├── transformed/
│   │   └── {sample}_transformed.csv
│   ├── combined/
│   │   └── {sample}_combined.csv
│   └── unified/
│       └── {sample}_unified_variants.csv  ← FINAL OUTPUT
└── logs/
    └── *.log
```

## Final Output Format

`{sample}_unified_variants.csv` columns:

| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| unified_pos | Position in unified reference |
| ref | Reference allele |
| alt | Alternate allele |
| af | Allele frequency |
| depth | Total depth |
| source | DUAL / IV_ONLY / CLAIR3_ONLY |
| tier | HC / MC / LF |
| decision | ACCEPT / CANDIDATE / REJECT |
| blacklist | TRUE if flagged as artifact |
| orig_consensus_pos | Original position in consensus |
| view_a_diff | TRUE if also in View A |

## Troubleshooting

### No variants detected
- Check coverage (need ≥30× for most variants)
- Verify consensus quality
- Sample may have no within-host variation (common for acute infection)

### Many blacklisted variants
- Normal for cohort-level artifacts
- Review blacklist.tsv for patterns
- Homopolymer regions commonly flagged

### Clair3 fails
- Check TensorFlow/Keras versions (need 2.x, not 3.x)
- On macOS ARM: ensure GNU getopt is installed
- Verify model path matches sequencing chemistry

### Coordinate mismatch
- Verify offset parameter (1823 for this reference)
- Check consensus naming matches expected pattern
- Review transformation logs for errors

## Dependencies

| Tool | Version | Environment |
|------|---------|-------------|
| minimap2 | ≥2.24 | hbv_base |
| samtools | ≥1.15 | hbv_base |
| ivar | ≥1.3.1 | hbv_variants |
| Clair3 | ≥1.0.0 | hbv_clair3 |
| bcftools | ≥1.15 | hbv_clair3 |
| pandas | ≥1.4 | hbv_variants |
| pysam | ≥0.19 | hbv_variants |

