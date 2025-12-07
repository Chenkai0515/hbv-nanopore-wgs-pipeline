# Variant Calling Module

## Overview

The variant calling module detects within-sample variation using a dual-tool approach (iVar + Clair3), applies tiered quality filtering, and transforms coordinates to a unified reference system.

## Strategy: Why Dual Tools?

### iVar Strengths
- **Amplicon-aware**: Handles strand bias from PCR
- **Low-frequency variants**: Designed for viral quasispecies
- **Simple output**: TSV format, easy to parse

### Clair3 Strengths
- **Deep learning**: Better indel detection
- **No reference bias**: Neural network approach
- **Complex variants**: Handles MNPs

### Combined Approach

```
                    Filtered BAM
                         │
         ┌───────────────┴───────────────┐
         │                               │
         ▼                               ▼
    ┌─────────┐                     ┌─────────┐
    │  iVar   │                     │ Clair3  │
    │  (TSV)  │                     │  (VCF)  │
    └────┬────┘                     └────┬────┘
         │                               │
         └───────────────┬───────────────┘
                         │
                         ▼
              ┌─────────────────────┐
              │ Joint Filtering     │
              │ - Merge calls       │
              │ - Assign sources    │
              │ - Apply tiers       │
              └─────────────────────┘
                         │
              ┌──────────┼──────────┐
              │          │          │
              ▼          ▼          ▼
           DUAL      IV_ONLY   CLAIR3_ONLY
         (highest)   (review)   (review)
```

## Quality Tiers

### Tier Definitions

| Tier | AF Range | Criteria | Use Case |
|------|----------|----------|----------|
| **HC** | ≥20% | High depth, strict QC, DUAL preferred | Clinical reporting |
| **MC** | 5-20% | Good depth, standard QC | Standard analysis |
| **LF** | 1-5% | Adequate depth, flag for validation | Exploratory |

### Tier Assignment Logic

```python
def assign_tier(af, depth, source, qc_pass):
    if af >= 0.20 and depth >= 100 and qc_pass:
        return "HC"
    elif af >= 0.05 and depth >= 50:
        return "MC"
    elif af >= 0.01 and depth >= 30:
        return "LF"
    else:
        return "REJECT"
```

## Filter Profiles

Three profiles for different use cases:

### Strict Profile
```yaml
# Maximum specificity
snp:
  min_af: 0.05
  min_depth: 50
  min_qual: 30
indel:
  min_af: 0.10
  min_depth: 100
  min_qual: 40
```

### Balanced Profile (Default)
```yaml
# Good sensitivity/specificity balance
snp:
  min_af: 0.03
  min_depth: 30
  min_qual: 20
indel:
  min_af: 0.05
  min_depth: 50
  min_qual: 25
```

### Lenient Profile
```yaml
# Maximum sensitivity
snp:
  min_af: 0.01
  min_depth: 10
  min_qual: 15
indel:
  min_af: 0.02
  min_depth: 20
  min_qual: 20
```

## Cohort Blacklist (Panel of Normals)

### Purpose

Filter recurrent low-frequency artifacts that appear across samples:
- Systematic sequencing errors
- Reference bias artifacts
- Homopolymer-associated errors

### Construction

```python
for variant in all_cohort_variants:
    sample_count = count_samples_with_variant(variant)
    max_af = max_af_across_samples(variant)
    
    if sample_count / total_samples >= 0.30 and max_af < 0.05:
        add_to_blacklist(variant)
```

### Application

Blacklisted variants are **flagged, not removed**:
- Allows manual review
- Can be overridden for known mutations
- Preserves full call set

## Workflow Scripts

### Step 11.1: Mapping to Per-Sample Consensus

```bash
python 11.1_mapping_batch_individual.py \
    --fastq-dir host_deconv_out \
    --consensus-dir Medaka_consensus/r2 \
    --output-dir variants_call/mapping
```

**Why per-sample consensus?**
- Removes reference bias
- Variants relative to sample's own sequence
- Better for low-frequency detection

### Step 11.2: BAM Filtering

```bash
python 11.2combined_bam_filter_batch.py \
    --input-dir variants_call/mapping \
    --output-dir variants_call/filtering \
    --min-mapq 20 \
    --max-soft-clip-frac 0.3
```

**Filters:**
- MAPQ ≥ 20
- Soft-clipping < 30% of read
- Aligned length ≥ 500 bp

### Step 11.3a: iVar Calling

```bash
python 11.3a_variant_call_ivar_batch.py \
    --bam-dir variants_call/filtering \
    --output-dir variants_call/variants/ivar \
    --min-quality 20 \
    --min-frequency 0.01
```

### Step 11.3b: Clair3 Calling

```bash
conda activate hbv_clair3

bash 11.3b_batch_clair3_variants.sh \
    --all \
    --threads 16
```

**Experimental Parameters for Viral Genomes:**
```bash
--var_pct_full=1.0          # Full pileup processing
--ref_pct_full=1.0          # Full reference coverage
--haploid_sensitive         # Haploid mode
--enable_variant_calling_at_sequence_head_and_tail
```

### Step 11.4: Joint Filtering

```bash
python 11.4_ivar_clair3_filter_batch.py \
    --all \
    --profile lenient
```

**Process:**
1. Load iVar TSV and Clair3 VCF
2. Merge overlapping calls
3. Assign sources (DUAL/IV_ONLY/CLAIR3_ONLY)
4. Apply profile thresholds
5. Assign quality tiers
6. Build cohort blacklist

### Step 11.5: Summarize Mutations

```bash
python 11.5_summarize_mutations_batch.py \
    --input-dir variants_call/variants/filtered \
    --all
```

**Output:**
- `accepted.tsv`: Final variant calls
- `candidates.tsv`: Variants for validation

### Steps 11.6-11.8: Coordinate Transformation

See [Coordinate Systems](coordinate_systems.md) for details.

## View A + View B System

### View A: Consensus vs Reference

```
Source: Step 10 (consensus comparison)
Content: Fixed differences between sample consensus and unified reference
Example: Genotype-specific SNPs, consensus-level variants
```

### View B: Reads vs Consensus

```
Source: Steps 11.3-11.5 (variant calling)
Content: Within-sample variation (quasispecies)
Example: Drug resistance mutations, mixed infections
```

### Unified Output

Combines both views in unified reference coordinates:

```
Final Variants = View B (transformed) + View A (reference differences)
                 └─────────────────────────────────────────────────┘
                            All in unified reference coords
```

## Output Files

```
variants_call_10/
├── mapping/
│   └── {sample}.bam
├── filtering/
│   └── {sample}.filtered.bam
├── variants/
│   ├── ivar/
│   │   └── {sample}/variants.tsv
│   ├── clair3/
│   │   └── {sample}_clair3/merge_output.vcf.gz
│   ├── filtered/
│   │   ├── {sample}/
│   │   │   ├── ivar_only_final.tsv
│   │   │   ├── accepted.tsv       ← Accepted variants
│   │   │   └── candidates.tsv     ← For validation
│   │   ├── blacklist.tsv          ← Cohort artifacts
│   │   └── summary_counts.tsv
│   └── unified/
│       └── {sample}_unified_variants.csv  ← FINAL
└── logs/
```

## Final Output Format

`{sample}_unified_variants.csv`:

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

## Special Sites

### BCP Mutations (positions 1762, 1764)

```
A1762T + G1764A: Associated with HCC risk
```

### PreCore Mutations (positions 1896, 1899)

```
G1896A: Creates stop codon, HBeAg-negative phenotype
```

### Drug Resistance Sites

RT domain positions for antiviral resistance:
- rtL169, rtV173, rtL180, rtA181, rtS184
- rtM204, rtN236, rtM250

## Troubleshooting

### No variants detected

**Possible causes:**
- Very uniform sample (no quasispecies)
- Coverage too low
- Consensus quality issues

**Solutions:**
- Check depth analysis
- Review consensus quality
- Consider acute vs chronic infection

### Many blacklisted variants

**Normal behavior** for:
- Homopolymer regions
- Low-complexity sequences
- Known systematic errors

**Action:**
- Review blacklist.tsv
- Override for known true variants

### Clair3 memory errors

**Solution:**
```bash
# Reduce parallel jobs
--threads 4

# Or increase memory
export CLAIR3_MAX_MEMORY=32
```

