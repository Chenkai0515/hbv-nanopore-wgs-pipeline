# Consensus Sequence Generation

## Overview

This module generates high-accuracy consensus sequences using Medaka's neural network-based approach with a two-round polishing strategy.

## Why Two Rounds?

### The Reference Bias Problem

Single-round consensus (reads → reference → consensus) suffers from reference bias:

```
Reference:  ...ATCGATCG...
True Sample: ...ATCAATCG...  (A→A at position 4)
                  ^
Biased Result: ...ATCGATCG...  (Reference allele favored)
```

**Cause:** Reads with the true variant align less well to the reference, biasing the pileup.

### Two-Round Solution

```
Round 1: Reads → Unified Reference → R1 Consensus
         ─────────────────────────────────────────
         Creates sample-specific consensus
         Some reference bias remains

Round 2: Reads → R1 Consensus → R2 Consensus
         ─────────────────────────────────────────
         Reads now align to sample-like reference
         Removes reference bias
         Captures true variants
```

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                 MEDAKA TWO-ROUND CONSENSUS                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   Step 9.1: Round 1 Alignment                                   │
│   ─────────────────────────────                                 │
│   Viral Reads ──► minimap2 ──► Unified Reference                │
│                       │                                         │
│                       ▼                                         │
│   Output: r1/minimap2_align/{sample}.bam                        │
│                       │                                         │
│                       ▼                                         │
│   Step 9.2: Round 1 Medaka                                      │
│   ─────────────────────────                                     │
│   BAM ──► medaka inference ──► medaka sequence                  │
│                       │                                         │
│                       ▼                                         │
│   Output: r1/{sample}/consensus.fasta                           │
│                       │                                         │
│                       ▼                                         │
│   Step 9.3: Round 2 Polish                                      │
│   ─────────────────────────                                     │
│   Reads ──► R1 Consensus ──► medaka_consensus                   │
│                       │                                         │
│                       ▼                                         │
│   Output: r2/{sample}/consensus.fasta  ← FINAL                  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Scripts

### Step 9.1: Round 1 Alignment

```bash
conda activate hbv_medaka

export FASTQ_ROOT="/path/to/host_deconv_out"
export REF_FASTA="/path/to/unified_reference.fasta"
export OUT_ROOT="/path/to/Medaka_consensus/r1"

bash 9.1_run_r1_minimap2_and_qc_v2.sh
```

**Output:**
- `r1/minimap2_align/bam/{sample}.primary.sorted.bam`
- `r1/minimap2_align/stats/{sample}.stats.tsv`

### Step 9.2: Round 1 Medaka

```bash
bash 9.2_run_medaka_from_bam_v2.sh
```

**Process:**
1. `medaka inference` - Neural network prediction
2. `medaka sequence` - Generate consensus FASTA

**Output:**
- `r1/{sample}/consensus.fasta`

### Step 9.3: Round 2 Polish

```bash
bash 9.3_run_medaka_round2_v2.sh
```

**Process:**
1. Align reads to R1 consensus
2. Run full `medaka_consensus` pipeline

**Output:**
- `r2/{sample}/consensus.fasta` ← **FINAL CONSENSUS**
- `r2/{sample}/consensus.fasta.fai`

## Medaka Models

### Model Selection

Choose model based on sequencing chemistry:

| Chemistry | Model |
|-----------|-------|
| R10.4.1 SUP | `r1041_e82_400bps_sup_v500` |
| R10.4.1 HAC | `r1041_e82_400bps_hac_v500` |
| R9.4.1 SUP | `r941_sup_plant_variant_v100` |
| R9.4.1 HAC | `r941_min_hac_g507` |

### Model Path

Models are typically in:
```bash
# Conda installation
$CONDA_PREFIX/lib/python3.x/site-packages/medaka/data/

# Manual installation
/path/to/medaka_models/
```

## Consensus Comparison (Step 10)

After generating R2 consensus, compare to unified reference:

```bash
python 10_consensus_comparison_with_transform.py \
    --consensus-dir Medaka_consensus_8/r2 \
    --reference unified_reference.fasta \
    --output-dir consensus_ref_viewa_9 \
    --offset 1823
```

### View A Output

The comparison produces **View A**: consensus vs reference differences.

| Column | Description |
|--------|-------------|
| position | Position in rotated reference |
| ref_base | Reference allele |
| consensus_base | Consensus allele |
| original_position | Position in original (unrotated) reference |

**View A captures:**
- Fixed genotype differences
- Consensus-level mutations
- Used to complement View B (read variants)

## Output Structure

```
Medaka_consensus_8/
├── r1/
│   ├── minimap2_align/
│   │   ├── indices/
│   │   ├── bam/
│   │   │   └── {sample}.primary.sorted.bam
│   │   ├── stats/
│   │   └── logs/
│   └── {sample}/
│       ├── consensus.fasta
│       └── medaka.log
└── r2/
    └── {sample}/
        ├── consensus.fasta      ← FINAL
        ├── consensus.fasta.fai
        └── medaka.log

consensus_ref_viewa_9/
├── {sample}_comparison_transformed.csv
├── integrated_summary.csv
└── variant_position_distribution.txt
```

## Quality Assessment

### Consensus Quality Indicators

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| Coverage | >500× | >100× | <50× |
| N bases | 0 | <10 | >50 |
| Length | ~3221 bp | ±50 bp | >100 bp diff |

### Common Issues

**N bases in consensus:**
- Cause: Low coverage regions
- Solution: Increase sequencing depth

**Consensus shorter than expected:**
- Cause: Deletion or gap
- Solution: Verify with read alignments

**Many differences from reference:**
- Cause: Different genotype (expected)
- Cause: Sample mix-up (investigate)

## Best Practices

1. **Verify input coverage** before running Medaka
   - Minimum 100× recommended
   - Optimal 500-1000×

2. **Check model compatibility**
   - Model must match sequencing chemistry
   - Wrong model = poor consensus

3. **Always use R2 consensus** for downstream
   - R1 is intermediate only
   - R2 has reduced reference bias

4. **Index all consensus files**
   ```bash
   samtools faidx consensus.fasta
   ```

## Troubleshooting

### Medaka fails with memory error

**Solution:** Reduce batch size
```bash
medaka inference --batch_size 100 ...
```

### Consensus quality poor

**Causes:**
- Low input coverage
- Wrong model selected
- Poor read quality

**Solutions:**
- Check coverage with depth analysis
- Verify model matches chemistry
- Review QC reports

### Coordinate mismatch with reference

**Cause:** Insertions/deletions shift coordinates

**Solution:** Use coordinate transformation scripts (11.6-11.8)

