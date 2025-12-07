# Host Decontamination Module

## Overview

The host decontamination module separates HBV viral reads from human host sequences. This is critical because:

1. Clinical samples contain both viral and human DNA
2. HBV can integrate into the human genome
3. False host assignments can occur without proper handling

## Methodology

### The HBV Masking Approach

Standard host removal (map to human, discard mapped reads) fails for HBV because:
- HBV integration sites exist in the human genome
- Reads spanning integrations map to both references
- Without masking, viral reads are incorrectly classified as host

**Solution: Mask-then-map strategy**

```
┌─────────────────────────────────────────────────────────────────┐
│                    MASK-THEN-MAP STRATEGY                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Original Read:                                                 │
│  5'──────[HBV Sequence]──────3'                                │
│                                                                 │
│  Step 1: Map to HBV panel                                       │
│          Identify HBV-aligned regions                           │
│                                                                 │
│  Masked Read:                                                   │
│  5'──────[NNNNNNNNNNNN]──────3'                                │
│           (HBV → N's)                                           │
│                                                                 │
│  Step 2: Map masked read to human                               │
│          Only non-HBV portions can match                        │
│                                                                 │
│  Result: Correct classification                                 │
│          - True host reads: map to human                        │
│          - True viral reads: don't map (masked regions)         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Handling Chimeric Reads

Integration events create chimeric reads (human-HBV junctions):

```
Chimeric Read Example:
5'──[Human Sequence]──[HBV Sequence]──[Human Tail]──3'
      (150 bp)           (2800 bp)       (250 bp)

After Masking:
5'──[Human Sequence]──[NNNNNNNNNNNNN]──[Human Tail]──3'

After Host Mapping:
- Human portions align to GRCh38
- Masked HBV portion doesn't align

Decision:
- If tails are small (<300 bp, <20% of read): TRIM and keep HBV
- If internal human: DISCARD (complex chimera)
```

## Decision Rules

### Read Classification Flowchart

```
                         Input Read
                              │
                              ▼
                    ┌─────────────────┐
                    │ Map to HBV Panel│
                    └────────┬────────┘
                             │
              ┌──────────────┴──────────────┐
              │                             │
              ▼                             ▼
        HBV Alignment                 No HBV Alignment
              │                             │
              ▼                             │
        Mask HBV Regions                    │
              │                             │
              ▼                             ▼
        Map Masked Read ◄───────────────────┘
        to Human Genome
              │
    ┌─────────┴─────────────┐
    │                       │
    ▼                       ▼
Human Mapped            No Human Hit
    │                       │
    │                       ▼
    │               ┌───────────────┐
    │               │ VIRAL_ENRICHED│
    │               └───────────────┘
    │
    ▼
┌───────────────────────────────┐
│ Check Human Alignment Pattern │
└───────────┬───────────────────┘
            │
  ┌─────────┴─────────┐
  │                   │
  ▼                   ▼
Hits at Ends      Hits Throughout
  │                   │
  ▼                   ▼
┌─────────┐     ┌───────────┐
│TRIM TAILS│    │ HOST_ONLY │
└────┬────┘     └───────────┘
     │
     ▼
Length After Trim ≥ MIN_KEEP?
     │
   ┌─┴─┐
   │   │
   ▼   ▼
  Yes  No
   │   │
   ▼   ▼
VIRAL  DISCARD
ENRICHED
```

### Parameter Thresholds

#### HBV Masking

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `HBV_MIN_MAPQ` | 20 | Confident alignment |
| `HBV_MIN_MATCH_LEN` | 100 | Avoid spurious short matches |
| `HBV_MERGE_DISTANCE` | 10 | Merge nearby HBV regions |

#### Human Classification

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `HUMAN_MIN_MAPQ` | 20 | Confident alignment |
| `HUMAN_MIN_MATCH_LEN` | 100 | Avoid spurious matches |
| `HUMAN_MIN_FRAC` | 0.05 | At least 5% of read aligns |
| `HUMAN_MIN_PID` | 0.80 | 80% identity minimum |
| `HUMAN_MIN_TOTAL_ALN` | 200 | At least 200 bp total alignment |

#### Tail Trimming

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `TAIL_NEAR_END` | 80 | Consider "near end" if within 80 bp |
| `TAIL_MAX_LEN` | 300 | Maximum tail to trim |
| `TAIL_MAX_FRAC` | 0.20 | Maximum 20% of read as tail |
| `MIN_KEEP_LEN` | 800 | Keep reads ≥800 bp after trim |
| `DISCARD_INTERNAL` | 1 | Discard reads with internal host |

## Scripts

### Main Pipeline (`6.2_host_deconv_mask_hbv.v2.sh`)

Orchestrates the complete host decontamination workflow.

```bash
# Set environment variables
export HBV_PANEL="/path/to/hbv_panel.fasta"
export RAW_DIR="/path/to/porechop_output"
export OUTDIR="/path/to/host_deconv_output"
export HUMAN_REF_FASTA="/path/to/GRCh38.primary_assembly.fa.gz"
export THREADS=16

# Run
bash 6.2_host_deconv_mask_hbv.v2.sh
```

### Helper Scripts

| Script | Purpose |
|--------|---------|
| `6.1_mask_hbv_regions.py` | Replace HBV-aligned bases with N |
| `6.3_trim_host_tails.py` | Remove host-derived read ends |
| `6.4_paf_partition_ids.py` | Classify reads as host/viral |

## Output Files

```
host_deconv_out/
├── {sample}/
│   ├── {sample}.vsHBV.paf.gz           # HBV alignments
│   ├── {sample}.maskHBV.fastq.gz       # Masked reads
│   ├── {sample}.maskHBV.vsHuman.paf.gz # Human alignments
│   ├── host.ids.txt                     # Host read IDs
│   ├── {sample}.host_only.fastq.gz     # Host reads (masked)
│   ├── {sample}.viral_enriched.fastq.gz     # Viral (masked)
│   └── {sample}.viral_enriched.unmasked.fastq.gz  # ← USE THIS
├── summary.tsv                          # Per-sample stats
└── logs/
```

### Summary Statistics

The `summary.tsv` file contains:

| Column | Description |
|--------|-------------|
| total_reads | Input read count |
| hbv_hit_reads | Reads with HBV alignment |
| masked_reads | Reads with masked bases |
| masked_bases | Total bases masked |
| host_reads | Reads classified as host |
| viral_enriched_reads | Reads kept as viral |
| trimmed_reads | Reads with tails trimmed |
| trimmed_bp | Total bases trimmed |

## Expected Results

### Typical Distribution

For a successful HBV extraction:

| Category | Expected % | Notes |
|----------|------------|-------|
| Viral enriched | 40-80% | Target reads |
| Host only | 10-40% | Normal carryover |
| Trimmed (kept) | 5-15% | Chimeric reads |
| Discarded | <5% | Complex chimeras |

### Warning Signs

| Observation | Possible Issue |
|-------------|----------------|
| >50% host | Extraction failure |
| <20% viral | Low viral load |
| >20% trimmed | Integration-heavy sample |

## Troubleshooting

### High host fraction

**Causes:**
- Inefficient viral enrichment
- Low viral load in sample
- Extraction protocol issues

**Solutions:**
- Review extraction QC
- Consider deeper sequencing
- Adjust enrichment protocol

### Many reads discarded

**Causes:**
- Complex integration events
- `MIN_KEEP_LEN` too high
- Chimeric library

**Solutions:**
- Lower `MIN_KEEP_LEN` if short fragments acceptable
- Review integration patterns
- Consider paired-end or targeted sequencing

## Integration Analysis

The host decontamination outputs enable integration analysis:

1. **PAF files** contain precise alignment coordinates
2. **Trimmed reads** indicate junction positions
3. **Host-mapped regions** reveal integration sites

For dedicated integration analysis, examine:
- `{sample}.orig.vsHuman.paf.gz` - Original reads vs human
- Reads that have both HBV and Human hits
- Junction positions in trimmed reads

