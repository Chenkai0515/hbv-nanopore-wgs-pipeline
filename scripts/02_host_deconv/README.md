# 02_host_deconv - Host Decontamination Module

This module separates HBV viral reads from human host sequences using a multi-step approach: HBV masking, human mapping, and tail trimming.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      Host Decontamination Pipeline                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Clean FASTQ (from QC)                                                     │
│       │                                                                     │
│       ▼                                                                     │
│   ┌───────────────────┐                                                     │
│   │ 6.1 Mask HBV      │  Mask HBV regions to prevent false host calls       │
│   │     Regions       │  Output: *.maskHBV.fastq.gz                         │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐                                                     │
│   │ 6.2 Host Deconv   │  Map masked reads to human reference                │
│   │     (Main Script) │  Partition: host_only vs viral_enriched             │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│       ┌─────┴─────┐                                                         │
│       │           │                                                         │
│       ▼           ▼                                                         │
│   host_only   viral_enriched                                                │
│       │           │                                                         │
│       │           ▼                                                         │
│       │   ┌───────────────────┐                                             │
│       │   │ 6.3 Trim Host     │  Remove host-derived tails                  │
│       │   │     Tails         │  from chimeric reads                        │
│       │   └─────────┬─────────┘                                             │
│       │             │                                                       │
│       ▼             ▼                                                       │
│   ┌─────────┐   ┌─────────────────┐                                        │
│   │ Host    │   │ Viral Enriched  │                                        │
│   │ Reads   │   │ (HBV + minimal  │                                        │
│   │         │   │  host tails)    │                                        │
│   └─────────┘   └─────────────────┘                                        │
│                                                                             │
│   Output files per sample:                                                  │
│   - *.host_only.fastq.gz        (pure host reads)                          │
│   - *.viral_enriched.fastq.gz   (HBV reads, masked version)                │
│   - *.viral_enriched.unmasked.fastq.gz (HBV reads, original sequence)      │
│   - *.host_only.unmasked.fastq.gz (host reads, original)                   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Design Rationale

### Why Mask HBV First?

HBV integrates into the human genome at various sites. Without masking:
1. Integrated HBV sequences in human ref would capture viral reads
2. HBV reads mapping to integration sites would be incorrectly classified as "host"

By masking HBV regions before human mapping, we ensure clean separation.

### Why Trim Host Tails?

Clinical samples often contain chimeric reads (human-HBV junctions):
```
5'──[Human Sequence]──[HBV Sequence]──[Human Tail]──3'
```

These arise from:
- Integration breakpoints
- Library prep artifacts
- Ligation chimeras

Tail trimming preserves the HBV portion while removing flanking host sequences.

## Scripts

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `6.1_mask_hbv_regions.py` | FASTQ + PAF | Masked FASTQ | Replace HBV-aligned bases with 'N' |
| `6.2_host_deconv_mask_hbv.v2.sh` | Clean FASTQ | Partitioned FASTQs | Main orchestration script |
| `6.3_trim_host_tails.py` | FASTQ + PAF | Trimmed FASTQ | Remove host-derived read tails |
| `6.4_paf_partition_ids.py` | PAF | Read ID lists | Classify reads as host/viral |

## Usage

### Main Script (Recommended)

```bash
conda activate hbv_base

# Set environment variables
export HBV_PANEL="/path/to/hbv_panel.fasta"
export RAW_DIR="/path/to/porechop_output"
export OUTDIR="/path/to/host_deconv_output"
export HUMAN_REF_FASTA="/path/to/GRCh38.primary_assembly.fa.gz"
export THREADS=16

# Run host deconvolution
bash 6.2_host_deconv_mask_hbv.v2.sh
```

### Individual Scripts

#### 6.1 Mask HBV Regions

```bash
python 6.1_mask_hbv_regions.py \
    --paf sample.vsHBV.paf.gz \
    --fastq sample.fastq.gz \
    --out sample.maskHBV.fastq.gz \
    --min-mapq 20 \
    --merge-distance 10 \
    --min-match-len 100
```

#### 6.3 Trim Host Tails

```bash
python 6.3_trim_host_tails.py \
    --paf sample.vsHuman.paf.gz \
    --fastq sample.viral_enriched.fastq.gz \
    --out sample.trimmed.fastq.gz \
    --near-end 80 \
    --max-tail-len 300 \
    --max-tail-frac 0.20 \
    --min-keep-len 800
```

#### 6.4 Partition Read IDs

```bash
python 6.4_paf_partition_ids.py \
    --paf sample.vsHuman.paf.gz \
    --out-host-ids host.ids.txt \
    --min-mapq 20 \
    --min-match-len 100 \
    --min-frac 0.05 \
    --min-pid 0.80
```

## Key Parameters

### HBV Masking Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-mapq` | 20 | Minimum MAPQ for HBV alignment |
| `--merge-distance` | 10 | Merge HBV hits within N bp |
| `--min-match-len` | 100 | Minimum alignment match length |

### Host Classification Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HUMAN_MIN_MAPQ` | 20 | Minimum MAPQ for human classification |
| `HUMAN_MIN_MATCH_LEN` | 100 | Minimum match length |
| `HUMAN_MIN_FRAC` | 0.05 | Min fraction of read aligned to human |
| `HUMAN_MIN_PID` | 0.80 | Minimum percent identity |
| `HUMAN_MIN_TOTAL_ALN` | 200 | Minimum total aligned bases |
| `HUMAN_MIN_TOTAL_FRAC` | 0.05 | Min fraction of read with human hits |

### Tail Trimming Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TAIL_NEAR_END` | 80 | Consider alignment "near end" if within N bp |
| `TAIL_MAX_LEN` | 300 | Maximum tail length to trim |
| `TAIL_MAX_FRAC` | 0.20 | Maximum tail as fraction of read |
| `MIN_KEEP_LEN` | 800 | Minimum remaining read length |
| `DISCARD_INTERNAL` | 1 | Discard reads with internal host |

## Decision Rules Flowchart

```
                    Read Classification Decision Tree
                    ─────────────────────────────────
                                    │
                                    ▼
                        ┌───────────────────┐
                        │ Align to Human    │
                        │ (masked reads)    │
                        └─────────┬─────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    │                           │
                    ▼                           ▼
            Human hits ≥ threshold?       No human hits
                    │                           │
                    ▼                           ▼
        ┌───────────────────┐           ┌───────────────┐
        │ Check hit pattern │           │ VIRAL_ENRICHED│
        └─────────┬─────────┘           └───────────────┘
                  │
        ┌─────────┴─────────┐
        │                   │
        ▼                   ▼
    Hits at ends?      Hits throughout?
        │                   │
        ▼                   ▼
    ┌─────────┐       ┌───────────┐
    │ TRIM    │       │ HOST_ONLY │
    │ TAILS   │       └───────────┘
    └────┬────┘
         │
         ▼
    Remaining ≥ MIN_KEEP_LEN?
         │
    ┌────┴────┐
    │         │
    ▼         ▼
   Yes        No
    │         │
    ▼         ▼
┌─────────┐ ┌─────────┐
│ VIRAL_  │ │ DISCARD │
│ ENRICHED│ │ (too    │
│ (trimmed)│ │ short)  │
└─────────┘ └─────────┘
```

## Output Directory Structure

```
host_deconv_out/
├── ref/
│   ├── hbv_panel.mmi         # HBV minimap2 index
│   └── human.mmi             # Human minimap2 index
├── logs/
│   └── run_YYYYMMDD_HHMMSS.log
├── summary.tsv               # Per-sample statistics
└── {sample}_subsampled.trimmed_filtered.porechop/
    ├── {sample}.vsHBV.paf.gz              # HBV alignments
    ├── {sample}.maskHBV.fastq.gz          # Masked reads
    ├── {sample}.maskHBV.vsHuman.paf.gz    # Human alignments (masked)
    ├── {sample}.orig.vsHuman.paf.gz       # Human alignments (original)
    ├── host.ids.txt                        # Host read IDs
    ├── {sample}.host_only.fastq.gz        # Pure host reads (masked)
    ├── {sample}.host_only.unmasked.fastq.gz
    ├── {sample}.viral_enriched.fastq.gz   # Viral reads (masked)
    ├── {sample}.viral_enriched.unmasked.fastq.gz  # ← For downstream
    ├── {sample}.maskHBV.log
    └── {sample}.trimhost.log
```

## Summary Statistics

The `summary.tsv` file contains per-sample metrics:

| Column | Description |
|--------|-------------|
| sample | Sample ID |
| total_reads | Input read count |
| hbv_hit_reads | Reads with HBV alignments |
| masked_reads | Reads with masked regions |
| masked_bases | Total bases masked |
| host_reads | Reads classified as host |
| viral_enriched_reads | Reads classified as viral |
| viral_unmasked_reads | Viral reads (unmasked version) |
| host_unmasked_reads | Host reads (unmasked version) |
| trimmed_reads | Reads with tails trimmed |
| trimmed_bp | Total bases trimmed |

## Performance Notes

- **Memory**: ~16GB for human reference indexing
- **Time**: ~2-5 min per 10,000 reads (16 threads)
- **Disk**: Temporary PAF files can be large; ~2x input FASTQ size

## Troubleshooting

### High host fraction (>50%)
- Check extraction protocol (did viral enrichment work?)
- Review library prep (target capture efficiency)
- Lower `HUMAN_MIN_FRAC` if integration analysis is intended

### Low viral enrichment
- Verify HBV panel matches sample genotype
- Check original sample viral load
- Review QC metrics from previous step

### Many reads discarded after tail trimming
- Consider increasing `MAX_TAIL_LEN` for samples with longer integrations
- Lower `MIN_KEEP_LEN` if accepting shorter fragments is acceptable

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | ≥2.24 | Read alignment |
| samtools | ≥1.15 | Index operations |
| seqkit | ≥2.3.0 | FASTQ filtering |
| pigz | ≥2.6 | Parallel compression |
| Python | ≥3.9 | Helper scripts |

