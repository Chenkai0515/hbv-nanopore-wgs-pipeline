# Pipeline Overview

## Introduction

The HBV Nanopore WGS Pipeline is designed to process Oxford Nanopore long-read sequencing data from Hepatitis B Virus (HBV) samples. It provides a complete workflow from raw basecalled reads to annotated variant calls.

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         HBV NANOPORE WGS PIPELINE                               │
│                              Complete Workflow                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  ╔═══════════════════════════════════════════════════════════════════════════╗ │
│  ║                        PHASE 1: QUALITY CONTROL                           ║ │
│  ╠═══════════════════════════════════════════════════════════════════════════╣ │
│  ║                                                                           ║ │
│  ║   Raw FASTQ ──► Dorado Trim ──► Length/Q Filter ──► Porechop              ║ │
│  ║       │             │                │                  │                 ║ │
│  ║       │             │                │                  │                 ║ │
│  ║       └─────────────┴────────────────┴──────────────────┘                 ║ │
│  ║                              │                                            ║ │
│  ║                              ▼                                            ║ │
│  ║                    FastQC / NanoPlot / Kraken2                            ║ │
│  ║                    (Quality & Contamination Reports)                      ║ │
│  ║                                                                           ║ │
│  ╚═══════════════════════════════════════════════════════════════════════════╝ │
│                                      │                                          │
│                                      ▼                                          │
│  ╔═══════════════════════════════════════════════════════════════════════════╗ │
│  ║                    PHASE 2: HOST DECONTAMINATION                          ║ │
│  ╠═══════════════════════════════════════════════════════════════════════════╣ │
│  ║                                                                           ║ │
│  ║   Clean FASTQ ──► Mask HBV ──► Map to Human ──► Partition Reads           ║ │
│  ║                                      │                                    ║ │
│  ║                            ┌─────────┴─────────┐                          ║ │
│  ║                            ▼                   ▼                          ║ │
│  ║                      Host Reads          Viral Enriched                   ║ │
│  ║                      (discard)                 │                          ║ │
│  ║                                               │                          ║ │
│  ║                                    Trim Host Tails                        ║ │
│  ║                                               │                          ║ │
│  ║                                               ▼                          ║ │
│  ║                                    Viral Reads (clean)                    ║ │
│  ║                                                                           ║ │
│  ╚═══════════════════════════════════════════════════════════════════════════╝ │
│                                      │                                          │
│                                      ▼                                          │
│  ╔═══════════════════════════════════════════════════════════════════════════╗ │
│  ║                   PHASE 3: MAPPING & CONSENSUS                            ║ │
│  ╠═══════════════════════════════════════════════════════════════════════════╣ │
│  ║                                                                           ║ │
│  ║   Viral Reads ──► minimap2 ──► Unified HBV Reference                      ║ │
│  ║                       │                                                   ║ │
│  ║                       ▼                                                   ║ │
│  ║               Coverage Analysis                                           ║ │
│  ║               (mosdepth, homopolymer masks)                               ║ │
│  ║                       │                                                   ║ │
│  ║                       ▼                                                   ║ │
│  ║   ┌─────────────────────────────────────────────────────────────────┐    ║ │
│  ║   │              MEDAKA TWO-ROUND CONSENSUS                         │    ║ │
│  ║   │                                                                 │    ║ │
│  ║   │   Reads ──► Unified Ref ──► R1 Consensus                        │    ║ │
│  ║   │                                   │                             │    ║ │
│  ║   │   Reads ──► R1 Consensus ──► R2 Consensus (FINAL)              │    ║ │
│  ║   │                                                                 │    ║ │
│  ║   └─────────────────────────────────────────────────────────────────┘    ║ │
│  ║                       │                                                   ║ │
│  ║                       ▼                                                   ║ │
│  ║   Consensus vs Reference Comparison (View A)                              ║ │
│  ║                                                                           ║ │
│  ╚═══════════════════════════════════════════════════════════════════════════╝ │
│                                      │                                          │
│                                      ▼                                          │
│  ╔═══════════════════════════════════════════════════════════════════════════╗ │
│  ║                      PHASE 4: VARIANT CALLING                             ║ │
│  ╠═══════════════════════════════════════════════════════════════════════════╣ │
│  ║                                                                           ║ │
│  ║   Reads ──► Per-Sample Consensus ──► Filtered BAM                         ║ │
│  ║                                           │                               ║ │
│  ║                        ┌──────────────────┴──────────────────┐            ║ │
│  ║                        │                                     │            ║ │
│  ║                        ▼                                     ▼            ║ │
│  ║                   ┌─────────┐                           ┌─────────┐       ║ │
│  ║                   │  iVar   │                           │ Clair3  │       ║ │
│  ║                   │ (TSV)   │                           │ (VCF)   │       ║ │
│  ║                   └────┬────┘                           └────┬────┘       ║ │
│  ║                        │                                     │            ║ │
│  ║                        └──────────────┬──────────────────────┘            ║ │
│  ║                                       │                                   ║ │
│  ║                                       ▼                                   ║ │
│  ║                        Joint Filtering & Tiering                          ║ │
│  ║                        (HC / MC / LF tiers)                               ║ │
│  ║                                       │                                   ║ │
│  ║                                       ▼                                   ║ │
│  ║                        Cohort Blacklist (PoN)                             ║ │
│  ║                                       │                                   ║ │
│  ║                                       ▼                                   ║ │
│  ║                        Coordinate Transformation                          ║ │
│  ║                        (View B + View A → Unified)                        ║ │
│  ║                                                                           ║ │
│  ╚═══════════════════════════════════════════════════════════════════════════╝ │
│                                      │                                          │
│                                      ▼                                          │
│  ┌───────────────────────────────────────────────────────────────────────────┐ │
│  │                           FINAL OUTPUTS                                   │ │
│  │                                                                           │ │
│  │   • QC Reports (HTML): FastQC, NanoPlot, MultiQC, Krona                  │ │
│  │   • Consensus Sequences: {sample}/consensus.fasta                         │ │
│  │   • Coverage Analysis: depth BEDs, summary TSVs                           │ │
│  │   • Variant Calls: {sample}_unified_variants.csv                          │ │
│  │   • Cohort Summary: blacklist.tsv, summary_counts.tsv                     │ │
│  │                                                                           │ │
│  └───────────────────────────────────────────────────────────────────────────┘ │
│                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
```

## Design Principles

### 1. Modular Architecture

Each phase is self-contained and can be run independently:
- Facilitates debugging and development
- Allows reprocessing from any stage
- Supports incremental analysis

### 2. Dual Variant Calling

Using both iVar and Clair3 maximizes both sensitivity and specificity:
- iVar: Optimized for amplicon data, handles strand bias
- Clair3: Deep learning approach, good for complex variants

### 3. Tiered Quality Assessment

Variants are categorized by confidence level:
- HC (High Confidence): Clinical reporting quality
- MC (Medium Confidence): Standard analysis
- LF (Low Frequency): Requires validation

### 4. Coordinate Standardization

All outputs use a unified coordinate system:
- Enables cross-sample comparison
- Facilitates cohort-level analysis
- Consistent with standard HBV annotation

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Raw reads | FASTQ/FASTQ.GZ | Basecalled ONT reads |
| HBV reference | FASTA | Unified rotated reference |
| Human reference | FASTA.GZ | GRCh38 for host decontamination |
| Kraken2 DB | k2d | Standard or Standard-8 database |
| Clair3 models | — | ONT-specific model files |

## Output Summary

| Output | Location | Description |
|--------|----------|-------------|
| QC Reports | `results/qc/` | HTML reports from multiple tools |
| Host-depleted reads | `results/host_deconv/` | Viral-enriched FASTQ files |
| Alignments | `results/mapping/` | BAM files and stats |
| Consensus | `results/consensus/` | Final FASTA sequences |
| Variants | `results/variants/unified/` | Annotated variant calls |

## Resource Requirements

### Computational

| Phase | CPU | Memory | Time (8 samples) |
|-------|-----|--------|------------------|
| QC | 16 cores | 8 GB | ~30 min |
| Host deconv | 16 cores | 32 GB | ~1 hour |
| Mapping/Consensus | 16 cores | 16 GB | ~2 hours |
| Variant calling | 16 cores | 16 GB | ~1 hour |

### Storage

| Data Type | Size (per sample) |
|-----------|-------------------|
| Raw FASTQ | ~1-5 GB |
| Intermediate BAM | ~2-10 GB |
| Final outputs | ~50-100 MB |

## Conda Environments

The pipeline uses separate environments for different tools:

```bash
# Base tools (minimap2, samtools, seqkit)
mamba env create -f environment/env_base.yml

# QC tools (FastQC, NanoPlot, Kraken2)
mamba env create -f environment/env_qc.yml

# Medaka consensus
mamba env create -f environment/env_medaka.yml

# Clair3 variant calling
mamba env create -f environment/env_clair3.yml

# iVar and analysis
mamba env create -f environment/env_variants.yml
```

## Configuration

All paths and parameters are configured in `config/config.yaml`:

```yaml
project:
  work_dir: "/path/to/project"
  
reference:
  hbv_unified: "/path/to/ref/rotated_reference.fasta"
  human_ref: "/path/to/ref/GRCh38.fa.gz"
  
resources:
  threads: 16
  java_memory_gb: 8
```

See `config/config.example.yaml` for all available options.

## Running the Pipeline

### Single Sample

```bash
./run_pipeline.sh config/config.yaml --sample SAMPLE_ID
```

### Batch Processing

```bash
# Prepare sample sheet
cp config/samplesheet.example.csv config/samplesheet.csv
# Edit with your samples

# Run batch
./run_pipeline.sh config/config.yaml --batch config/samplesheet.csv
```

### Step-by-Step

```bash
# Run specific step
./run_pipeline.sh config/config.yaml --step 3

# Resume from step
./run_pipeline.sh config/config.yaml --from 4
```

## Next Steps

- [QC Module](qc.md): Detailed QC documentation
- [Host Decontamination](host_deconv.md): Host removal methodology
- [Consensus](consensus.md): Medaka two-round approach
- [Variants](variants.md): Dual-tool calling strategy
- [Coordinate Systems](coordinate_systems.md): Reference transformation

