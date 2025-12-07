# HBV Nanopore WGS Pipeline

A comprehensive bioinformatics pipeline for **Hepatitis B Virus (HBV)** whole-genome sequencing analysis using Oxford Nanopore long-read data.

## Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        HBV Nanopore WGS Pipeline                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  Raw FASTQ ──► QC & Filtering ──► Host Deconv ──► Mapping ──► Consensus    │
│                    │                  │              │            │         │
│               ┌────┴────┐        ┌────┴────┐    ┌───┴───┐   ┌────┴────┐   │
│               │ Dorado  │        │ Minimap2│    │Medaka │   │ Variant │   │
│               │ Porechop│        │ SeqKit  │    │ R1+R2 │   │ Calling │   │
│               │ FastQC  │        │ Masking │    │       │   │iVar+Clair3│ │
│               └─────────┘        └─────────┘    └───────┘   └─────────┘   │
│                                                                             │
│  Output: High-quality consensus sequences + Annotated variant calls        │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Features

- **Quality Control**: Adapter trimming (Dorado), length/quality filtering, mid-adapter removal (Porechop)
- **Host Decontamination**: HBV masking + human read separation + tail trimming
- **Reference Alignment**: Minimap2 mapping with coverage analysis and homopolymer masking
- **Consensus Generation**: Two-round Medaka polishing for high-accuracy consensus
- **Variant Calling**: Dual-tool approach (iVar + Clair3) with tiered confidence levels
- **Coordinate System**: Rotated reference (offset=1823) with unified coordinate transformation

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- Reference genomes (HBV panel + Human GRCh38)
- ~100GB disk space for databases (Kraken2, human reference)

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/hbv-nanopore-wgs-pipeline.git
cd hbv-nanopore-wgs-pipeline

# Create conda environments (choose modules as needed)
mamba env create -f environment/env_base.yml
mamba env create -f environment/env_qc.yml
mamba env create -f environment/env_medaka.yml
mamba env create -f environment/env_clair3.yml
mamba env create -f environment/env_variants.yml

# Copy and edit configuration
cp config/config.example.yaml config/config.yaml
# Edit config/config.yaml with your paths
```

## Quick Start

### Single Sample Processing

```bash
# 1. Edit your configuration
cp config/config.example.yaml config/config.yaml
vim config/config.yaml  # Set your paths

# 2. Run the complete pipeline
./run_pipeline.sh config/config.yaml --sample SAMPLE_ID

# Or run step by step:
# Step 1: QC
conda activate hbv_qc
bash scripts/01_qc/1_trim_adapters_dorado_fastq.sh --input raw_data/ --output trimmed/

# Step 2: Host decontamination
conda activate hbv_base
bash scripts/02_host_deconv/6.2_host_deconv_mask_hbv.v2.sh

# Step 3: Mapping & Consensus
conda activate hbv_medaka
bash scripts/03_mapping_consensus/9.1_run_r1_minimap2_and_qc_v2.sh
bash scripts/03_mapping_consensus/9.2_run_medaka_from_bam_v2.sh
bash scripts/03_mapping_consensus/9.3_run_medaka_round2_v2.sh

# Step 4: Variant calling
conda activate hbv_variants
bash scripts/04_variants/11_run_variants_call_10.sh
```

### Batch Processing

```bash
# Prepare sample sheet
cp config/samplesheet.example.csv config/samplesheet.csv
# Edit samplesheet.csv with your samples

# Run batch pipeline
./run_pipeline.sh config/config.yaml --batch config/samplesheet.csv
```

## Pipeline Modules

| Module | Scripts | Description |
|--------|---------|-------------|
| **01_qc** | `1_*.sh`, `2_*.py`, `3_*.sh`, `4_*.sh`, `5_*.sh` | Raw read QC, adapter trimming, filtering |
| **02_host_deconv** | `6.1_*.py`, `6.2_*.sh`, `6.3_*.py`, `6.4_*.py` | Host separation, HBV masking, tail trimming |
| **03_mapping_consensus** | `7_*.sh`, `8_*.sh`, `9.1-9.3_*.sh`, `10_*.py` | Alignment, coverage, Medaka consensus |
| **04_variants** | `11_*.sh`, `11.1-11.8_*.py/sh` | iVar+Clair3 calling, filtering, annotation |

See [docs/overview.md](docs/overview.md) for detailed workflow documentation.

## Output Structure

```
results/
├── qc/                    # FastQC, NanoPlot reports
├── host_deconv/           # Host-depleted reads
├── mapping/               # BAM files, coverage stats
├── consensus/             # Medaka consensus sequences
└── variants/
    ├── ivar/              # iVar raw calls
    ├── clair3/            # Clair3 raw calls
    ├── filtered/          # Filtered variants with tiers
    └── unified/           # Final unified-coordinate variants
```

## Variant Quality Tiers

| Tier | AF Range | Description |
|------|----------|-------------|
| **HC** (High Confidence) | ≥20% | High coverage, strict QC passed |
| **MC** (Medium Confidence) | 5-20% | Good coverage, standard QC |
| **LF** (Low Frequency) | 1-5% | Candidate variants for validation |

## Coordinate System

This pipeline uses a **rotated HBV reference** starting from DR1:
- **Offset**: 1823 bp
- **Genome Length**: 3221 bp
- **Transform**: `original_pos = ((current_pos - 1) + 1823) % 3221 + 1`

See [docs/coordinate_systems.md](docs/coordinate_systems.md) for details.

## Requirements

### Software Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Dorado | ≥1.1.0 | Adapter trimming |
| Porechop_ABI | latest | Mid-adapter removal |
| Minimap2 | ≥2.24 | Read alignment |
| Samtools | ≥1.15 | BAM manipulation |
| Medaka | ≥1.8.0 | Consensus polishing |
| Clair3 | ≥1.0.0 | Deep learning variant calling |
| iVar | ≥1.3 | Amplicon-aware variant calling |
| Kraken2 | ≥2.1 | Contamination screening |

### Reference Data

- HBV reference panel (genotypes A-H)
- Human reference genome (GRCh38)
- Kraken2 database (Standard or Standard-8)

## Citation

If you use this pipeline in your research, please cite:

```
Jiang, C. (2025). HBV Nanopore WGS Pipeline: A comprehensive workflow for 
hepatitis B virus whole-genome sequencing analysis. GitHub repository: 
https://github.com/YOUR_USERNAME/hbv-nanopore-wgs-pipeline
```

See [CITATION.cff](CITATION.cff) for machine-readable citation information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions and issues, please open a GitHub Issue or contact the author.

---

**Author**: Chenkai Jiang  
**Version**: 1.0.0  
**Last Updated**: December 2025

