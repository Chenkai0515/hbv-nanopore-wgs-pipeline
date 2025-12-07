# 01_qc - Quality Control Module

This module handles raw FASTQ processing from basecalling output to clean, filtered reads ready for downstream analysis.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         QC & Filtering Pipeline                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Raw FASTQ                                                                 │
│       │                                                                     │
│       ▼                                                                     │
│   ┌───────────────────┐                                                     │
│   │ 1. Dorado Trim    │  Adapter trimming (sequencing kit-specific)         │
│   │    (Step 1)       │  Output: *.trimmed.fastq.gz                         │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐                                                     │
│   │ 2. Length/Quality │  Filter: 2-4kb length, Q17+ quality                 │
│   │    Filter (Step 2)│  Output: *_filtered.fastq.gz                        │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐                                                     │
│   │ 3. Porechop       │  Mid-read adapter removal (chimeric reads)          │
│   │    (Step 3)       │  Output: *.porechop.fastq.gz                        │
│   └─────────┬─────────┘                                                     │
│             │                                                               │
│             ▼                                                               │
│   ┌───────────────────┐     ┌───────────────────┐                          │
│   │ 4. Multi-tool QC  │ ──► │ 5. Kraken2 QC     │  Contamination check      │
│   │    (Step 4)       │     │    (Step 5)       │  Output: *.kreport        │
│   └───────────────────┘     └───────────────────┘                          │
│   FastQC, NanoPlot,         Kraken2 + Krona                                │
│   nanoQC, SeqKit            visualization                                  │
│                                                                             │
│             ▼                                                               │
│   Clean FASTQ (ready for host decontamination)                              │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Scripts

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `1_trim_adapters_dorado_fastq.sh` | Raw FASTQ | `*.trimmed.fastq.gz` | Dorado adapter trimming with kit-specific settings |
| `2_filter_fastq_gz_2-4k_Q17.py` | Trimmed FASTQ | `*_filtered.fastq.gz` | Length (2-4kb) and quality (Q17+) filtering |
| `3_run_porechop_midcheck.sh` | Filtered FASTQ | `*.porechop.fastq.gz` | Mid-read adapter removal using Porechop_ABI |
| `4_ont_multi_qc.sh` | Porechop FASTQ | QC reports | Comprehensive QC: FastQC, NanoPlot, nanoQC, SeqKit |
| `5_kraken2_nanopore_qc.sh` | Porechop FASTQ | Kraken reports | Contamination screening with Kraken2 + Krona |

## Usage Examples

### Step 1: Adapter Trimming (Dorado)

```bash
# Activate dorado environment
conda activate dorado

# Test mode (first file, limited reads)
bash 1_trim_adapters_dorado_fastq.sh --test 2000 \
    --input /path/to/raw_fastq \
    --output /path/to/trimmed \
    --kit SQK-NBD114-24

# Full run
bash 1_trim_adapters_dorado_fastq.sh --run \
    --input /path/to/raw_fastq \
    --output /path/to/trimmed \
    --kit SQK-NBD114-24 \
    --jobs 4
```

**Key Parameters:**
- `--kit`: Sequencing kit name (e.g., SQK-NBD114-24)
- `--jobs`: Number of parallel jobs
- `--threads-per-job`: CPU threads per job

### Step 2: Length and Quality Filtering

```bash
conda activate hbv_base

python 2_filter_fastq_gz_2-4k_Q17.py \
    --input /path/to/trimmed \
    --output /path/to/filtered \
    --min-len 2000 \
    --max-len 4000 \
    --min-quality 17
```

**Key Parameters:**
- `--min-len`: Minimum read length (default: 2000)
- `--max-len`: Maximum read length (default: 4000)
- `--min-quality`: Minimum mean Q score (default: 17)

**Rationale:** HBV genome is ~3.2kb. Reads 2-4kb capture full-length plus some flanking. Q17 (~2% error rate) balances sensitivity with quality.

### Step 3: Porechop Mid-Adapter Removal

```bash
conda activate hbv_qc

bash 3_run_porechop_midcheck.sh \
    --input /path/to/filtered \
    --output /path/to/porechop \
    --threads 8
```

**Why Porechop?** Nanopore reads occasionally contain mid-read adapters from chimeric molecules. Porechop_ABI detects and splits these, preventing alignment artifacts.

### Step 4: Multi-Tool QC

```bash
conda activate hbv_qc

# Set paths
export IN_DIR="/path/to/porechop"
export OUT_DIR="/path/to/qc_results"

bash 4_ont_multi_qc.sh
```

**Generated Reports:**
- `fastqc/`: Per-base quality, sequence length distribution
- `nanoplot/`: ONT-specific metrics (N50, read length vs quality)
- `nanoqc/`: Per-position quality heatmaps
- `seqkit/`: Summary statistics (TSV format)
- `multiqc/`: Aggregated HTML report

### Step 5: Kraken2 Contamination Check

```bash
conda activate hbv_qc

# Ensure Kraken2 database is set up
# export KRAKEN2_DB=/path/to/kraken2_db

bash 5_kraken2_nanopore_qc.sh \
    --input /path/to/porechop \
    --output /path/to/kraken_results \
    --db /path/to/kraken2_db
```

**Expected Results:**
- Most reads: Unclassified (HBV not in standard Kraken DB)
- Some reads: Human (host carryover)
- Flag: Bacterial/reagent taxa (see `config/kraken_reagent_taxa.tsv`)

## Output Directory Structure

```
multi_tool_qc_4/
├── fastqc/
│   ├── sample1_fastqc.html
│   └── sample1_fastqc.zip
├── nanoplot/
│   └── sample1/
│       ├── NanoPlot-report.html
│       ├── LengthvsQualityScatterPlot_dot.png
│       └── ...
├── nanoqc/
│   └── sample1_nanoQC.html
├── seqkit/
│   └── stats.tsv
├── kraken_raw/
│   └── sample1.kraken
├── reports/
│   └── sample1.kreport
├── krona/
│   └── sample1.html
├── multiqc/
│   └── multiqc_report.html
├── logs/
│   └── *.log
└── versions.txt
```

## Quality Thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| Read length | 2000-4000 bp | HBV genome ~3.2kb |
| Mean Q score | ≥17 | ~2% error rate |
| N50 | ≥3000 bp | Full-length coverage |
| Adapter content | <1% | Clean trimming |

## Troubleshooting

### Low read count after filtering
- Check raw data quality with NanoPlot first
- Consider relaxing Q threshold to 15 for exploratory analysis
- Verify sequencing kit parameter matches actual kit used

### High adapter content in FastQC
- Ensure correct `--kit` parameter in Dorado
- Re-run Porechop with `--discard_middle` for strict chimera removal

### Kraken showing unexpected taxa
- Compare against `config/kraken_reagent_taxa.tsv`
- Low-abundance reagent contaminants (<0.1%) are typically acceptable
- High human fraction (>50%) indicates extraction issues

## Dependencies

| Tool | Version | Environment |
|------|---------|-------------|
| Dorado | ≥1.1.0 | dorado |
| FastQC | ≥0.12.1 | hbv_qc |
| NanoPlot | ≥1.40.0 | hbv_qc |
| nanoQC | ≥0.9.4 | hbv_qc |
| SeqKit | ≥2.3.0 | hbv_qc |
| MultiQC | ≥1.14 | hbv_qc |
| Porechop_ABI | ≥0.5.0 | hbv_qc |
| Kraken2 | ≥2.1.2 | hbv_qc |
| Krona | ≥2.8.1 | hbv_qc |

