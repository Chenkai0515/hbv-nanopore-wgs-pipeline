# Examples

This directory contains example configuration files and a mini dataset for testing the pipeline.

## Contents

```
examples/
├── README.md                    # This file
├── example_config.yaml          # Complete configuration example
├── samplesheet.csv              # Example sample sheet
├── mini_dataset/                # Subsampled test data (8 samples)
│   └── README.md                # Instructions for obtaining test data
└── expected_outputs/            # Reference outputs for validation
    └── README.md                # Expected output descriptions
```

## Quick Start

### 1. Set Up Configuration

```bash
# Copy example config to project config
cp examples/example_config.yaml config/config.yaml

# Edit with your paths
vim config/config.yaml
```

### 2. Prepare Sample Sheet

```bash
# Copy example sample sheet
cp examples/samplesheet.csv config/samplesheet.csv

# Edit with your samples
vim config/samplesheet.csv
```

### 3. Run Pipeline

```bash
# Single sample test
./run_pipeline.sh config/config.yaml --sample 10090

# Batch processing
./run_pipeline.sh config/config.yaml --batch config/samplesheet.csv
```

## Mini Dataset

The `mini_dataset/` directory is designed to hold subsampled FASTQ files for testing. 

### Test Samples

| Sample ID | Description |
|-----------|-------------|
| 10090 | Control sample |
| 10893 | Treatment sample with variants |
| 21024 | Standard sample |
| 23146 | Standard sample |
| 33117 | Standard sample |
| 40154 | Validation sample with variants |
| 40625 | Sample with variants |
| 40750 | Sample with variants |

### Obtaining Test Data

See `mini_dataset/README.md` for instructions on:
- Downloading pre-subsampled data
- Creating your own subsampled dataset
- Expected file sizes and read counts

## Example Configuration

The `example_config.yaml` file provides a complete configuration template with:

- All path settings
- Default parameter values
- Comments explaining each option

Key sections:
- `project`: Working directories
- `reference`: Reference genomes
- `databases`: Kraken2, Clair3 models
- `resources`: Thread counts, memory
- `qc/host_deconv/mapping/variants`: Module-specific parameters

## Sample Sheet Format

The `samplesheet.csv` file format:

```csv
sample_id,fastq_path,barcode,genotype,notes
10090,/path/to/10090.fastq.gz,barcode01,B,Control
10893,/path/to/10893.fastq.gz,barcode02,C,
```

Required columns:
- `sample_id`: Unique sample identifier
- `fastq_path`: Path to input FASTQ file

Optional columns:
- `barcode`: Sequencing barcode
- `genotype`: HBV genotype (A-H)
- `notes`: Any additional notes

## Expected Outputs

The `expected_outputs/` directory contains reference results for validation:

- QC statistics
- Coverage metrics
- Variant counts

Use these to verify your installation and identify any issues.

## Troubleshooting

### Config file not found
```bash
# Ensure you're in the pipeline root directory
cd /path/to/hbv-nanopore-wgs-pipeline
ls config/config.yaml
```

### Sample not found
```bash
# Verify sample ID matches samplesheet
grep "SAMPLE_ID" config/samplesheet.csv

# Check FASTQ path exists
ls -la /path/to/sample.fastq.gz
```

### Memory issues
```bash
# Reduce thread count in config
resources:
  threads: 8  # Reduce from default 16
```

## Contact

For questions about examples or test data, please open a GitHub Issue.

