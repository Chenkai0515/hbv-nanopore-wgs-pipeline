# Mini Dataset

This directory contains a subsampled HBV sequencing dataset for testing the pipeline.

## Contents

```
mini_dataset/
├── raw_fastq/
│   └── 10090_subsampled.fastq    # ~2000 reads, ~7MB
├── work/                          # Pipeline outputs (created on first run)
└── README.md                      # This file
```

## Sample Information

| Field | Value |
|-------|-------|
| Sample ID | 10090 |
| Reads | ~2000 |
| Mean Read Length | ~3.2 kb |
| Mean Quality | Q18+ |
| Size | ~7 MB |
| Format | FASTQ (uncompressed) |

## Origin

This sample was subsampled from a real HBV whole-genome sequencing run performed on an Oxford Nanopore MinION device using R10.4.1 flow cells with the SUP basecalling model.

The subsampling was performed to:
- Reduce file size for quick testing
- Maintain realistic data characteristics
- Enable pipeline validation in minutes rather than hours

## Usage

```bash
# Run pipeline on this sample
./run_pipeline.sh examples/example_config.yaml --sample 10090

# Dry-run to see what would execute
./run_pipeline.sh examples/example_config.yaml --sample 10090 --dry-run
```

## Expected Results

After running the full pipeline, this sample should produce:
- High-coverage consensus sequence (~3221 bp)
- QC reports showing good quality reads
- Variant calls (if any variants present)
