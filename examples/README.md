# Example Data and Configuration

This directory contains a minimal example dataset to demonstrate the pipeline.

## Quick Start (for Reviewers)

You can test the pipeline in **~10-20 minutes** using the included mini dataset:

```bash
# 1. Check the pipeline without running (dry-run)
./run_pipeline.sh examples/example_config.yaml --sample 10090 --dry-run

# 2. Run QC steps only (steps 1-4, no external databases needed)
./run_pipeline.sh examples/example_config.yaml --sample 10090 --from 1 --to 4
```

## Contents

```
examples/
├── example_config.yaml     # Pre-configured config for mini dataset
├── samplesheet.csv         # Sample sheet for batch processing
├── README.md               # This file
├── mini_dataset/
│   ├── raw_fastq/          # Sample FASTQ (10090_subsampled.fastq)
│   └── README.md           # Dataset description
└── expected_outputs/
    ├── qc_stats.tsv        # Expected QC statistics
    ├── consensus.fasta     # Expected consensus sequence (truncated)
    └── variants.csv        # Expected variant calls
```

## Mini Dataset

The mini dataset includes a **subsampled** HBV sequencing run:

| Sample | Reads | Size | Description |
|--------|-------|------|-------------|
| 10090  | ~2000 | ~7MB | Subsampled from full ONT run |

This sample has been subsampled to enable quick testing while preserving realistic data characteristics.

## Running the Pipeline

### Prerequisites

1. **Create conda environments** (one-time setup):
   ```bash
   make env-qc        # For steps 1-5
   make env-medaka    # For steps 9-10
   make env-variants  # For step 11
   ```

2. **Edit configuration** (update reference paths):
   ```bash
   # Edit examples/example_config.yaml
   # Update paths marked with EDIT_THIS_PATH
   ```

### Run Options

```bash
# Dry-run (show what would execute)
./run_pipeline.sh examples/example_config.yaml --sample 10090 --dry-run

# Run QC only (no external references needed)
./run_pipeline.sh examples/example_config.yaml --sample 10090 --from 1 --to 4

# Run full pipeline (requires all references)
./run_pipeline.sh examples/example_config.yaml --sample 10090

# Run specific step
./run_pipeline.sh examples/example_config.yaml --sample 10090 --step 2
```

### Expected Runtime

| Steps | Description | Time (4 threads) |
|-------|-------------|------------------|
| 1-4   | QC pipeline | ~5 minutes |
| 5     | Kraken2 | ~2 minutes |
| 6     | Host deconv | ~3 minutes |
| 7-8   | Mapping | ~2 minutes |
| 9     | Medaka | ~10 minutes |
| 10    | Comparison | ~1 minute |
| 11    | Variants | ~5 minutes |

**Total**: ~30 minutes with all steps

## Expected Outputs

After running the full pipeline, you should see:

```
examples/mini_dataset/work/
├── fastq_dorado/               # Step 1: Trimmed FASTQ
├── fastq_filter_2/             # Step 2: Filtered FASTQ
├── fastq_porechop_3/           # Step 3: Porechop output
├── multi_tool_qc_4/            # Step 4-5: QC reports
├── host_deconv_out_5/          # Step 6: Decontaminated reads
├── TA1_map_6_V2/               # Step 7: BAM files
├── TA2_depth_7/                # Step 8: Coverage analysis
├── Medaka_consensus_8/         # Step 9: Consensus sequences
├── consensus_ref_viewa_9/      # Step 10: Comparison results
├── variants_call_10/           # Step 11: Variant calls
│   └── variants/unified/       # Final variant results
└── logs/                       # Pipeline logs
```

## Sample Sheet Format

For batch processing, use `samplesheet.csv`:

```csv
sample_id,fastq_path
10090,examples/mini_dataset/raw_fastq/10090_subsampled.fastq
```

Run with:
```bash
./run_pipeline.sh examples/example_config.yaml --batch examples/samplesheet.csv
```

## Verifying Results

Compare your outputs with `expected_outputs/`:

```bash
# Check QC stats match
diff -y examples/mini_dataset/work/multi_tool_qc_4/seqkit/seqkit_stats.tsv \
       examples/expected_outputs/qc_stats.tsv

# Check consensus generated
head examples/mini_dataset/work/Medaka_consensus_8/r2/10090/consensus.fasta
```

## Troubleshooting

### Missing references

If you don't have all reference files, run only the steps that don't require them:

```bash
# Steps 1-4 don't need external references
./run_pipeline.sh examples/example_config.yaml --sample 10090 --from 1 --to 4
```

### Conda environment issues

```bash
# Check available environments
conda info --envs

# Create missing environment
mamba env create -f environment/env_qc.yml
```

### Log files

Check logs for detailed error messages:
```bash
cat examples/mini_dataset/work/logs/pipeline_10090_*.log
```
