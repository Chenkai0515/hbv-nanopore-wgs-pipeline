# Mini Dataset

This directory is designed to hold subsampled FASTQ files for testing the pipeline.

## Expected Files

```
mini_dataset/
├── README.md
├── 10090_subsampled.fastq
├── 10893_subsampled.fastq
├── 21024_subsampled.fastq
├── 23146_subsampled.fastq
├── 33117_subsampled.fastq
├── 40154_subsampled.fastq
├── 40625_subsampled.fastq
└── 40750_subsampled.fastq
```

## Obtaining Test Data

### Option 1: Download Pre-Subsampled Data

If available, download from:
```bash
# Example (replace with actual URL)
wget https://example.com/hbv_test_data.tar.gz
tar -xzf hbv_test_data.tar.gz -C mini_dataset/
```

### Option 2: Create Your Own Subsample

Use seqkit to subsample from full FASTQ files:

```bash
# Activate environment with seqkit
conda activate hbv_base

# Subsample to ~1000 reads per sample
for fq in /path/to/full/fastq/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    seqkit sample -n 1000 "$fq" > "mini_dataset/${sample}_subsampled.fastq"
done
```

### Option 3: Generate Synthetic Data

For CI/CD testing, generate minimal synthetic data:

```python
# generate_test_data.py
import random

def generate_fastq(filename, n_reads=100, read_length=3000):
    """Generate synthetic FASTQ for testing."""
    bases = "ATCG"
    with open(filename, 'w') as f:
        for i in range(n_reads):
            seq = ''.join(random.choice(bases) for _ in range(read_length))
            qual = 'I' * read_length  # Q40 quality
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")

# Generate for each sample
samples = ['10090', '10893', '21024', '23146', '33117', '40154', '40625', '40750']
for sample in samples:
    generate_fastq(f"mini_dataset/{sample}_subsampled.fastq")
```

## Expected File Sizes

For the subsampled test dataset:

| Sample | Reads | Approx Size |
|--------|-------|-------------|
| 10090 | ~1000 | ~10 MB |
| 10893 | ~1000 | ~10 MB |
| 21024 | ~1000 | ~10 MB |
| 23146 | ~1000 | ~10 MB |
| 33117 | ~1000 | ~10 MB |
| 40154 | ~1000 | ~10 MB |
| 40625 | ~1000 | ~10 MB |
| 40750 | ~1000 | ~10 MB |

**Total: ~80 MB** (uncompressed)

## Data Not Included

The actual sequencing data is not included in this repository because:

1. **Size**: Full FASTQ files can be several GB each
2. **Privacy**: Clinical data may have restrictions
3. **Reproducibility**: Users should validate on their own data

## Verification

After obtaining test data, verify with:

```bash
# Check file existence
ls -la mini_dataset/*.fastq

# Check read counts
seqkit stats mini_dataset/*.fastq

# Expected output format:
# file                              format  type  num_seqs  sum_len  min_len  avg_len  max_len
# mini_dataset/10090_subsampled.fastq  FASTQ   DNA     1000  3000000     2500     3000     3500
```

## Running Tests

Once data is in place:

```bash
# Run single sample test
./run_pipeline.sh config/config.yaml --sample 10090 --test

# Run quick validation
./run_pipeline.sh config/config.yaml --batch examples/samplesheet.csv --dry-run
```

