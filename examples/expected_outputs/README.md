# Expected Outputs

This directory contains reference outputs for validating pipeline results.

## Purpose

Use these expected outputs to:
1. Verify your installation is working correctly
2. Validate results after code changes
3. Debug issues by comparing to known-good outputs

## Expected Metrics

### QC Statistics

After running the QC module, expect approximately:

| Sample | Input Reads | After Filter | % Retained |
|--------|-------------|--------------|------------|
| 10090 | ~1000 | ~800-900 | 80-90% |
| 10893 | ~1000 | ~800-900 | 80-90% |
| ... | ... | ... | ... |

### Host Decontamination

| Sample | Viral Enriched | Host Only | % Viral |
|--------|----------------|-----------|---------|
| 10090 | ~600-800 | ~100-200 | 60-80% |
| 10893 | ~600-800 | ~100-200 | 60-80% |
| ... | ... | ... | ... |

### Coverage

| Sample | Mean Depth | % ≥10× | % ≥100× |
|--------|------------|--------|---------|
| 10090 | ~500× | >99% | >95% |
| 10893 | ~500× | >99% | >95% |
| ... | ... | ... | ... |

### Variants Detected

| Sample | Total Variants | HC | MC | LF |
|--------|----------------|----|----|-----|
| 10090 | 0-2 | 0 | 0-1 | 0-1 |
| 10893 | 3-8 | 1-3 | 1-3 | 1-2 |
| 21024 | 0-2 | 0 | 0-1 | 0-1 |
| 23146 | 0-2 | 0 | 0-1 | 0-1 |
| 33117 | 0-2 | 0 | 0-1 | 0-1 |
| 40154 | 3-8 | 1-3 | 1-3 | 1-2 |
| 40625 | 3-8 | 1-3 | 1-3 | 1-2 |
| 40750 | 3-8 | 1-3 | 1-3 | 1-2 |

## Reference Files

When full validation data is available, this directory will contain:

```
expected_outputs/
├── qc_summary.tsv           # Expected QC metrics
├── host_deconv_summary.tsv  # Expected host decontamination metrics
├── coverage_summary.tsv     # Expected coverage statistics
└── variants_summary.tsv     # Expected variant counts per sample
```

## Validation Script

Use this script to compare your results to expected outputs:

```bash
#!/bin/bash
# validate_outputs.sh

RESULTS_DIR="$1"
EXPECTED_DIR="examples/expected_outputs"

echo "Comparing outputs..."

# Compare QC summary
if [ -f "$EXPECTED_DIR/qc_summary.tsv" ]; then
    diff <(sort "$RESULTS_DIR/qc/summary.tsv") \
         <(sort "$EXPECTED_DIR/qc_summary.tsv") && \
    echo "QC: PASS" || echo "QC: DIFFERS"
fi

# Compare variant counts
if [ -f "$EXPECTED_DIR/variants_summary.tsv" ]; then
    diff <(sort "$RESULTS_DIR/variants/summary_counts.tsv") \
         <(sort "$EXPECTED_DIR/variants_summary.tsv") && \
    echo "Variants: PASS" || echo "Variants: DIFFERS"
fi

echo "Validation complete."
```

## Notes

- Expected values are approximate due to:
  - Stochastic elements in some algorithms
  - Minor version differences in tools
  - Platform-specific floating-point variations

- Results within ±10% of expected values are generally acceptable

- For exact reproducibility:
  - Use the same tool versions
  - Set random seeds where available
  - Use identical reference files

