# Expected Outputs

This directory contains example output files to help verify your pipeline runs correctly.

## Files

| File | Step | Description |
|------|------|-------------|
| `qc_stats.tsv` | 4 | SeqKit statistics for the mini dataset |
| `consensus.fasta` | 9 | Truncated consensus sequence (~500bp shown) |
| `variants.csv` | 11 | Example variant table with annotations |

## How to Verify

After running the pipeline, compare your outputs:

```bash
# Check QC stats
diff examples/expected_outputs/qc_stats.tsv \
     examples/mini_dataset/work/multi_tool_qc_4/seqkit/seqkit_stats.tsv

# Check consensus exists and has correct length
seqkit stats examples/mini_dataset/work/Medaka_consensus_8/r2/10090/consensus.fasta

# Check variant output format
head examples/mini_dataset/work/variants_call_10/variants/unified/*.tsv
```

## Notes

- Exact values may vary slightly between software versions
- The consensus sequence should be ~3221 bp (HBV genome length)
- Variant calls depend on read quality and coverage
- Some steps (5, 6) require external databases and may be skipped

## Output Format Details

### QC Stats (SeqKit)

Key columns:
- `num_seqs`: Number of reads after filtering
- `avg_len`: Average read length (expect ~3100-3200 for HBV)
- `AvgQual`: Mean quality score (expect >20 for good data)
- `GC(%)`: GC content (HBV is typically ~48-49%)

### Consensus FASTA

- Header includes reference ID and sample ID
- Sequence length should match HBV genome (~3221 bp)
- Positions with low coverage will be 'N'

### Variant Table

Key columns:
- `pos`: Position in unified reference coordinates
- `ref`/`alt`: Reference and alternate alleles
- `AF_primary`: Allele frequency (0-1)
- `tier`: Confidence tier (HC/MC/LF for SNPs, INDEL-* for indels)
- `decision`: KEEP or CANDIDATE based on filtering rules
