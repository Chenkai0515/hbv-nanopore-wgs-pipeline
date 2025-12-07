# Quality Control Module

## Overview

The QC module processes raw basecalled FASTQ files through adapter trimming, length/quality filtering, and comprehensive quality assessment.

## Workflow

```
Raw FASTQ → Dorado Trim → Length Filter → Q Filter → Porechop → Clean FASTQ
                                                          │
                                                          ▼
                                         FastQC / NanoPlot / nanoQC / Kraken2
```

## Scripts

### 1. Dorado Adapter Trimming (`1_trim_adapters_dorado_fastq.sh`)

Removes sequencing adapters using ONT's Dorado trimmer.

**Key Features:**
- Kit-specific adapter detection
- Parallel processing support
- Test mode for validation

**Usage:**
```bash
conda activate dorado

bash 1_trim_adapters_dorado_fastq.sh --run \
    --input /path/to/raw_fastq \
    --output /path/to/trimmed \
    --kit SQK-NBD114-24 \
    --jobs 4
```

### 2. Length and Quality Filtering (`2_filter_fastq_gz_2-4k_Q17.py`)

Filters reads by length and mean quality score.

**Default Thresholds:**
- Length: 2000-4000 bp (captures full HBV genome ~3.2kb)
- Quality: Q17+ (≤2% error rate)

**Usage:**
```bash
python 2_filter_fastq_gz_2-4k_Q17.py \
    --input trimmed/ \
    --output filtered/ \
    --min-len 2000 \
    --max-len 4000 \
    --min-quality 17
```

### 3. Porechop Mid-Adapter Removal (`3_run_porechop_midcheck.sh`)

Detects and removes adapters found in the middle of reads (chimeric reads).

**Why This Matters:**
Nanopore sequencing occasionally produces reads with mid-read adapters from:
- Ligation of multiple fragments
- Chimeric library molecules
- Incomplete adapter removal

**Usage:**
```bash
bash 3_run_porechop_midcheck.sh \
    --input filtered/ \
    --output porechop/ \
    --threads 8
```

### 4. Multi-Tool QC (`4_ont_multi_qc.sh`)

Runs comprehensive quality assessment using multiple tools.

**Tools Used:**
| Tool | Output | Purpose |
|------|--------|---------|
| FastQC | HTML/ZIP | Per-base quality, GC content |
| NanoPlot | HTML/PNG | ONT-specific metrics, N50 |
| nanoQC | HTML | Per-position quality heatmap |
| SeqKit | TSV | Summary statistics |
| MultiQC | HTML | Aggregated report |

**Usage:**
```bash
export IN_DIR="/path/to/porechop"
export OUT_DIR="/path/to/qc_results"

bash 4_ont_multi_qc.sh
```

### 5. Kraken2 Contamination Screening (`5_kraken2_nanopore_qc.sh`)

Screens for contamination and classifies reads taxonomically.

**Expected Results for HBV:**
- Majority: Unclassified (HBV not in standard database)
- Some: Human (host carryover)
- Flag: Bacterial taxa (potential contamination)

**Usage:**
```bash
bash 5_kraken2_nanopore_qc.sh \
    --input porechop/ \
    --output kraken_results/ \
    --db /path/to/kraken2_db
```

## Quality Metrics

### Read Length Distribution

For HBV sequencing, expect:
- Peak around 3000-3200 bp (full-length genome)
- Shoulder at ~1600 bp (potential half-genomes)
- Few reads >4000 bp (concatenated or chimeric)

### Quality Score Interpretation

| Q Score | Error Rate | Interpretation |
|---------|------------|----------------|
| Q10 | 10% | Minimum acceptable |
| Q15 | 3% | Good quality |
| Q17 | 2% | **Recommended threshold** |
| Q20 | 1% | High quality |
| Q30 | 0.1% | Excellent |

### Contamination Thresholds

| Category | Threshold | Action |
|----------|-----------|--------|
| Human reads | <50% | Normal |
| Human reads | >50% | Check extraction |
| Bacterial taxa | <1% | Acceptable |
| Bacterial taxa | >5% | Investigate |
| Reagent taxa | <0.1% | Normal background |

## Output Files

```
multi_tool_qc_4/
├── fastqc/
│   ├── {sample}_fastqc.html
│   └── {sample}_fastqc.zip
├── nanoplot/
│   └── {sample}/
│       ├── NanoPlot-report.html
│       ├── LengthvsQualityScatterPlot_dot.png
│       └── NanoStats.txt
├── nanoqc/
│   └── {sample}_nanoQC.html
├── seqkit/
│   └── stats.tsv
├── kraken_raw/
│   └── {sample}.kraken
├── reports/
│   └── {sample}.kreport
├── krona/
│   └── {sample}.html
├── multiqc/
│   └── multiqc_report.html
└── versions.txt
```

## Troubleshooting

### Low read count after filtering

**Possible causes:**
1. Sequencing quality issues
2. Over-stringent filtering
3. Sample degradation

**Solutions:**
- Check raw NanoPlot report first
- Consider relaxing Q threshold to Q15
- Verify library prep quality

### High adapter content

**Possible causes:**
1. Incomplete adapter trimming
2. Wrong kit specified
3. Library prep issues

**Solutions:**
- Verify `--kit` parameter matches actual kit
- Run Porechop with `--discard_middle`
- Review adapter sequences

### Unexpected taxonomy in Kraken

**Interpretation:**
- Low-level reagent contaminants are common
- Compare to `config/kraken_reagent_taxa.tsv`
- Investigate if any taxa >1%

## Best Practices

1. **Always run test mode first** before full processing
2. **Check QC reports** before proceeding to host decontamination
3. **Document any parameter changes** from defaults
4. **Keep raw data** until analysis is complete and validated

