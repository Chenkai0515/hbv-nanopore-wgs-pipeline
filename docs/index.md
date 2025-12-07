# HBV Nanopore WGS Pipeline Documentation

Welcome to the documentation for the HBV Nanopore WGS Pipeline.

## Quick Navigation

| Document | Description |
|----------|-------------|
| [Overview](overview.md) | Pipeline architecture and workflow summary |
| [QC Module](qc.md) | Quality control and filtering |
| [Host Decontamination](host_deconv.md) | Host read removal methodology |
| [Mapping & Depth](mapping_depth.md) | Alignment and coverage analysis |
| [Consensus](consensus.md) | Medaka consensus generation |
| [Variants](variants.md) | iVar + Clair3 variant calling |
| [Coordinate Systems](coordinate_systems.md) | Rotated reference and transformations |

## Getting Started

1. **Installation**: See the main [README](../README.md) for setup instructions
2. **Configuration**: Copy and edit `config/config.example.yaml`
3. **Quick Start**: Run `./run_pipeline.sh config/config.yaml --sample SAMPLE_ID`

## Pipeline at a Glance

```
Raw FASTQ
    │
    ├─► 01_qc/         ─► Clean, filtered reads
    │
    ├─► 02_host_deconv/ ─► Viral-enriched reads
    │
    ├─► 03_mapping_consensus/
    │       ├─► Alignment to unified reference
    │       ├─► Coverage analysis
    │       └─► Medaka R1+R2 consensus
    │
    └─► 04_variants/
            ├─► iVar variant calling
            ├─► Clair3 variant calling
            ├─► Joint filtering & tiering
            └─► Unified coordinate output
```

## Key Concepts

### Rotated Reference

This pipeline uses an HBV reference that starts at **DR1** (position 1824 in standard numbering):
- Offset: 1823 bp
- Genome length: 3221 bp

See [Coordinate Systems](coordinate_systems.md) for conversion formulas.

### Dual Variant Calling

We use both **iVar** (amplicon-aware) and **Clair3** (deep learning) for maximum sensitivity:
- DUAL calls: Both tools agree (highest confidence)
- Single-tool calls: Subject to tier-specific filtering

See [Variants](variants.md) for details.

### Quality Tiers

| Tier | AF Range | Confidence |
|------|----------|------------|
| HC | ≥20% | High - suitable for clinical reporting |
| MC | 5-20% | Medium - standard analysis |
| LF | 1-5% | Low - candidates for validation |

## Module Dependencies

```
env_qc ──────────────► 01_qc/
env_base ────────────► 02_host_deconv/
env_base + env_medaka ► 03_mapping_consensus/
env_variants + env_clair3 ► 04_variants/
```

## Support

For issues and questions, please open a GitHub Issue.

