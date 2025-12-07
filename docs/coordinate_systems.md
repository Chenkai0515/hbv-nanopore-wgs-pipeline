# Coordinate Systems

## Overview

This pipeline uses a **rotated HBV reference** to handle the circular genome. Understanding the coordinate system is essential for interpreting results and comparing with external databases.

## The HBV Circular Genome Problem

HBV has a circular DNA genome (~3.2 kb). For computational analysis, we must linearize it:

```
        Circular HBV Genome
        ───────────────────

              ┌─────┐
             ╱       ╲
            │  HBV   │
            │ 3221bp │
             ╲       ╱
              └─────┘

        Where to "cut" for linearization?
```

### Standard Numbering (EcoRI site)

Traditional HBV coordinates start at the EcoRI restriction site:
- Position 1 = EcoRI site
- Core gene starts around position 1814
- PreS1 starts around position 2848

### This Pipeline (DR1 start)

We start at **Direct Repeat 1 (DR1)** for several reasons:

1. **Natural breakpoint**: DR1 is the replication origin
2. **No gene splitting**: Core gene is not split across the boundary
3. **Integration analysis**: DR1 is a common integration site

## Coordinate Transformation

### Key Parameters

```yaml
OFFSET: 1823        # Bases to add for conversion
GENOME_LENGTH: 3221  # Total genome size
```

### Conversion Formulas

**Rotated → Original (Standard):**
```python
original_pos = ((rotated_pos - 1) + OFFSET) % GENOME_LENGTH + 1
```

**Original → Rotated:**
```python
rotated_pos = ((original_pos - 1) - OFFSET) % GENOME_LENGTH + 1
# If result ≤ 0, add GENOME_LENGTH
```

### Example Conversions

| Rotated Position | Calculation | Original Position |
|------------------|-------------|-------------------|
| 1 | (0 + 1823) % 3221 + 1 | 1824 |
| 100 | (99 + 1823) % 3221 + 1 | 1923 |
| 1500 | (1499 + 1823) % 3221 + 1 | 101 |
| 3000 | (2999 + 1823) % 3221 + 1 | 1601 |

## Visual Representation

```
Original HBV Reference (Standard Numbering):
Position:  1────500────1000────1500────1824────2500────3000────3221
           │                      │      ▲                        │
           │      PreS/S          │      │        Core            │
           │                      │   DR1 (our start)             │
           └──────────────────────┴──────┴────────────────────────┘

Rotated HBV Reference (This Pipeline):
Position:  1────500────1000────1398────1500────2000────2500────3221
           │                      │                               │
           │      Core            │      PreS/S                   │
           │  (no longer split)   │                               │
           └──────────────────────┴───────────────────────────────┘
           ▲
           └── Starts at DR1 (position 1824 in original)
```

## Gene Positions

### In Rotated Reference (This Pipeline)

| Gene | Start | End | Notes |
|------|-------|-----|-------|
| Core | 1 | 638 | Complete, not split |
| X | 639 | 1097 | Regulatory protein |
| PreS1 | 1025 | 1389 | Overlaps X |
| PreS2 | 1390 | 1554 | Surface antigen |
| S | 1555 | 2241 | Major surface antigen |
| Pol | 1 | 2592 | Spans most of genome |

### Conversion Table for Key Sites

| Feature | Original Position | Rotated Position |
|---------|-------------------|------------------|
| DR1 | 1824 | 1 |
| DR2 | 1590 | 2988 |
| EcoRI site | 1 | 1399 |
| Core ATG | 1814 | 3212 |
| PreCore ATG | 1785 | 3183 |
| BCP (1762) | 1762 | 3160 |
| BCP (1764) | 1764 | 3162 |
| PreCore (1896) | 1896 | 73 |

## View A + View B Transformation

### The Two Views

```
┌─────────────────────────────────────────────────────────────────┐
│                     COORDINATE VIEWS                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   VIEW A: Consensus vs Unified Reference                        │
│   ──────────────────────────────────────                        │
│   Coordinates: Rotated reference positions                      │
│   Content: Fixed differences (genotype, consensus-level)        │
│                                                                 │
│   VIEW B: Reads vs Per-Sample Consensus                         │
│   ──────────────────────────────────────                        │
│   Coordinates: Consensus positions (may differ from ref)        │
│   Content: Within-sample variants (quasispecies)                │
│                                                                 │
│   UNIFIED: Final Combined Output                                │
│   ──────────────────────────────                                │
│   Coordinates: Rotated reference positions                      │
│   Content: All variants in consistent coordinate system         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Transformation Pipeline

```
Step 11.6: Consensus coords → Rotated ref coords
           ─────────────────────────────────────
           Uses alignment of consensus to reference
           Handles insertions/deletions in consensus

Step 11.7: Merge View A + Transformed View B
           ─────────────────────────────────────
           Combines all variant sources
           Resolves coordinate overlaps

Step 11.8: Flip to unified output format
           ─────────────────────────────────────
           Final coordinate validation
           Standard output columns
```

## Handling Insertions/Deletions

### The Coordinate Shift Problem

When consensus has indels relative to reference:

```
Reference:  ATG---GATCGA  (positions 1-9)
Consensus:  ATGAAAGAT-GA  (positions 1-11)
                 ^^^  ^
                 ins  del

Variant at consensus position 8:
- Consensus coords: position 8
- Reference coords: position 5 (after accounting for insertion)
```

### Solution: Alignment-Based Mapping

The transformation scripts use pairwise alignment to map positions:

```python
def transform_position(cons_pos, alignment):
    """Map consensus position to reference position."""
    ref_pos = 0
    cons_pos_current = 0
    
    for ref_char, cons_char in alignment:
        if cons_char != '-':
            cons_pos_current += 1
        if ref_char != '-':
            ref_pos += 1
        if cons_pos_current == cons_pos:
            return ref_pos
    
    return None  # Position beyond alignment
```

## Output Coordinate Columns

### Unified Variants CSV

| Column | Description | Example |
|--------|-------------|---------|
| `unified_pos` | Position in rotated reference | 1762 |
| `original_pos` | Position in standard numbering | 3160 |
| `consensus_pos` | Original position in consensus | 1758 |
| `mapping_note` | Transformation notes | "ins_shift:+4" |

## Practical Examples

### Example 1: BCP Mutation

```
Detecting A1762T mutation:

Standard position: 1762
Rotated position: (1762 - 1) - 1823 + 3221 + 1 = 3160

In pipeline output:
- unified_pos: 3160
- original_pos: 1762
- ref: A
- alt: T
```

### Example 2: PreCore Mutation

```
Detecting G1896A mutation:

Standard position: 1896
Rotated position: (1896 - 1) - 1823 + 1 = 73

In pipeline output:
- unified_pos: 73
- original_pos: 1896
- ref: G
- alt: A
```

### Example 3: Variant with Consensus Indel

```
Sample has 3bp insertion at rotated position 500

Read variant at consensus position 600:
- Consensus position: 600
- After transformation: 600 - 3 = 597 (rotated ref)
- Original position: (597 - 1 + 1823) % 3221 + 1 = 2420

Output includes mapping_note: "ins_adjustment:-3"
```

## Comparison with Databases

When comparing pipeline results with external databases:

1. **Check the reference**: What coordinate system does the database use?
2. **Convert if needed**: Use the formulas above
3. **Verify gene positions**: Cross-check known mutations

### Common Databases

| Database | Coordinate System | Conversion Needed |
|----------|-------------------|-------------------|
| HBVdb | Standard (EcoRI) | Yes |
| NCBI RefSeq | Standard (EcoRI) | Yes |
| Geno2Pheno | Standard (EcoRI) | Yes |

## Troubleshooting

### Position doesn't match expected

1. Check which coordinate system you're using
2. Verify the offset (1823 for this pipeline)
3. Account for any indels in consensus

### Variant not found at expected position

1. May be in a different reading frame due to indels
2. Check the `mapping_note` column
3. Review the consensus alignment

### Cross-sample positions don't align

1. Different samples may have different indel profiles
2. Always compare using unified_pos, not consensus_pos
3. Check View A differences for each sample

