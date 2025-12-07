#!/usr/bin/env bash
# Kraken2 contamination screening for Nanopore reads
# Computes human/bacterial/reagent fractions with visualization
set -euo pipefail

# =========== Configuration ===========
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

INDIR="${INDIR:-${PROJECT_DIR}/fastq_porechop_3}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/multi_tool_qc_4}"

KRAKEN_DB="${KRAKEN_DB:-/path/to/kraken2_db}"
THREADS="${THREADS:-16}"
USE_MEMORY_MAPPING="${USE_MEMORY_MAPPING:-0}"
KRAKEN_CONFIDENCE="${KRAKEN_CONFIDENCE:-0.00}"
KRAKEN_ADDITIONAL_OPTS="${KRAKEN_ADDITIONAL_OPTS:-}"

# =========== Directory setup ===========
REPORT_DIR="$OUTDIR/reports"
KRAKEN_DIR="$OUTDIR/kraken_raw"
KRONA_DIR="$OUTDIR/krona"
SUMMARY_DIR="$OUTDIR/summary"

mkdir -p "$REPORT_DIR" "$KRAKEN_DIR" "$KRONA_DIR" "$SUMMARY_DIR"

# =========== Dependency check ===========
command -v kraken2 >/dev/null 2>&1 || { echo "kraken2 not found"; exit 1; }
python3 - <<'PY'
try:
    import pandas, matplotlib
except Exception as e:
    import sys
    sys.stderr.write("Missing Python packages: pandas, matplotlib\n")
    sys.exit(1)
PY

if [[ ! -d "$KRAKEN_DB" ]]; then
  echo "Kraken2 database not found: $KRAKEN_DB"
  echo "Set KRAKEN_DB environment variable or edit this script"
  exit 1
fi

# =========== Check for Krona ===========
HAVE_KRONA=0
if command -v ktImportTaxonomy >/dev/null 2>&1; then
  HAVE_KRONA=1
fi

# =========== Reagent taxa list ===========
REAGENT_TSV="$SUMMARY_DIR/reagent_taxa.tsv"
if [[ ! -s "$REAGENT_TSV" ]]; then
cat > "$REAGENT_TSV" <<'EOF'
#rank	taxid	name	notes
genus	48736	Ralstonia	kit/reagent common
genus	13687	Sphingomonas	kit/reagent common
genus	374	Bradyrhizobium	kit/reagent common
genus	407	Methylobacterium	kit/reagent common
genus	80865	Delftia	kit/reagent common
genus	32008	Burkholderia	kit/reagent common
genus	286	Pseudomonas	lab/kit common
genus	469	Acinetobacter	lab/kit common
genus	1386	Bacillus	lab/kit common
genus	40323	Stenotrophomonas	lab/kit common
genus	165696	Novosphingobium	lab/kit common
genus	106589	Cupriavidus	lab/kit common
genus	1716	Corynebacterium	skin/lab common
genus	1269	Micrococcus	skin/lab common
genus	1743	Propionibacterium	skin/lab common (old)
genus	1912216	Cutibacterium	skin/lab common (new)
species	10710	Enterobacteria phage lambda	ONT/other controls
species	10847	Enterobacteria phage phiX174	Illumina control
EOF
fi

# =========== Collect input files ===========
echo "Scanning: $INDIR"
mapfile -d '' FASTQS < <(find "$INDIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -print0 | sort -z)
if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "No fastq.gz files found in $INDIR"
  exit 1
fi
echo "Found ${#FASTQS[@]} fastq.gz files"

# =========== Kraken2 classification ===========
for fq in "${FASTQS[@]}"; do
  base="$(basename "$fq")"
  sample="${base%.fastq.gz}"
  sample="${sample%.fq.gz}"

  report="$REPORT_DIR/${sample}.kreport"
  kraken_out="$KRAKEN_DIR/${sample}.kraken"

  if [[ -s "$report" && -s "$kraken_out" ]]; then
    echo "[SKIP] Already exists: $sample"
    continue
  fi

  echo "[Kraken2] Processing: $sample"
  mmopt=()
  [[ "$USE_MEMORY_MAPPING" == "1" ]] && mmopt+=(--memory-mapping)

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --use-names \
    --gzip-compressed \
    --report "$report" \
    --output "$kraken_out" \
    --confidence "$KRAKEN_CONFIDENCE" \
    "${mmopt[@]}" \
    $KRAKEN_ADDITIONAL_OPTS \
    --classified-out /dev/null \
    --unclassified-out /dev/null \
    -- \
    "$fq"

  if [[ "$HAVE_KRONA" -eq 1 ]]; then
    echo "[Krona] Generating: $sample.krona.html"
    cut -f2,3 "$kraken_out" | ktImportTaxonomy -o "$KRONA_DIR/${sample}.krona.html" - || true
  fi
done

# =========== Summary statistics (Python) ===========
python3 - "$REPORT_DIR" "$SUMMARY_DIR" "$REAGENT_TSV" <<'PY'
import sys, os, csv
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

report_dir, summary_dir, reagent_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

TAX = {
    "root": 1, "unclassified": 0, "bacteria": 2, "archaea": 2157,
    "eukaryota": 2759, "viruses": 10239, "human": 9606,
}

# Load reagent taxids
reagent = []
with open(reagent_tsv, newline='') as fh:
    for row in fh:
        row = row.strip()
        if not row or row.startswith("#"):
            continue
        parts = row.split("\t")
        if len(parts) >= 3:
            taxid = parts[1].strip()
            if taxid.isdigit():
                reagent.append(int(taxid))
REAGENT_SET = set(reagent)

def parse_kreport(path):
    tax_clade = {}
    unclassified = 0
    classified = 0
    with open(path, newline='') as fh:
        tsv = csv.reader(fh, delimiter='\t')
        for row in tsv:
            if len(row) < 6:
                continue
            _, clade, _, _, taxid, _ = row
            try:
                taxid = int(taxid)
                clade = int(clade)
            except:
                continue
            tax_clade[taxid] = clade
            if taxid == 0:
                unclassified = clade
            elif taxid == 1:
                classified = clade
    total = unclassified + classified
    return tax_clade, unclassified, classified, total

rows = []
report_files = sorted([p for p in os.listdir(report_dir) if p.endswith(".kreport")])
if not report_files:
    sys.stderr.write("No .kreport files found\n")
    sys.exit(1)

for rep in report_files:
    sample = rep.replace(".kreport","")
    path = os.path.join(report_dir, rep)
    tax_clade, uncls, clsf, total = parse_kreport(path)

    get = lambda tid: tax_clade.get(tid, 0)
    bacteria = get(TAX["bacteria"])
    archaea  = get(TAX["archaea"])
    euk     = get(TAX["eukaryota"])
    viruses = get(TAX["viruses"])
    human   = get(TAX["human"])

    reagent_reads = sum(tax_clade.get(tid, 0) for tid in REAGENT_SET)
    bacterial_non_reagent = max(0, bacteria - reagent_reads)
    viral_non_reagent     = max(0, viruses - reagent_reads)
    other_euk             = max(0, euk - human)
    known = human + bacterial_non_reagent + reagent_reads + viral_non_reagent + other_euk + archaea
    other_classified = max(0, clsf - known)

    d = OrderedDict()
    d["sample"] = sample
    d["total_reads"] = total
    d["classified_reads"] = clsf
    d["unclassified_reads"] = uncls
    d["human"] = human
    d["bacterial_non_reagent"] = bacterial_non_reagent
    d["reagent"] = reagent_reads
    d["viral_non_reagent"] = viral_non_reagent
    d["other_eukaryote"] = other_euk
    d["archaea"] = archaea
    d["other_classified"] = other_classified

    for k in list(d.keys()):
        if k in ("sample","total_reads","classified_reads"):
            continue
        d[k + "_frac"] = (d[k] / total) if total > 0 else 0.0
    rows.append(d)

df = pd.DataFrame(rows).sort_values("sample")

if "unclassified_reads_frac" not in df.columns and "unclassified_reads" in df.columns:
    df["unclassified_reads_frac"] = df["unclassified_reads"] / df["total_reads"]

counts_cols = ["human","bacterial_non_reagent","reagent","viral_non_reagent",
               "other_eukaryote","archaea","other_classified","unclassified_reads"]
frac_cols = [c + "_frac" for c in counts_cols]

counts_path = os.path.join(summary_dir, "kraken_category_counts.tsv")
frac_path   = os.path.join(summary_dir, "kraken_category_fractions.tsv")
df[["sample","total_reads","classified_reads"] + counts_cols].to_csv(counts_path, sep="\t", index=False)
df[["sample"] + [c for c in frac_cols if c in df.columns]].to_csv(frac_path, sep="\t", index=False)

# Stacked bar plot
plot_order = ["human_frac","bacterial_non_reagent_frac","reagent_frac",
              "viral_non_reagent_frac","other_eukaryote_frac","archaea_frac",
              "other_classified_frac","unclassified_reads_frac"]
plot_cols = [c for c in plot_order if c in df.columns]
plot_df = df[["sample"] + plot_cols].set_index("sample")

plt.figure(figsize=(max(8, len(plot_df)*0.5), 6))
cum = None
for col in plot_df.columns:
    vals = plot_df[col].values
    if cum is None:
        plt.bar(plot_df.index, vals, label=col.replace("_frac",""))
        cum = vals
    else:
        plt.bar(plot_df.index, vals, bottom=cum, label=col.replace("_frac",""))
        cum = cum + vals

plt.xticks(rotation=75, ha='right')
plt.ylabel("Fraction of reads")
plt.title("Kraken2 category fractions per sample")
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.tight_layout()
png_path = os.path.join(summary_dir, "kraken_category_fractions_stacked.png")
plt.savefig(png_path, dpi=180)
plt.close()

print(f"Done. Summary in: {summary_dir}")
PY

echo "Kraken2 QC complete. Output: $SUMMARY_DIR"
