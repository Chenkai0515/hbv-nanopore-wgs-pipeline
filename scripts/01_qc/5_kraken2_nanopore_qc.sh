#!/usr/bin/env bash
# kraken2_nanopore_qc.sh
# Run Kraken2 on Nanopore fastq.gz, compute human/bacterial/reagent fractions,
# and visualize per-sample and cohort summaries.

set -euo pipefail

### ----------- 用户可修改的参数 -----------
INDIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_porechop_3"
OUTDIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/multi_tool_qc_4"

# Kraken2 数据库路径（必须：请改为你的本地 Kraken2 DB 路径，需含 human）
KRAKEN_DB="${KRAKEN_DB:-/Users/jck/Desktop/workflow-qc/raw-high-accuracy-fastq/sc/standard_8gb}"

# 计算资源
THREADS="${THREADS:-16}"        # 使用 16 核
USE_MEMORY_MAPPING="${USE_MEMORY_MAPPING:-0}"  # 0=关闭(更快，需足够内存)；1=开启(更省内存)

# Kraken2 其它常用参数（可按需调整）
KRAKEN_CONFIDENCE="${KRAKEN_CONFIDENCE:-0.00}"  # Nanopore 建议保持默认或略低
KRAKEN_ADDITIONAL_OPTS="${KRAKEN_ADDITIONAL_OPTS:-}"

### ----------- 目录布局 -----------
REPORT_DIR="$OUTDIR/reports"
KRAKEN_DIR="$OUTDIR/kraken_raw"
KRONA_DIR="$OUTDIR/krona"
SUMMARY_DIR="$OUTDIR/summary"

mkdir -p "$REPORT_DIR" "$KRAKEN_DIR" "$KRONA_DIR" "$SUMMARY_DIR"

### ----------- 依赖检查（只检查名称是否可用）-----------
command -v kraken2 >/dev/null 2>&1 || { echo "未找到 kraken2，请先安装。"; exit 1; }
python3 - <<'PY'
try:
    import pandas, matplotlib
except Exception as e:
    import sys
    sys.stderr.write("缺少 Python 包：需要 pandas、matplotlib\n")
    sys.exit(1)
PY

if [[ ! -d "$KRAKEN_DB" ]]; then
  echo "找不到 Kraken2 数据库目录：$KRAKEN_DB"
  echo "请设置正确的 KRAKEN_DB 环境变量或在脚本中修改 KRAKEN_DB。"
  exit 1
fi

### ----------- 试探 Krona 工具（可选）-----------
HAVE_KRONA=0
if command -v ktImportTaxonomy >/dev/null 2>&1; then
  HAVE_KRONA=1
fi

### ----------- 内置 reagent（试剂污染/kitome）候选清单 -----------
# 采用常见试剂/水系/实验室常见污染属/对照的 NCBI TaxID（以属为主以避免重复计数）
# 你可以根据自己的项目编辑/扩充此清单（位于 $SUMMARY_DIR/reagent_taxa.tsv）
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
# 常见对照/载体噬菌体等
species	10710	Enterobacteria phage lambda	ONT/other controls
species	10847	Enterobacteria phage phiX174	Illumina control
EOF
fi

### ----------- 汇总输入列表 -----------
echo "扫描输入目录：$INDIR"
mapfile -d '' FASTQS < <(find "$INDIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -print0 | sort -z)
if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "在 $INDIR 未找到 fastq.gz 文件。"
  exit 1
fi
echo "共发现 ${#FASTQS[@]} 个 fastq.gz 文件。"

### ----------- Kraken2 分类循环 -----------
for fq in "${FASTQS[@]}"; do
  base="$(basename "$fq")"
  sample="$base"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fq.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"

  report="$REPORT_DIR/${sample}.kreport"
  kraken_out="$KRAKEN_DIR/${sample}.kraken"

  if [[ -s "$report" && -s "$kraken_out" ]]; then
    echo "[跳过] 已存在：$sample"
    continue
  fi

  echo "[Kraken2] 处理样本：$sample"
  mmopt=()
  if [[ "$USE_MEMORY_MAPPING" == "1" ]]; then
    mmopt+=(--memory-mapping)
  fi

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
    echo "[Krona] 生成：$sample.krona.html"
    # Krona 可直接读取 Kraken 输出
    cut -f2,3 "$kraken_out" | ktImportTaxonomy -o "$KRONA_DIR/${sample}.krona.html" - || true

  fi
done

### ----------- 统计/可视化（Python）-----------
python3 - "$REPORT_DIR" "$SUMMARY_DIR" "$REAGENT_TSV" <<'PY'
import sys, os, math, csv
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

report_dir, summary_dir, reagent_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

TAX = {
    "root": 1,
    "unclassified": 0,
    "bacteria": 2,
    "archaea": 2157,
    "eukaryota": 2759,
    "viruses": 10239,
    "human": 9606,
}

# 读取 reagent taxids
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
    # Kraken2 report: pct, clade_reads, direct_reads, rank, taxid, name
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
    sys.stderr.write("未发现 .kreport 文件。\n")
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

    d["cat3_human"] = human
    d["cat3_bacterial"] = bacterial_non_reagent
    d["cat3_reagent"] = reagent_reads

    # 计算占比（不跳过 unclassified_reads）
    for k in list(d.keys()):
        if k in ("sample","total_reads","classified_reads"):
            continue
        d[k + "_frac"] = (d[k] / total) if total > 0 else 0.0

    rows.append(d)

df = pd.DataFrame(rows).sort_values("sample")

# 兜底：确保存在 unclassified_reads_frac
if "unclassified_reads_frac" not in df.columns and "unclassified_reads" in df.columns and "total_reads" in df.columns:
    df["unclassified_reads_frac"] = df["unclassified_reads"] / df["total_reads"]

counts_cols = [
    "human","bacterial_non_reagent","reagent","viral_non_reagent",
    "other_eukaryote","archaea","other_classified","unclassified_reads"
]
frac_cols = [c + "_frac" for c in counts_cols]

counts_path = os.path.join(summary_dir, "kraken_category_counts.tsv")
frac_path   = os.path.join(summary_dir, "kraken_category_fractions.tsv")
df_out_counts = df[["sample","total_reads","classified_reads"] + counts_cols]
df_out_fracs  = df[["sample"] + [c for c in frac_cols if c in df.columns]]
df_out_counts.to_csv(counts_path, sep="\t", index=False)
df_out_fracs.to_csv(frac_path, sep="\t", index=False)

# 三大类
cat3_cols = ["cat3_human","cat3_bacterial","cat3_reagent",
             "cat3_human_frac","cat3_bacterial_frac","cat3_reagent_frac"]
have_cat3 = [c for c in cat3_cols if c in df.columns]
df_cat3 = df[["sample"] + have_cat3]
df_cat3.to_csv(os.path.join(summary_dir, "fractions_3categories.tsv"), sep="\t", index=False)

# 画图（堆叠占比）
plot_order = ["human_frac","bacterial_non_reagent_frac","reagent_frac",
              "viral_non_reagent_frac","other_eukaryote_frac","archaea_frac","other_classified_frac","unclassified_reads_frac"]
plot_cols = [c for c in plot_order if c in df_out_fracs.columns]
plot_df = df[["sample"] + plot_cols].set_index("sample")

plt.figure(figsize=(max(8, len(plot_df)*0.5), 6))
cum = None
for i, col in enumerate(plot_df.columns):
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
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", borderaxespad=0.)
plt.tight_layout()
png_path = os.path.join(summary_dir, "kraken_category_fractions_stacked.png")
plt.savefig(png_path, dpi=180)
plt.close()

with open(os.path.join(summary_dir, "README_summary.txt"), "w") as fh:
    fh.write("Files written:\n")
    fh.write(f"- {counts_path}\n")
    fh.write(f"- {frac_path}\n")
    fh.write(f"- {os.path.join(summary_dir, 'fractions_3categories.tsv')}\n")
    fh.write(f"- {png_path}\n")
    fh.write("- Krona HTML per sample (if KronaTools installed) in ../krona/\n")

print("完成统计与可视化。汇总见：", summary_dir)

