#!/usr/bin/env bash
# Multi-tool QC for ONT reads: FastQC, NanoPlot, nanoQC, SeqKit
set -Eeuo pipefail

# ====== Path Configuration ======
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

IN_DIR="${IN_DIR:-${PROJECT_DIR}/fastq_porechop_3}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/multi_tool_qc_4}"

# ====== Resource Configuration ======
TOTAL_CORES="${TOTAL_CORES:-16}"
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"
FASTQC_THREADS="${FASTQC_THREADS:-$TOTAL_CORES}"

NANOPLOT_THREADS_PER_JOB="${NANOPLOT_THREADS_PER_JOB:-2}"
NANOPLOT_JOBS=$(( TOTAL_CORES / NANOPLOT_THREADS_PER_JOB ))
if (( NANOPLOT_JOBS < 1 )); then NANOPLOT_JOBS=1; fi

NANOQC_JOBS="${NANOQC_JOBS:-$TOTAL_CORES}"
SEQKIT_THREADS="${SEQKIT_THREADS:-$TOTAL_CORES}"
NANOQC_MINLEN="${NANOQC_MINLEN:-0}"

# ====== Tool binaries ======
FASTQC_BIN="${FASTQC_BIN:-fastqc}"
NANOPLOT_BIN="${NANOPLOT_BIN:-NanoPlot}"
NANOQC_BIN="${NANOQC_BIN:-nanoQC}"
SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"

mkdir -p "$OUT_DIR"/{logs,fastqc,nanoplot,nanoqc,seqkit}

# ====== Dependency check ======
for exe in "$FASTQC_BIN" "$NANOPLOT_BIN" "$NANOQC_BIN" "$SEQKIT_BIN"; do
  if ! command -v "$exe" >/dev/null 2>&1; then
    echo "[ERROR] Executable not found: $exe" >&2
    exit 1
  fi
done

# ====== Collect input FASTQ files ======
FASTQS=()
while IFS= read -r -d '' f; do FASTQS+=("$f"); done < <(
  find "$IN_DIR" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -print0
)
if (( ${#FASTQS[@]} == 0 )); then
  echo "[ERROR] No FASTQ files found in $IN_DIR" >&2
  exit 1
fi

# ====== Record versions ======
{
  echo "Started: $(date)"
  "$FASTQC_BIN" --version || true
  "$NANOPLOT_BIN" --version || true
  "$NANOQC_BIN" --version || true
  "$SEQKIT_BIN" version || true
} | tee "$OUT_DIR/versions.txt"

# ====== Step 1: FastQC ======
echo "[INFO] Running FastQC with $FASTQC_THREADS threads ..."
env _JAVA_OPTIONS="-Xmx${JAVA_MEM_GB}g -Xms512m" \
  "$FASTQC_BIN" -t "$FASTQC_THREADS" -o "$OUT_DIR/fastqc" "${FASTQS[@]}" \
  2>&1 | tee "$OUT_DIR/logs/fastqc.log"

# ====== Step 2: NanoPlot ======
echo "[INFO] Running NanoPlot: $NANOPLOT_JOBS jobs, $NANOPLOT_THREADS_PER_JOB threads/job ..."
printf '%s\0' "${FASTQS[@]}" | xargs -0 -n 1 -P "$NANOPLOT_JOBS" -I {} bash -c '
  set -Eeuo pipefail
  in="$1"
  bn=$(basename "$in")
  sample=${bn%.fastq.gz}; sample=${sample%.fq.gz}
  od="'"$OUT_DIR"'/nanoplot/${sample}"
  log="'"$OUT_DIR"'/logs/nanoplot_${sample}.log"
  mkdir -p "$od"
  "'"$NANOPLOT_BIN"'" -t "'"$NANOPLOT_THREADS_PER_JOB"'" --fastq "$in" \
    --tsv_stats --N50 --loglength --huge \
    -o "$od" \
    > "$log" 2>&1
' _ {}

# ====== Step 3: nanoQC ======
echo "[INFO] Running nanoQC: $NANOQC_JOBS jobs ..."
printf '%s\0' "${FASTQS[@]}" | xargs -0 -n 1 -P "$NANOQC_JOBS" -I {} bash -c '
  set -Eeuo pipefail
  in="$1"
  bn=$(basename "$in")
  sample=${bn%.fastq.gz}; sample=${sample%.fq.gz}
  od="'"$OUT_DIR"'/nanoqc/${sample}"
  log="'"$OUT_DIR"'/logs/nanoqc_${sample}.log"
  mkdir -p "$od"
  if [ "'"$NANOQC_MINLEN"'" -gt 0 ]; then
    "'"$NANOQC_BIN"'" -o "$od" -l "'"$NANOQC_MINLEN"'" "$in" > "$log" 2>&1
  else
    "'"$NANOQC_BIN"'" -o "$od" "$in" > "$log" 2>&1
  fi
' _ {}

# ====== Step 4: SeqKit stats ======
echo "[INFO] Running SeqKit stats with $SEQKIT_THREADS threads ..."
"$SEQKIT_BIN" stats -a -T -j "$SEQKIT_THREADS" "${FASTQS[@]}" \
  | tee "$OUT_DIR/seqkit/seqkit_stats.tsv" > /dev/null

# ====== Step 5: MultiQC (optional) ======
if command -v multiqc >/dev/null 2>&1; then
  echo "[INFO] Running MultiQC to aggregate reports ..."
  multiqc "$OUT_DIR" -o "$OUT_DIR/multiqc" 2>&1 | tee "$OUT_DIR/logs/multiqc.log"
else
  echo "[INFO] MultiQC not found, skipping (install: conda install -c bioconda multiqc)"
fi

echo "[DONE] All QC complete. Output: $OUT_DIR"
