#!/usr/bin/env bash
# Round 2 Medaka consensus polishing
# Uses R1 consensus as reference for further refinement
set -euo pipefail

# ========= Configuration =========
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v5.0.0}"
R1_CONS_DIR="${R1_CONS_DIR:-${PROJECT_DIR}/Medaka_consensus_8/r1/medaka_R1_consensus}"
READS_ROOT="${READS_ROOT:-${PROJECT_DIR}/host_deconv_out_5}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_DIR}/Medaka_consensus_8/r2}"

THREADS_PER_JOB="${THREADS_PER_JOB:-1}"

if command -v nproc >/dev/null 2>&1; then
  CORES=$(nproc)
else
  CORES=$(sysctl -n hw.ncpu)
fi
JOBS="${JOBS:-$CORES}"

# ========= Dependency check =========
command -v medaka_consensus >/dev/null 2>&1 || {
  echo "[ERROR] medaka_consensus not found" >&2
  exit 1
}
mkdir -p "$OUT_ROOT"

# ========= Generate commands =========
CMDS_FILE="$OUT_ROOT/medaka_r2_commands.txt"
: > "$CMDS_FILE"

while IFS= read -r -d '' SAMPLE_DIR; do
  SAMPLE_ID="$(basename "$SAMPLE_DIR")"
  REF="$SAMPLE_DIR/consensus.fasta"
  
  READS_SUBDIR="${READS_ROOT}/${SAMPLE_ID}_subsampled.trimmed_filtered.porechop"
  READS="${READS_SUBDIR}/${SAMPLE_ID}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz"
  
  OUTDIR="$OUT_ROOT/$SAMPLE_ID"

  if [[ ! -s "$REF" ]]; then
    echo "echo '[SKIP] $SAMPLE_ID: missing R1 consensus $REF' 1>&2" >> "$CMDS_FILE"
    continue
  fi
  if [[ ! -s "$READS" ]]; then
    echo "echo '[SKIP] $SAMPLE_ID: missing reads $READS' 1>&2" >> "$CMDS_FILE"
    continue
  fi

  mkdir -p "$OUTDIR"
  printf "medaka_consensus -i '%s' -d '%s' -o '%s' -m '%s' -t %d && \
find '%s' -maxdepth 1 -type f \\( -name '*.bam' -o -name '*.bai' \\) -delete\n" \
    "$READS" "$REF" "$OUTDIR" "$MODEL" "$THREADS_PER_JOB" "$OUTDIR" >> "$CMDS_FILE"
done < <(find "$R1_CONS_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

# ========= Parallel execution =========
echo "[INFO] Running $JOBS parallel jobs, $THREADS_PER_JOB threads each"
if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" --bar < "$CMDS_FILE"
else
  cat "$CMDS_FILE" | xargs -I CMD -P "$JOBS" bash -lc CMD
fi

echo "[DONE] Medaka Round 2 complete. Output: $OUT_ROOT"
