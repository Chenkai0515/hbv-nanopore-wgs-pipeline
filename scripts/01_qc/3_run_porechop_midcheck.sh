#!/usr/bin/env bash
# Mid-adapter removal using Porechop/Porechop_ABI
set -euo pipefail

# ========== Configuration ==========
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

INPUT_DIR="${INPUT_DIR:-${PROJECT_DIR}/fastq_filter_2}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/fastq_porechop_3}"

TOTAL_THREADS="${TOTAL_THREADS:-16}"
THREADS_PER_JOB="${THREADS_PER_JOB:-8}"

CHECK_ONLY="${CHECK_ONLY:-0}"      # 1=check only, 0=discard middle adapters
USE_ABI="${USE_ABI:-1}"            # 1=use porechop_abi, 0=use porechop

MIDDLE_THRESHOLD="${MIDDLE_THRESHOLD:-85}"
END_THRESHOLD="${END_THRESHOLD:-75}"
CHECK_READS="${CHECK_READS:-10000}"

# ========== Environment setup ==========
umask 002
ulimit -n 8192 || true
export LC_ALL=C

LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_DIR="${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${REPORT_DIR}"

# Select executable: prefer porechop_abi if USE_ABI=1
if [[ "${USE_ABI}" -eq 1 ]]; then
  if command -v porechop_abi >/dev/null 2>&1; then
    PORECHOP_BIN="porechop_abi"
    ABI_FLAGS=( -abi -ddb )
  else
    echo "WARN: porechop_abi not found, falling back to porechop" >&2
    PORECHOP_BIN="porechop"
    ABI_FLAGS=()
  fi
else
  PORECHOP_BIN="porechop"
  ABI_FLAGS=()
fi

if ! command -v "${PORECHOP_BIN}" >/dev/null 2>&1; then
  echo "ERROR: ${PORECHOP_BIN} not found. Install with:"
  echo "       mamba install -c bioconda porechop  # or porechop_abi"
  exit 1
fi

MAX_JOBS=$(( TOTAL_THREADS / THREADS_PER_JOB ))
[[ "${MAX_JOBS}" -lt 1 ]] && MAX_JOBS=1

echo "Using ${PORECHOP_BIN}; ${MAX_JOBS} jobs, ${THREADS_PER_JOB} threads/job"
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "-----------------------------------------"

# Summary table
SUMMARY_TSV="${REPORT_DIR}/porechop_summary.tsv"
echo -e "sample\tmode\ttotal_reads\tstart_trimmed\tend_trimmed\tmiddle_events" > "${SUMMARY_TSV}"

# Parse log to extract statistics
parse_log() {
  local log="$1"; local sample="$2"
  local total="" start_trim="" end_trim="" middle_events="" mode=""

  start_trim=$(grep -E 'reads had adapters trimmed from their start' "$log" | tail -n1 | awk '{gsub(",",""); print $1}' || true)
  total=$(grep -E 'reads had adapters trimmed from their start' "$log" | tail -n1 | awk '{gsub(",",""); print $3}' || true)
  end_trim=$(grep -E 'reads had adapters trimmed from their end' "$log" | tail -n1 | awk '{gsub(",",""); print $1}' || true)

  local mid_line
  mid_line=$(grep -E 'reads were (split|discarded) based on middle adapters' "$log" | tail -n1 || true)
  if [[ -n "${mid_line}" ]]; then
    middle_events=$(echo "$mid_line" | sed -E 's/^ *([0-9,]+) \/ ([0-9,]+).*/\1/' | tr -d ',')
    [[ -z "${total}" ]] && total=$(echo "$mid_line" | sed -E 's/^ *([0-9,]+) \/ ([0-9,]+).*/\2/' | tr -d ',')
    mode=$(echo "$mid_line" | sed -nE 's/.*reads were (split|discarded) based on middle adapters.*/\1/p')
  else
    mode="nosplit"
    middle_events="NA"
  fi

  [[ -z "${start_trim}" ]] && start_trim="0"
  [[ -z "${end_trim}"   ]] && end_trim="0"
  [[ -z "${total}"      ]] && total="NA"
  [[ -z "${mode}"       ]] && mode="split"

  echo -e "${sample}\t${mode}\t${total}\t${start_trim}\t${end_trim}\t${middle_events}" >> "${SUMMARY_TSV}"
}

# Process single file
process_one() {
  local in_f="$1"
  local base; base=$(basename "$in_f")
  local sample="${base%.fastq.gz}"
  local out_f="${OUTPUT_DIR}/${sample}.porechop.fastq.gz"
  local log_f="${LOG_DIR}/${sample}.porechop.log"

  common_args=( -i "$in_f" -t "${THREADS_PER_JOB}" -v 1 \
                --end_threshold "${END_THRESHOLD}" \
                --middle_threshold "${MIDDLE_THRESHOLD}" \
                --check_reads "${CHECK_READS}" )

  if [[ "${CHECK_ONLY}" -eq 1 ]]; then
    "${PORECHOP_BIN}" "${ABI_FLAGS[@]}" "${common_args[@]}" --no_split -o "$out_f" --format fastq.gz >"$log_f" 2>&1
  else
    "${PORECHOP_BIN}" "${ABI_FLAGS[@]}" "${common_args[@]}" --discard_middle -o "$out_f" --format fastq.gz >"$log_f" 2>&1
  fi

  parse_log "$log_f" "$sample"
}

# Parallel processing
active=0
while IFS= read -r -d '' f; do
  process_one "$f" &
  active=$((active+1))
  if (( active % MAX_JOBS == 0 )); then
    wait
  fi
done < <(find "${INPUT_DIR}" -type f -name '*.fastq.gz' -print0)

wait

echo ""
echo "Done."
echo "  Trimmed FASTQ: ${OUTPUT_DIR}"
echo "  Logs: ${LOG_DIR}"
echo "  Summary: ${SUMMARY_TSV}"
