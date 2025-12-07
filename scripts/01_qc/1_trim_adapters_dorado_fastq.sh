#!/usr/bin/env bash
# Adapter trimming for ONT FASTQ with dorado trim
# Supports .fastq/.fq and .gz variants
set -Eeuo pipefail

# -------- User defaults --------
SRC_DIR="${SRC_DIR:-/path/to/input}"
OUT_DIR="${OUT_DIR:-/path/to/output}"
KIT="${KIT:-SQK-NBD114-24}"
JOBS="${JOBS:-2}"
TPJ="${TPJ:-}"
MAX_READS_TEST="${MAX_READS_TEST:-2000}"
MODE=""
FORCE=0
DRYRUN=0

usage() {
  cat <<EOF
Usage:
  $0 --test [N]              # Test mode: process first N reads (default: $MAX_READS_TEST)
  $0 --run                   # Full run: process all FASTQ files
Options:
  -i, --input  DIR           # Input directory (default: $SRC_DIR)
  -o, --output DIR           # Output directory (default: $OUT_DIR)
  -k, --kit    NAME          # Sequencing kit (default: $KIT)
  -j, --jobs   INT           # Parallel jobs (default: $JOBS)
  -t, --threads-per-job INT  # Threads per job (default: auto)
      --force                # Overwrite existing files
      --dry-run              # Print commands without executing
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test) MODE="test"; shift; [[ "${1-}" =~ ^[0-9]+$ ]] && MAX_READS_TEST="$1" && shift || true ;;
    --run) MODE="run"; shift ;;
    -i|--input) SRC_DIR="$2"; shift 2 ;;
    -o|--output) OUT_DIR="$2"; shift 2 ;;
    -k|--kit) KIT="$2"; shift 2 ;;
    -j|--jobs) JOBS="$2"; shift 2 ;;
    -t|--threads-per-job) TPJ="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$MODE" ]] && { echo "Please specify --test or --run"; usage; exit 1; }
[[ -d "$SRC_DIR" ]] || { echo "Input directory not found: $SRC_DIR"; exit 1; }
mkdir -p "$OUT_DIR" "$OUT_DIR/_logs" "$OUT_DIR/_test"

# ---- Check dorado availability ----
DORADO_AVAILABLE=0
if command -v dorado >/dev/null 2>&1; then
  DORADO_AVAILABLE=1
else
  if [[ ${DRYRUN:-0} -eq 1 ]]; then
    echo "[WARN] dorado not found, dry-run mode will use default flags"
  else
    echo "dorado not found. Please run: conda activate dorado"; exit 1
  fi
fi

# Version check (requires >= 1.1.0)
version_ge() {
  local IFS=.
  local a b i
  read -r -a a <<< "${1:-0.0.0}"
  read -r -a b <<< "${2:-0.0.0}"
  for ((i=${#a[@]}; i<3; i++)); do a[i]=0; done
  for ((i=${#b[@]}; i<3; i++)); do b[i]=0; done
  for i in 0 1 2; do
    if ((10#${a[i]} > 10#${b[i]})); then return 0; fi
    if ((10#${a[i]} < 10#${b[i]})); then return 1; fi
  done
  return 0
}

if [[ $DORADO_AVAILABLE -eq 1 ]]; then
  DORADO_VER_LINE="$(dorado --version 2>&1 | tail -n 1)"
  DORADO_SEMVER="$(printf '%s' "$DORADO_VER_LINE" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true)"
  echo "[INFO] dorado version: ${DORADO_SEMVER:-unknown} (requires >= 1.1.0)"
  if [[ -n "$DORADO_SEMVER" ]]; then
    if ! version_ge "$DORADO_SEMVER" "1.1.0"; then
      echo "[ERROR] dorado version too low: $DORADO_SEMVER, need at least 1.1.0"; exit 1
    fi
  else
    echo "[WARN] Cannot parse dorado version: $DORADO_VER_LINE"
  fi
else
  echo "[INFO] dorado not ready (dry-run mode), skipping version check"
fi

# CPU & threading setup
CORES="$( (sysctl -n hw.ncpu 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 8) | awk '{print $1}')"
if [[ -z "$TPJ" ]]; then
  TPJ=$(( CORES / JOBS ))
  (( TPJ < 1 )) && TPJ=1
fi
echo "[INFO] CPU cores=$CORES; jobs=$JOBS; threads-per-job=$TPJ"

# Compressor selection
if command -v pigz >/dev/null 2>&1; then
  COMPRESS_CMD=(pigz -p "$TPJ" -1)
  COMPRESS="pigz -p ${TPJ} -1"
else
  COMPRESS_CMD=(gzip -1)
  COMPRESS="gzip -1"
fi
echo "[INFO] compressor: ${COMPRESS}"

# ---- Detect dorado flags ----
KIT_FLAG=""
EMIT_FASTQ_FLAG=""
NO_PRIMERS_FLAG=""
MAX_READS_FLAG=""

detect_dorado_flags() {
  if [[ $DORADO_AVAILABLE -eq 0 ]]; then
    KIT_FLAG="--sequencing-kit"
    EMIT_FASTQ_FLAG="--emit-fastq"
    NO_PRIMERS_FLAG="--no-trim-primers"
    MAX_READS_FLAG="-n"
    echo "[INFO] dry-run: using default flags"
    return 0
  fi

  local help
  if ! help="$(dorado trim --help 2>&1)"; then
    echo "[ERROR] Cannot get 'dorado trim --help' output"; exit 1
  fi

  if grep -q -- "--sequencing-kit" <<< "$help"; then
    KIT_FLAG="--sequencing-kit"
  elif grep -q -- "--kit-name" <<< "$help"; then
    KIT_FLAG="--kit-name"
  else
    echo "[ERROR] Kit flag not found in dorado trim --help"; exit 1
  fi

  [[ $(grep -c -- "--emit-fastq" <<< "$help") -gt 0 ]] && EMIT_FASTQ_FLAG="--emit-fastq"
  [[ $(grep -c -- "--no-trim-primers" <<< "$help") -gt 0 ]] && NO_PRIMERS_FLAG="--no-trim-primers"
  
  if grep -qE '(^|[,[:space:]])-n([,[:space:]]|$)' <<< "$help"; then
    MAX_READS_FLAG="-n"
  elif grep -q -- "--max-reads" <<< "$help"; then
    MAX_READS_FLAG="--max-reads"
  fi
}
detect_dorado_flags

# ---- Helper functions ----
list_inputs() {
  find "$SRC_DIR" -type f \
    \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fq" \) \
    -not -name "*.trimmed.fastq.gz" -not -name "*.trimmed.fastq" \
    -not -path '*/.*' -not -name '._*' \
    -not -path "$OUT_DIR/_logs/*" -not -path "$OUT_DIR/_test/*" \
    -print0
}

stem_from() {
  local base
  base="$(basename "$1")"
  case "$base" in
    *.fastq.gz) printf '%s' "${base%.fastq.gz}" ;;
    *.fq.gz)    printf '%s' "${base%.fq.gz}" ;;
    *.fastq)    printf '%s' "${base%.fastq}" ;;
    *.fq)       printf '%s' "${base%.fq}" ;;
    *)          printf '%s' "${base%.*}" ;;
  esac
}

first_fastq() {
  while IFS= read -r -d '' f; do
    printf '%s' "$f"
    return 0
  done < <(list_inputs)
  return 1
}

process_file() {
  local in="$1"
  local base stem out log
  base="$(basename "$in")"
  stem="$(stem_from "$in")"
  out="$OUT_DIR/${stem}.trimmed.fastq.gz"
  log="$OUT_DIR/_logs/${stem}.log"

  if [[ $FORCE -eq 0 && -s "$out" ]]; then
    echo "[SKIP] Already exists: $out"
    return 0
  fi

  local dorado_cmd=(dorado trim "$in" "$KIT_FLAG" "$KIT" -t "$TPJ")
  [[ -n "$EMIT_FASTQ_FLAG" ]] && dorado_cmd+=("$EMIT_FASTQ_FLAG")
  [[ -n "$NO_PRIMERS_FLAG" ]] && dorado_cmd+=("$NO_PRIMERS_FLAG")

  echo "[RUN ] $base -> $(basename "$out")"
  if [[ $DRYRUN -eq 1 ]]; then
    echo "CMD: ${dorado_cmd[*]} 2> '$log' | ${COMPRESS} > '$out'"
    return 0
  fi

  set -o pipefail
  "${dorado_cmd[@]}" 2> "$log" | "${COMPRESS_CMD[@]}" > "$out"
  echo "[DONE] $base"
}

# ---- Main ----
if [[ "$MODE" == "test" ]]; then
  TEST_FILE="$(first_fastq)"
  [[ -z "$TEST_FILE" ]] && { echo "No FASTQ files found"; exit 1; }

  base="$(basename "$TEST_FILE")"
  stem="$(stem_from "$TEST_FILE")"
  out="$OUT_DIR/_test/${stem}.test.trimmed.fastq.gz"
  log="$OUT_DIR/_logs/${stem}.test.log"

  dorado_cmd=(dorado trim "$TEST_FILE" "$KIT_FLAG" "$KIT" -t "$TPJ")
  [[ -n "$EMIT_FASTQ_FLAG" ]] && dorado_cmd+=("$EMIT_FASTQ_FLAG")
  [[ -n "$NO_PRIMERS_FLAG" ]] && dorado_cmd+=("$NO_PRIMERS_FLAG")
  [[ -n "$MAX_READS_FLAG" ]] && dorado_cmd+=("$MAX_READS_FLAG" "$MAX_READS_TEST")

  echo "[TEST] $base -> $(basename "$out") (max_reads=$MAX_READS_TEST)"
  if [[ $DRYRUN -eq 1 ]]; then
    echo "CMD: ${dorado_cmd[*]} 2> \"$log\" | ${COMPRESS} > \"$out\""; exit 0
  fi

  set -o pipefail
  if [[ -n "$MAX_READS_FLAG" ]]; then
    "${dorado_cmd[@]}" 2> "$log" | "${COMPRESS_CMD[@]}" > "$out"
  else
    "${dorado_cmd[@]}" 2> "$log" | head -n $((MAX_READS_TEST*4)) | "${COMPRESS_CMD[@]}" > "$out"
  fi
  echo "[TEST] Done: $out"
  exit 0
fi

if [[ "$MODE" == "run" ]]; then
  echo "[INFO] Starting full trimming ..."
  cnt=0
  while IFS= read -r -d '' f; do
    process_file "$f" &
    cnt=$((cnt + 1))
    if (( cnt % JOBS == 0 )); then
      wait
    fi
  done < <(list_inputs)
  wait
  echo "[INFO] All done. Output: $OUT_DIR"
fi
