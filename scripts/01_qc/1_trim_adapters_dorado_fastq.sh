#!/usr/bin/env bash
# Adapter trimming for ONT FASTQ(.gz) with dorado trim
# Now also supports plain .fastq/.fq (and .gz variants)
# Fixed: robust version detection & safe first_fastq (no SIGPIPE)
set -Eeuo pipefail

# -------- User defaults --------
SRC_DIR="${SRC_DIR:-/Volumes/OWC/Analysis/trimmed}"
OUT_DIR="${OUT_DIR:-/Volumes/OWC/Analysis/trimmed_out}"
KIT="${KIT:-SQK-NBD114-24}"              # 你的试剂盒
JOBS="${JOBS:-2}"                         # 并发作业数
TPJ="${TPJ:-}"                            # 每作业线程数（为空则自动均分）
MAX_READS_TEST="${MAX_READS_TEST:-2000}"  # 测试模式最多处理的 reads
MODE=""
FORCE=0
DRYRUN=0

usage() {
  cat <<EOF
Usage:
  $0 --test [N]           # 测试模式：只跑一个文件的前 N 条 reads（默认 $MAX_READS_TEST）
  $0 --run                # 正式运行：处理目录下所有 .fastq(.gz) / .fq(.gz)
Options:
  -i, --input  DIR        # 输入目录（默认 $SRC_DIR）
  -o, --output DIR        # 输出目录（默认 $OUT_DIR）
  -k, --kit    NAME       # sequencing kit（默认 $KIT）
  -j, --jobs   INT        # 并发作业数（默认 $JOBS）
  -t, --threads-per-job INT  # 每作业线程数（默认自动=CPU核数/JOBS）
      --force             # 覆盖已存在的输出文件
      --dry-run           # 仅打印命令，不实际运行
  -h, --help
EOF
}

# ---- parse args ----
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

[[ -z "$MODE" ]] && { echo "请指定 --test 或 --run"; usage; exit 1; }
[[ -d "$SRC_DIR" ]] || { echo "输入目录不存在: $SRC_DIR"; exit 1; }
mkdir -p "$OUT_DIR" "$OUT_DIR/_logs" "$OUT_DIR/_test"

# ---- env checks ----
DORADO_AVAILABLE=0
if command -v dorado >/dev/null 2>&1; then
  DORADO_AVAILABLE=1
else
  if [[ ${DRYRUN:-0} -eq 1 ]]; then
    echo "[WARN ] 未找到 dorado，但处于 --dry-run 模式，将仅打印命令"
  else
    echo "未找到 dorado，请先 'conda activate dorado'"; exit 1
  fi
fi

# 版本解析与校验 (要求 >= 1.1.0) —— 仅在 dorado 可用时检查
version_ge() {
  # 语义化版本比较：a >= b 返回 0
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
  echo "[INFO] dorado version: ${DORADO_SEMVER:-unknown} (要求 >= 1.1.0)"
  if [[ -n "$DORADO_SEMVER" ]]; then
    if ! version_ge "$DORADO_SEMVER" "1.1.0"; then
      echo "[ERROR] dorado 版本过低：$DORADO_SEMVER，至少需要 1.1.0"; exit 1
    fi
  else
    echo "[WARN ] 无法解析 dorado 版本字符串：$DORADO_VER_LINE"
  fi
else
  echo "[INFO] dorado 未就绪（dry-run 模式），跳过版本检查"
fi

# CPU cores & threading
CORES="$( (sysctl -n hw.ncpu 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 8) | awk '{print $1}')"
if [[ -z "$TPJ" ]]; then
  TPJ=$(( CORES / JOBS ))
  (( TPJ < 1 )) && TPJ=1
fi
echo "[INFO] CPU cores=$CORES; jobs=$JOBS; threads-per-job=$TPJ"

# Compressor
if command -v pigz >/dev/null 2>&1; then
  COMPRESS_CMD=(pigz -p "$TPJ" -1)
  COMPRESS="pigz -p ${TPJ} -1"
else
  COMPRESS_CMD=(gzip -1)
  COMPRESS="gzip -1"
fi
echo "[INFO] compressor: ${COMPRESS}"

# ---- detect dorado flags (兼容不同版本命名) ----
KIT_FLAG=""
EMIT_FASTQ_FLAG=""
NO_PRIMERS_FLAG=""
MAX_READS_FLAG=""
detect_dorado_flags() {
  if [[ $DORADO_AVAILABLE -eq 0 ]]; then
    # dry-run 回退：采用默认/常见参数名
    KIT_FLAG="--sequencing-kit"
    EMIT_FASTQ_FLAG="--emit-fastq"
    NO_PRIMERS_FLAG="--no-trim-primers"
    MAX_READS_FLAG="-n"
    echo "[INFO] dry-run：使用默认参数名 (未探测)"
    return 0
  fi

  local help
  if ! help="$(dorado trim --help 2>&1)"; then
    echo "[ERROR] 无法获取 'dorado trim --help' 输出"; exit 1
  fi

  if grep -q -- "--sequencing-kit" <<< "$help"; then
    KIT_FLAG="--sequencing-kit"
  elif grep -q -- "--kit-name" <<< "$help"; then
    KIT_FLAG="--kit-name"
  else
    echo "[ERROR] 未在 'dorado trim --help' 中找到 kit 相关选项 (--sequencing-kit/--kit-name)"; exit 1
  fi

  if grep -q -- "--emit-fastq" <<< "$help"; then
    EMIT_FASTQ_FLAG="--emit-fastq"
  else
    EMIT_FASTQ_FLAG=""
  fi

  if grep -q -- "--no-trim-primers" <<< "$help"; then
    NO_PRIMERS_FLAG="--no-trim-primers"
  else
    NO_PRIMERS_FLAG=""
  fi

  if grep -qE '(^|[,[:space:]])-n([,[:space:]]|$)' <<< "$help"; then
    MAX_READS_FLAG="-n"
  elif grep -q -- "--max-reads" <<< "$help"; then
    MAX_READS_FLAG="--max-reads"
  else
    MAX_READS_FLAG=""
  fi
}
detect_dorado_flags

# ---- helpers ----

# 统一的输入文件枚举：支持 .fastq(.gz) 和 .fq(.gz)，并排除 trimmed 结果与隐藏/临时文件
list_inputs() {
  find "$SRC_DIR" -type f \
    \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fq" \) \
    -not -name "*.trimmed.fastq.gz" -not -name "*.trimmed.fastq" \
    -not -path '*/.*' -not -name '._*' \
    -not -path "$OUT_DIR/_logs/*" -not -path "$OUT_DIR/_test/*" \
    -print0
}

# 规范化去除扩展名，得到 stem（样本名）
stem_from() {
  local base
  base="$(basename "$1")"
  case "$base" in
    *.fastq.gz) printf '%s' "${base%.fastq.gz}" ;;
    *.fq.gz)    printf '%s' "${base%.fq.gz}" ;;
    *.fastq)    printf '%s' "${base%.fastq}" ;;
    *.fq)       printf '%s' "${base%.fq}" ;;
    *)          printf '%s' "${base%.*}" ;;  # fallback
  esac
}

first_fastq() {
  # 返回第一个非隐藏 fastq/fq（避免 SIGPIPE）
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
    echo "[SKIP] 已存在: $out"
    return 0
  fi

  # 组装 dorado 命令 (数组形式以避免 eval)
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

# ---- main ----
if [[ "$MODE" == "test" ]]; then
  TEST_FILE="$(first_fastq)"
  [[ -z "$TEST_FILE" ]] && { echo "未找到 *.fastq(.gz) 或 *.fq(.gz)"; exit 1; }

  base="$(basename "$TEST_FILE")"
  stem="$(stem_from "$TEST_FILE")"
  out="$OUT_DIR/_test/${stem}.test.trimmed.fastq.gz"
  log="$OUT_DIR/_logs/${stem}.test.log"

  dorado_cmd=(dorado trim "$TEST_FILE" "$KIT_FLAG" "$KIT" -t "$TPJ")
  [[ -n "$EMIT_FASTQ_FLAG" ]] && dorado_cmd+=("$EMIT_FASTQ_FLAG")
  [[ -n "$NO_PRIMERS_FLAG" ]] && dorado_cmd+=("$NO_PRIMERS_FLAG")
  if [[ -n "$MAX_READS_FLAG" ]]; then
    dorado_cmd+=("$MAX_READS_FLAG" "$MAX_READS_TEST")
  fi

  echo "[TEST] $base -> $(basename "$out") (max_reads=$MAX_READS_TEST)"
  if [[ $DRYRUN -eq 1 ]]; then
    if [[ -n "$MAX_READS_FLAG" ]]; then
      echo "CMD: ${dorado_cmd[*]} 2> \"$log\" | ${COMPRESS} > \"$out\""; exit 0
    else
      echo "CMD: ${dorado_cmd[*]} 2> \"$log\" | head -n $((MAX_READS_TEST*4)) | ${COMPRESS} > \"$out\""; exit 0
    fi
  fi

  set -o pipefail
  if [[ -n "$MAX_READS_FLAG" ]]; then
    "${dorado_cmd[@]}" 2> "$log" | "${COMPRESS_CMD[@]}" > "$out"
  else
    # 回退：若无 max-reads 选项，则用 head 近似限制到 N 条 FASTQ 记录（每条 4 行）
    "${dorado_cmd[@]}" 2> "$log" | head -n $((MAX_READS_TEST*4)) | "${COMPRESS_CMD[@]}" > "$out"
  fi
  echo "[TEST] 完成：$out"
  exit 0
fi

if [[ "$MODE" == "run" ]]; then
  echo "[INFO] 开始全量修剪 ..."
  cnt=0
  while IFS= read -r -d '' f; do
    process_file "$f" &
    cnt=$((cnt + 1))
    if (( cnt % JOBS == 0 )); then
      wait
    fi
  done < <(list_inputs)
  wait
  echo "[INFO] 全部完成。输出目录: $OUT_DIR"
fi
