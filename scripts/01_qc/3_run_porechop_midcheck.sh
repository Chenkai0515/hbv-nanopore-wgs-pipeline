#!/usr/bin/env bash
set -euo pipefail

# ========== 用户可改参数 ==========
# 输入与输出目录（与你提供的路径一致）
INPUT_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_filter_2"
OUTPUT_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_porechop_3"

# 资源配置：总核数与每个作业的线程数（自动计算并发作业数）
TOTAL_THREADS=16
THREADS_PER_JOB=8     # 建议设为 8（并发 2 个作业），更高并发请酌情调整内存

# 仅检查中间接头（不拆分/不丢弃）？0/1
CHECK_ONLY=0

# 是否使用 Porechop_ABI 做 ab initio 接头推断（更适合 Kit14/114 等新试剂盒）？0/1
USE_ABI=1

# 阈值可调（维持 Porechop 默认值最稳妥）
MIDDLE_THRESHOLD=85   # 中间接头识别阈值（越低越敏感）
END_THRESHOLD=75      # 端部接头修剪阈值
CHECK_READS=10000     # 参与“识别存在哪些接头”的采样 reads 数（默认 10000）

# ========== 环境设置（macOS 友好） ==========
umask 002
ulimit -n 8192 || true
export LC_ALL=C

LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_DIR="${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${REPORT_DIR}"

# 选择可执行文件：优先用用户指定的 ABI 模式，否则退回到标准 Porechop
if [[ "${USE_ABI}" -eq 1 ]]; then
  if command -v porechop_abi >/dev/null 2>&1; then
    PORECHOP_BIN="porechop_abi"
    ABI_FLAGS=( -abi -ddb )   # ab initio + 丢弃内置数据库，避免老数据库误报/加速
  else
    echo "WARN: 未找到 porechop_abi，将退回使用 porechop。" >&2
    PORECHOP_BIN="porechop"
    ABI_FLAGS=()
  fi
else
  PORECHOP_BIN="porechop"
  ABI_FLAGS=()
fi

if ! command -v "${PORECHOP_BIN}" >/dev/null 2>&1; then
  echo "ERROR: 未找到 ${PORECHOP_BIN}，请先用 conda/mamba 安装："
  echo "       mamba/conda install -c bioconda porechop    # 或安装 porechop_abi"
  exit 1
fi

# 计算并发作业数
MAX_JOBS=$(( TOTAL_THREADS / THREADS_PER_JOB ))
[[ "${MAX_JOBS}" -lt 1 ]] && MAX_JOBS=1

echo "使用 ${PORECHOP_BIN}；并发 ${MAX_JOBS} 作业，每作业 ${THREADS_PER_JOB} 线程。"
echo "输入：${INPUT_DIR}"
echo "输出：${OUTPUT_DIR}"
echo -e "-----------------------------------------\n"

# 汇总表
SUMMARY_TSV="${REPORT_DIR}/porechop_summary.tsv"
echo -e "sample\tmode(split_or_discard_or_nosplit)\ttotal_reads\tstart_trimmed\tend_trimmed\tmiddle_events" > "${SUMMARY_TSV}"

# 解析日志以提取统计
parse_log() {
  local log="$1"; local sample="$2"
  local total="" start_trim="" end_trim="" middle_events="" mode=""

  start_trim=$(grep -E 'reads had adapters trimmed from their start' "$log" | tail -n1 | awk '{gsub(",",""); print $1}' || true)
  total=$(grep -E 'reads had adapters trimmed from their start' "$log" | tail -n1 | awk '{gsub(",",""); print $3}' || true)
  end_trim=$(grep -E 'reads had adapters trimmed from their end' "$log"   | tail -n1 | awk '{gsub(",",""); print $1}' || true)

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

# 单文件处理
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
    # 仅检查中间接头，不做分割/丢弃（仍会输出与输入等长的 reads；若不想产出文件，可改 -o /dev/null）
    "${PORECHOP_BIN}" "${ABI_FLAGS[@]}" "${common_args[@]}" --no_split -o "$out_f" --format fastq.gz >"$log_f" 2>&1
  else
    # 推荐：对中间接头 reads 直接丢弃（在已 demux 的条码数据上更安全）:contentReference[oaicite:3]{index=3}
    "${PORECHOP_BIN}" "${ABI_FLAGS[@]}" "${common_args[@]}" --discard_middle -o "$out_f" --format fastq.gz >"$log_f" 2>&1
  fi

  parse_log "$log_f" "$sample"
}

# 遍历输入文件并并发执行（兼容 macOS Bash 3.x）
active=0
while IFS= read -r -d '' f; do
  process_one "$f" &
  active=$((active+1))
  if (( active % MAX_JOBS == 0 )); then
    wait
  fi
done < <(find "${INPUT_DIR}" -type f -name '*.fastq.gz' -print0)

wait

echo -e "\n✔ 完成。"
echo "  修剪后的 FASTQ： ${OUTPUT_DIR}"
echo "  逐样本日志：      ${LOG_DIR}"
echo "  汇总统计表：      ${SUMMARY_TSV}"
