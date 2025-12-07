#!/usr/bin/env bash
set -euo pipefail

########################################
# 用户参数（按需修改）
########################################
# 输入：host_deconv_out_5 目录，每个样本在子文件夹中
FASTQ_ROOT="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/host_deconv_out_5"
# 统一参考序列
REF_FASTA="/Users/jck/Desktop/ref/gold/reordered_1824_3221_reference_A2763.fasta"
# 输出根目录
OUT_ROOT="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r1"

THREADS_TOTAL="$( (command -v sysctl >/dev/null && sysctl -n hw.logicalcpu) || (command -v nproc >/dev/null && nproc) || echo 8 )"
THREADS_PER_JOB=4
JOBS=$(( THREADS_TOTAL / THREADS_PER_JOB )); [ "$JOBS" -lt 1 ] && JOBS=1

MM2_PRESET="map-ont"
MIN_MAPQ=20

########################################
# 目录与依赖
########################################
# 9.1 专用子目录
MINIMAP_DIR="${OUT_ROOT}/minimap2_align"
INDEX_DIR="${MINIMAP_DIR}/indices"
LOG_DIR="${MINIMAP_DIR}/logs"
TMP_DIR="${MINIMAP_DIR}/tmp"
BAM_DIR="${MINIMAP_DIR}/bam"
STAT_DIR="${MINIMAP_DIR}/stats"

mkdir -p "${INDEX_DIR}" "${LOG_DIR}" "${TMP_DIR}" "${BAM_DIR}" "${STAT_DIR}"

for exe in minimap2 samtools awk; do
  command -v "$exe" >/dev/null || { echo "[ERROR] 依赖未找到：$exe"; exit 1; }
done

########################################
# 准备统一参考序列索引
########################################
echo "[$(date)] Preparing FASTA index for ${REF_FASTA} ..."
[ -f "${REF_FASTA}.fai" ] || samtools faidx "${REF_FASTA}"

# 构建 minimap2 索引
REF_MMI="${INDEX_DIR}/$(basename "${REF_FASTA}" .fasta).mmi"
if [ ! -f "${REF_MMI}" ]; then
  echo "[$(date)] Building minimap2 index: ${REF_MMI}"
  minimap2 -d "${REF_MMI}" "${REF_FASTA}"
fi

########################################
# 扫描样本目录，生成样本列表
########################################
MAP_LIST="${MINIMAP_DIR}/sample_list.tsv"
: > "${MAP_LIST}"

echo "[$(date)] Scanning for samples in ${FASTQ_ROOT} ..."
shopt -s nullglob
for SAMPLE_DIR in "${FASTQ_ROOT}"/*_subsampled.trimmed_filtered.porechop; do
  [ -d "${SAMPLE_DIR}" ] || continue
  DIRNAME="$(basename "${SAMPLE_DIR}")"
  # 提取样本ID: 10090_subsampled.trimmed_filtered.porechop -> 10090
  SAMPLE_ID="${DIRNAME%%_subsampled.trimmed_filtered.porechop}"
  
  # 查找 fastq 文件
  FASTQ_FILE="${SAMPLE_DIR}/${DIRNAME}.viral_enriched.unmasked.fastq.gz"
  if [ -f "${FASTQ_FILE}" ]; then
    echo -e "${SAMPLE_ID}\t${FASTQ_FILE}" >> "${MAP_LIST}"
  else
    echo "[WARN] 找不到 fastq: ${FASTQ_FILE}"
  fi
done

N_TASKS="$(wc -l < "${MAP_LIST}" | tr -d ' ')"
echo "[$(date)] 将处理 ${N_TASKS} 个样本"

if [ "${N_TASKS}" -eq 0 ]; then
  echo "[ERROR] 未找到任何样本，请检查路径"
  exit 1
fi

########################################
# 单样本处理脚本（统一参考序列映射）
########################################
WORKER="${MINIMAP_DIR}/_r1_worker.sh"
cat > "${WORKER}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

sample="$1"          # e.g. 10090
fastq_file="$2"      # 完整路径
REF_FASTA="$3"
REF_MMI="$4"
BAM_DIR="$5"
STAT_DIR="$6"
LOG_FILE="$7"
THREADS="$8"
MM2_PRESET="$9"
MIN_MAPQ="${10}"

# 让脚本从一开始就把 stdout/stderr 都写到样本日志里，并打开命令回显
mkdir -p "$(dirname "${LOG_FILE}")"
exec >> "${LOG_FILE}" 2>&1
set -x

# 遇错时在日志里标注行号
trap 'echo "[ERROR] ${sample} failed at line $LINENO"; exit 1' ERR

if [ ! -f "${fastq_file}" ]; then
  echo "[WARN] 找不到 fastq：${fastq_file}"
  exit 0
fi

out_bam="${BAM_DIR}/${sample}.filtered.mapq${MIN_MAPQ}.bam"
flagstat_txt="${STAT_DIR}/${sample}.flagstat.txt"
idxstats_txt="${STAT_DIR}/${sample}.idxstats.tsv"
coverage_tsv="${STAT_DIR}/${sample}.coverage.tsv"

echo "[$(date)] [${sample}] Mapping to unified reference ..."

# === 映射 + 过滤 + 排序（仅排除 unmapped+secondary；保留 supplementary） ===
minimap2 -x "${MM2_PRESET}" -a -t "${THREADS}" \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ONT" \
  "${REF_MMI}" "${fastq_file}" | \
samtools view -@ "${THREADS}" -b -F 260 -q "${MIN_MAPQ}" - | \
samtools sort -@ "${THREADS}" -m 1G -o "${out_bam}" -

samtools index -@ "${THREADS}" "${out_bam}"

echo "[$(date)] [${sample}] QC stats on filtered BAM ..."
samtools flagstat -@ "${THREADS}" "${out_bam}" > "${flagstat_txt}"
samtools idxstats "${out_bam}" > "${idxstats_txt}"
samtools coverage -q "${MIN_MAPQ}" "${out_bam}" > "${coverage_tsv}"

echo "[$(date)] [${sample}] Done."
EOS
chmod +x "${WORKER}"

########################################
# 并行运行
########################################
echo "[$(date)] Starting mapping+QC with ${JOBS} concurrent jobs; ${THREADS_PER_JOB} threads/job."

if command -v parallel >/dev/null 2>&1; then
  parallel --jobs "${JOBS}" --colsep '\t' --halt soon,fail=1 --results "${LOG_DIR}/parallel_results" --bar \
    "${WORKER} {1} {2} ${REF_FASTA} ${REF_MMI} ${BAM_DIR} ${STAT_DIR} ${LOG_DIR}/{1}.log ${THREADS_PER_JOB} ${MM2_PRESET} ${MIN_MAPQ}" \
    :::: "${MAP_LIST}"
else
  echo "[WARN] 未检测到 GNU parallel，将串行运行（建议安装 parallel 以提速）"
  while IFS=$'\t' read -r sample fastq_file; do
    "${WORKER}" "${sample}" "${fastq_file}" "${REF_FASTA}" "${REF_MMI}" \
      "${BAM_DIR}" "${STAT_DIR}" "${LOG_DIR}/${sample}.log" \
      "${THREADS_PER_JOB}" "${MM2_PRESET}" "${MIN_MAPQ}"
  done < "${MAP_LIST}"
fi

echo "[$(date)] All done. Filtered BAMs in: ${BAM_DIR}"
echo "Stats in: ${STAT_DIR} ; Logs in: ${LOG_DIR}"

