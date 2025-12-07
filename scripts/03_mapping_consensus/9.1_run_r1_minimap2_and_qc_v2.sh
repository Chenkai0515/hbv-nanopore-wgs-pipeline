#!/usr/bin/env bash
# Round 1: Map reads to unified HBV reference for Medaka consensus
set -euo pipefail

########################################
# Path Configuration
########################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

FASTQ_ROOT="${FASTQ_ROOT:-${PROJECT_DIR}/host_deconv_out_5}"
REF_FASTA="${REF_FASTA:-${HBV_REF:-/path/to/reference.fasta}}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_DIR}/Medaka_consensus_8/r1}"

THREADS_TOTAL="$( (command -v sysctl >/dev/null && sysctl -n hw.logicalcpu) || (command -v nproc >/dev/null && nproc) || echo 8 )"
THREADS_PER_JOB=4
JOBS=$(( THREADS_TOTAL / THREADS_PER_JOB )); [ "$JOBS" -lt 1 ] && JOBS=1

MM2_PRESET="map-ont"
MIN_MAPQ=20

########################################
# Directory setup & dependencies
########################################
MINIMAP_DIR="${OUT_ROOT}/minimap2_align"
INDEX_DIR="${MINIMAP_DIR}/indices"
LOG_DIR="${MINIMAP_DIR}/logs"
TMP_DIR="${MINIMAP_DIR}/tmp"
BAM_DIR="${MINIMAP_DIR}/bam"
STAT_DIR="${MINIMAP_DIR}/stats"

mkdir -p "${INDEX_DIR}" "${LOG_DIR}" "${TMP_DIR}" "${BAM_DIR}" "${STAT_DIR}"

for exe in minimap2 samtools awk; do
  command -v "$exe" >/dev/null || { echo "[ERROR] Dependency not found: $exe"; exit 1; }
done

########################################
# Prepare reference index
########################################
echo "[$(date)] Preparing FASTA index for ${REF_FASTA} ..."
[ -f "${REF_FASTA}.fai" ] || samtools faidx "${REF_FASTA}"

REF_MMI="${INDEX_DIR}/$(basename "${REF_FASTA}" .fasta).mmi"
if [ ! -f "${REF_MMI}" ]; then
  echo "[$(date)] Building minimap2 index: ${REF_MMI}"
  minimap2 -d "${REF_MMI}" "${REF_FASTA}"
fi

########################################
# Scan samples
########################################
MAP_LIST="${MINIMAP_DIR}/sample_list.tsv"
: > "${MAP_LIST}"

echo "[$(date)] Scanning for samples in ${FASTQ_ROOT} ..."
shopt -s nullglob
for SAMPLE_DIR in "${FASTQ_ROOT}"/*_subsampled.trimmed_filtered.porechop; do
  [ -d "${SAMPLE_DIR}" ] || continue
  DIRNAME="$(basename "${SAMPLE_DIR}")"
  SAMPLE_ID="${DIRNAME%%_subsampled.trimmed_filtered.porechop}"
  
  FASTQ_FILE="${SAMPLE_DIR}/${DIRNAME}.viral_enriched.unmasked.fastq.gz"
  if [ -f "${FASTQ_FILE}" ]; then
    echo -e "${SAMPLE_ID}\t${FASTQ_FILE}" >> "${MAP_LIST}"
  else
    echo "[WARN] FASTQ not found: ${FASTQ_FILE}"
  fi
done

N_TASKS="$(wc -l < "${MAP_LIST}" | tr -d ' ')"
echo "[$(date)] Found ${N_TASKS} samples"

if [ "${N_TASKS}" -eq 0 ]; then
  echo "[ERROR] No samples found"
  exit 1
fi

########################################
# Worker script
########################################
WORKER="${MINIMAP_DIR}/_r1_worker.sh"
cat > "${WORKER}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

sample="$1"
fastq_file="$2"
REF_FASTA="$3"
REF_MMI="$4"
BAM_DIR="$5"
STAT_DIR="$6"
LOG_FILE="$7"
THREADS="$8"
MM2_PRESET="$9"
MIN_MAPQ="${10}"

mkdir -p "$(dirname "${LOG_FILE}")"
exec >> "${LOG_FILE}" 2>&1
set -x

trap 'echo "[ERROR] ${sample} failed at line $LINENO"; exit 1' ERR

if [ ! -f "${fastq_file}" ]; then
  echo "[WARN] FASTQ not found: ${fastq_file}"
  exit 0
fi

out_bam="${BAM_DIR}/${sample}.filtered.mapq${MIN_MAPQ}.bam"
flagstat_txt="${STAT_DIR}/${sample}.flagstat.txt"
idxstats_txt="${STAT_DIR}/${sample}.idxstats.tsv"
coverage_tsv="${STAT_DIR}/${sample}.coverage.tsv"

echo "[$(date)] [${sample}] Mapping to unified reference ..."

minimap2 -x "${MM2_PRESET}" -a -t "${THREADS}" \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ONT" \
  "${REF_MMI}" "${fastq_file}" | \
samtools view -@ "${THREADS}" -b -F 260 -q "${MIN_MAPQ}" - | \
samtools sort -@ "${THREADS}" -m 1G -o "${out_bam}" -

samtools index -@ "${THREADS}" "${out_bam}"

echo "[$(date)] [${sample}] QC stats ..."
samtools flagstat -@ "${THREADS}" "${out_bam}" > "${flagstat_txt}"
samtools idxstats "${out_bam}" > "${idxstats_txt}"
samtools coverage -q "${MIN_MAPQ}" "${out_bam}" > "${coverage_tsv}"

echo "[$(date)] [${sample}] Done."
EOS
chmod +x "${WORKER}"

########################################
# Run in parallel
########################################
echo "[$(date)] Starting mapping with ${JOBS} jobs, ${THREADS_PER_JOB} threads/job"

if command -v parallel >/dev/null 2>&1; then
  parallel --jobs "${JOBS}" --colsep '\t' --halt soon,fail=1 --results "${LOG_DIR}/parallel_results" --bar \
    "${WORKER} {1} {2} ${REF_FASTA} ${REF_MMI} ${BAM_DIR} ${STAT_DIR} ${LOG_DIR}/{1}.log ${THREADS_PER_JOB} ${MM2_PRESET} ${MIN_MAPQ}" \
    :::: "${MAP_LIST}"
else
  echo "[WARN] GNU parallel not found, running sequentially"
  while IFS=$'\t' read -r sample fastq_file; do
    "${WORKER}" "${sample}" "${fastq_file}" "${REF_FASTA}" "${REF_MMI}" \
      "${BAM_DIR}" "${STAT_DIR}" "${LOG_DIR}/${sample}.log" \
      "${THREADS_PER_JOB}" "${MM2_PRESET}" "${MIN_MAPQ}"
  done < "${MAP_LIST}"
fi

echo "[$(date)] All done. BAMs in: ${BAM_DIR}"
echo "Stats in: ${STAT_DIR} ; Logs in: ${LOG_DIR}"
