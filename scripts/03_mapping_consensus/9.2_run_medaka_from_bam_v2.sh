#!/usr/bin/env bash
# Round 1 Medaka consensus from BAM files
# Uses medaka inference + sequence workflow
set -euo pipefail

############################################
# Configuration
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

OUT_ROOT="${OUT_ROOT:-${PROJECT_DIR}/Medaka_consensus_8/r1}"
BAM_DIR="${BAM_DIR:-${OUT_ROOT}/minimap2_align/bam}"
REF_FASTA="${REF_FASTA:-${HBV_REF:-/path/to/reference.fasta}}"
MEDAKA_R1_DIR="${MEDAKA_R1_DIR:-${OUT_ROOT}/medaka_R1_consensus}"

MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v5.0.0}"
THREADS="${THREADS:-16}"

############################################
# Dependency check
############################################
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need_cmd samtools
need_cmd medaka

MEDAKA_VER=$(medaka --version | awk '{print $2}')
echo "[INFO] Detected medaka version: ${MEDAKA_VER}"
echo "[INFO] Using medaka >= 2.0 CLI (inference/sequence positional args)"

if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Reference not found: ${REF_FASTA}"
  exit 1
fi

mkdir -p "${MEDAKA_R1_DIR}"

if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "[INFO] Indexing reference"
  samtools faidx "${REF_FASTA}"
fi

############################################
# Main loop
############################################
shopt -s nullglob
bam_list=("${BAM_DIR}"/*.bam)
if (( ${#bam_list[@]} == 0 )); then
  echo "ERROR: No BAM files in ${BAM_DIR}"
  exit 1
fi

echo "[INFO] Found ${#bam_list[@]} BAM files"

for BAM in "${bam_list[@]}"; do
  base="$(basename "${BAM}")"
  sample_id="${base%%.*}"

  sample_dir="${MEDAKA_R1_DIR}/${sample_id}"
  mkdir -p "${sample_dir}"

  hdf_out="${sample_dir}/${sample_id}.medaka.hdf"
  consensus_out="${sample_dir}/consensus.fasta"
  log_file="${sample_dir}/medaka.log"

  echo "========================================================"
  echo "[INFO] Processing: ${sample_id}"
  echo "[INFO] Model: ${MEDAKA_MODEL} | Threads: ${THREADS}"
  echo "========================================================"

  # Step 1: medaka inference
  if [[ -s "${hdf_out}" ]]; then
    echo "[SKIP] HDF exists: ${hdf_out}"
  else
    echo "[RUN ] medaka inference -> ${hdf_out}"
    medaka inference \
      "${BAM}" \
      "${hdf_out}" \
      --model "${MEDAKA_MODEL}" \
      --threads "${THREADS}" \
      2>&1 | tee "${log_file}"
  fi

  # Step 2: medaka sequence
  if [[ -s "${consensus_out}" ]]; then
    echo "[SKIP] Consensus exists: ${consensus_out}"
  else
    echo "[RUN ] medaka sequence -> ${consensus_out}"
    medaka sequence \
        --threads "${THREADS}" \
        "${hdf_out}" \
        "${REF_FASTA}" \
        "${consensus_out}" \
        2>&1 | tee -a "${log_file}"
  fi

  echo "[DONE] ${sample_id} -> ${consensus_out}"
done

echo "All samples complete. Output: ${MEDAKA_R1_DIR}"
