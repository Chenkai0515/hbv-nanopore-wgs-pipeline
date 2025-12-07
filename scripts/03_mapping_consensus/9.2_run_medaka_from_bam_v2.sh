#!/usr/bin/env bash
set -euo pipefail

############################################
# 可调参数
############################################
# 输入：9.1 产生的 BAM 目录
OUT_ROOT="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r1"
BAM_DIR="${OUT_ROOT}/minimap2_align/bam"

# 统一参考序列（与9.1相同）
REF_FASTA="/Users/jck/Desktop/ref/gold/reordered_1824_3221_reference_A2763.fasta"

# 输出：9.2 专用子目录
MEDAKA_R1_DIR="${OUT_ROOT}/medaka_R1_consensus"

# 模型与线程
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"
THREADS=16

############################################
# 依赖检查与存在性检查
############################################
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' 未安装或不在 PATH"; exit 1; }; }
need_cmd samtools
need_cmd medaka

# Medaka 主版本检查（仅提示）
MEDAKA_VER=$(medaka --version | awk '{print $2}')
echo "[INFO] Detected medaka version: ${MEDAKA_VER}"
echo "[INFO] 本脚本已按 medaka ≥2.0 的 CLI 适配（inference/sequence 使用位置参数）。"

# 参考序列是否存在
if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: 参考序列不存在：${REF_FASTA}"
  exit 1
fi

mkdir -p "${MEDAKA_R1_DIR}"

# 为参考序列建立 faidx
if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "[INFO] samtools faidx 索引参考序列：${REF_FASTA}"
  samtools faidx "${REF_FASTA}"
fi

############################################
# 主循环
############################################
shopt -s nullglob
bam_list=("${BAM_DIR}"/*.bam)
if (( ${#bam_list[@]} == 0 )); then
  echo "ERROR: 未在 ${BAM_DIR} 找到 .bam 文件。"
  exit 1
fi

echo "[INFO] 找到 ${#bam_list[@]} 个 BAM 文件"

for BAM in "${bam_list[@]}"; do
  base="$(basename "${BAM}")"
  # 提取样本ID: 10090.filtered.mapq20.bam -> 10090
  sample_id="${base%%.*}"

  sample_dir="${MEDAKA_R1_DIR}/${sample_id}"
  mkdir -p "${sample_dir}"

  hdf_out="${sample_dir}/${sample_id}.medaka.hdf"
  consensus_out="${sample_dir}/consensus.fasta"
  log_file="${sample_dir}/medaka.log"

  echo "========================================================"
  echo "[INFO] 处理样本 ${sample_id}"
  echo "[INFO] 参考序列：${REF_FASTA}"
  echo "[INFO] BAM：${BAM}"
  echo "[INFO] 模型：${MEDAKA_MODEL} | 线程：${THREADS}"
  echo "========================================================"

  # 1) medaka inference（注意：2.x 为位置参数：<bam> <hdf>）
  if [[ -s "${hdf_out}" ]]; then
    echo "[SKIP] 已存在 HDF：${hdf_out}"
  else
    echo "[RUN ] medaka inference → ${hdf_out}"
    medaka inference \
      "${BAM}" \
      "${hdf_out}" \
      --model "${MEDAKA_MODEL}" \
      --threads "${THREADS}" \
      2>&1 | tee "${log_file}"
  fi

  # 2) medaka sequence（注意：2.x 只吃 HDF → 输出fasta；不再需要 -r）
  if [[ -s "${consensus_out}" ]]; then
    echo "[SKIP] 已存在 consensus：${consensus_out}"
  else
    echo "[RUN ] medaka sequence → ${consensus_out}"
    medaka sequence \
        --threads "${THREADS}" \
        "${hdf_out}" \
        "${REF_FASTA}" \
        "${consensus_out}" \
        2>&1 | tee -a "${log_file}"
  fi

  echo "[DONE] ${sample_id} → ${consensus_out}"
done

echo "全部样本处理完成。输出位于：${MEDAKA_R1_DIR}"

