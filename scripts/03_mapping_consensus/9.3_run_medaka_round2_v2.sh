#!/usr/bin/env bash
set -euo pipefail

# ========= 用户需确认的参数 =========
MODEL="r1041_e82_400bps_sup_v5.0.0"  # 预训练模型（与 Dorado v5 sup 匹配）

# 输入：9.2 产生的 R1 共识目录
R1_CONS_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r1/medaka_R1_consensus"

# 输入：原始 reads 目录（与 9.1 相同）
READS_ROOT="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/host_deconv_out_5"

# 输出：R2 共识目录
OUT_ROOT="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r2"

# 并行相关（可按需改）：每个样本使用的线程数；并发样本数（默认=CPU核心数）
THREADS_PER_JOB="${THREADS_PER_JOB:-1}"

if command -v nproc >/dev/null 2>&1; then
  CORES=$(nproc)
else
  CORES=$(sysctl -n hw.ncpu)
fi
JOBS="${JOBS:-$CORES}"

# ========= 预检 =========
command -v medaka_consensus >/dev/null 2>&1 || {
  echo "[ERROR] 未找到 medaka_consensus。请确认 Medaka 2.1.1 已正确安装到当前环境。" >&2
  exit 1
}
mkdir -p "$OUT_ROOT"

# ========= 生成每个样本的运行命令列表 =========
CMDS_FILE="$OUT_ROOT/medaka_r2_commands.txt"
: > "$CMDS_FILE"

# 遍历第一轮共识目录下的所有样本文件夹
while IFS= read -r -d '' SAMPLE_DIR; do
  SAMPLE_ID="$(basename "$SAMPLE_DIR")"
  REF="$SAMPLE_DIR/consensus.fasta"
  
  # 构建 reads 文件路径（适配新的目录结构）
  READS_SUBDIR="${READS_ROOT}/${SAMPLE_ID}_subsampled.trimmed_filtered.porechop"
  READS="${READS_SUBDIR}/${SAMPLE_ID}_subsampled.trimmed_filtered.porechop.viral_enriched.unmasked.fastq.gz"
  
  OUTDIR="$OUT_ROOT/$SAMPLE_ID"

  if [[ ! -s "$REF" ]]; then
    echo "echo '[SKIP] $SAMPLE_ID: 缺少R1共识 $REF' 1>&2" >> "$CMDS_FILE"
    continue
  fi
  if [[ ! -s "$READS" ]]; then
    echo "echo '[SKIP] $SAMPLE_ID: 缺少reads $READS' 1>&2" >> "$CMDS_FILE"
    continue
  fi

  # 组装单样本命令：第二轮 medaka 共识 + 清理BAM
  mkdir -p "$OUTDIR"
  printf "medaka_consensus -i '%s' -d '%s' -o '%s' -m '%s' -t %d && \
find '%s' -maxdepth 1 -type f \\( -name '*.bam' -o -name '*.bai' \\) -delete\n" \
    "$READS" "$REF" "$OUTDIR" "$MODEL" "$THREADS_PER_JOB" "$OUTDIR" >> "$CMDS_FILE"
done < <(find "$R1_CONS_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

# ========= 并行执行 =========
echo "[INFO] 将以并发 $JOBS 个样本、每个样本 $THREADS_PER_JOB 线程 运行 medaka 共识（Round 2）。"
if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" --bar < "$CMDS_FILE"
else
  # 回退到 xargs 并行
  # macOS 的 xargs 支持 -P 并行；每行一个完整命令
  cat "$CMDS_FILE" | xargs -I CMD -P "$JOBS" bash -lc CMD
fi

echo "[DONE] Medaka Round 2 全部完成。输出位于：$OUT_ROOT"

