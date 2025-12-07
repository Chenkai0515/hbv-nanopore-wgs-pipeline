#!/usr/bin/env bash
set -Eeuo pipefail

# ====== 用户路径（按需改）======
IN_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/fastq_porechop_3"
OUT_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/multi_tool_qc_4"

# ====== 资源与并行策略（已针对 16C/64G 调优，可按需改）======
TOTAL_CORES="${TOTAL_CORES:-16}"     # 总 CPU 核数
JAVA_MEM_GB="${JAVA_MEM_GB:-8}"      # FastQC 的 Java 内存上限
FASTQC_THREADS="${FASTQC_THREADS:-$TOTAL_CORES}"

# NanoPlot：每个任务 2 线程，并发任务数=总核/2（确保总线程≈16）
NANOPLOT_THREADS_PER_JOB="${NANOPLOT_THREADS_PER_JOB:-2}"
NANOPLOT_JOBS=$(( TOTAL_CORES / NANOPLOT_THREADS_PER_JOB ))
if (( NANOPLOT_JOBS < 1 )); then NANOPLOT_JOBS=1; fi

# nanoQC：单线程工具，按总核数并发
NANOQC_JOBS="${NANOQC_JOBS:-$TOTAL_CORES}"

# SeqKit：多线程统计
SEQKIT_THREADS="${SEQKIT_THREADS:-$TOTAL_CORES}"

# nanoQC 最小 read 长度（0 表示不限制；可按需调整，例如 200/500/1000）
NANOQC_MINLEN="${NANOQC_MINLEN:-0}"

# ====== 可执行程序名字（若用自定义路径可在这里改）======
FASTQC_BIN="${FASTQC_BIN:-fastqc}"
NANOPLOT_BIN="${NANOPLOT_BIN:-NanoPlot}"
NANOQC_BIN="${NANOQC_BIN:-nanoQC}"
SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"

# ====== 目录结构 ======
mkdir -p "$OUT_DIR"/{logs,fastqc,nanoplot,nanoqc,seqkit}

# ====== 依赖检查（有哪个缺就退出）======
for exe in "$FASTQC_BIN" "$NANOPLOT_BIN" "$NANOQC_BIN" "$SEQKIT_BIN"; do
  if ! command -v "$exe" >/dev/null 2>&1; then
    echo "[ERROR] 未找到可执行程序：$exe ；请先通过 conda/bioconda 或 pip 安装。" >&2
    exit 1
  fi
done

# ====== 收集输入 FASTQ 列表 ======
FASTQS=()
while IFS= read -r -d '' f; do FASTQS+=("$f"); done < <(
  find "$IN_DIR" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -print0
)
if (( ${#FASTQS[@]} == 0 )); then
  echo "[ERROR] 在 $IN_DIR 下未找到 *.fastq.gz / *.fq.gz 文件" >&2
  exit 1
fi

# ====== 记录版本信息 ======
{
  echo "Started: $(date)"
  "$FASTQC_BIN" --version || true
  "$NANOPLOT_BIN" --version || true
  "$NANOQC_BIN" --version || true
  "$SEQKIT_BIN" version || true
} | tee "$OUT_DIR/versions.txt"

# ====== Step 1: FastQC（多线程并发处理多个文件）======
# FastQC 的 -t/--threads 控制“同时处理的文件数”，每线程约 250MB 内存。
# 因此这里直接把全部文件交给 fastqc，让它用 16 线程分摊。 详见 fastqc man page。
# （如果文件非常多，也可拆成多批次；一般直接一次跑最简单）
echo "[INFO] Running FastQC with $FASTQC_THREADS threads ..."
env _JAVA_OPTIONS="-Xmx${JAVA_MEM_GB}g -Xms512m" \
  "$FASTQC_BIN" -t "$FASTQC_THREADS" -o "$OUT_DIR/fastqc" "${FASTQS[@]}" \
  2>&1 | tee "$OUT_DIR/logs/fastqc.log"

# ====== Step 2: NanoPlot（每样本单独目录；2 线程/样本，NANOPLOT_JOBS 并发）======
echo "[INFO] Running NanoPlot: $NANOPLOT_JOBS jobs in parallel, $NANOPLOT_THREADS_PER_JOB threads/job ..."
printf '%s\0' "${FASTQS[@]}" | xargs -0 -n 1 -P "$NANOPLOT_JOBS" -I {} bash -c '
  set -Eeuo pipefail
  in="$1"
  bn=$(basename "$in")
  sample=${bn%.fastq.gz}; sample=${sample%.fq.gz}
  od="'"$OUT_DIR"'/nanoplot/${sample}"
  log="'"$OUT_DIR"'/logs/nanoplot_${sample}.log"
  mkdir -p "$od"
  "'"$NANOPLOT_BIN"'" -t "'"$NANOPLOT_THREADS_PER_JOB"'" --fastq "$in" \
    --tsv_stats --N50 --loglength --huge \
    -o "$od" \
    > "$log" 2>&1
' _ {}

# ====== Step 3: nanoQC（每样本单独目录；按核数并发）======
echo "[INFO] Running nanoQC: $NANOQC_JOBS jobs in parallel ..."
printf '%s\0' "${FASTQS[@]}" | xargs -0 -n 1 -P "$NANOQC_JOBS" -I {} bash -c '
  set -Eeuo pipefail
  in="$1"
  bn=$(basename "$in")
  sample=${bn%.fastq.gz}; sample=${sample%.fq.gz}
  od="'"$OUT_DIR"'/nanoqc/${sample}"
  log="'"$OUT_DIR"'/logs/nanoqc_${sample}.log"
  mkdir -p "$od"
  if [ "'"$NANOQC_MINLEN"'" -gt 0 ]; then
    "'"$NANOQC_BIN"'" -o "$od" -l "'"$NANOQC_MINLEN"'" "$in" > "$log" 2>&1
  else
    "'"$NANOQC_BIN"'" -o "$od" "$in" > "$log" 2>&1
  fi
' _ {}

# ====== Step 4: SeqKit 统计（一次性 TSV 汇总，含 N50 等）======
echo "[INFO] Running SeqKit stats with $SEQKIT_THREADS threads ..."
"$SEQKIT_BIN" stats -a -T -j "$SEQKIT_THREADS" "${FASTQS[@]}" \
  | tee "$OUT_DIR/seqkit/seqkit_stats.tsv" > /dev/null

# ====== Step 5（可选）：MultiQC 汇总（若已安装则自动生成）======
if command -v multiqc >/dev/null 2>&1; then
  echo "[INFO] Running MultiQC to aggregate reports ..."
  multiqc "$OUT_DIR" -o "$OUT_DIR/multiqc" 2>&1 | tee "$OUT_DIR/logs/multiqc.log"
else
  echo "[INFO] MultiQC 未检测到，跳过该步骤（如需总览报告：conda install -c bioconda multiqc）"
fi

echo "[DONE] 全部完成！输出目录：$OUT_DIR"
