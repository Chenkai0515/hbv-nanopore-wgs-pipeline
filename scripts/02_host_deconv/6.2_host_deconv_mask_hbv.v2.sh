#!/usr/bin/env bash
# 6.2_host_deconv_mask_hbv.v2.sh
set -Eeuo pipefail

THREADS="${THREADS:-16}"

HBV_PANEL="${HBV_PANEL:-/path/to/ref_panel_B_all_subgt_primer_start.cleaned.fasta}"
RAW_DIR="${RAW_DIR:-/path/to/fastq_porechop_3}"
OUTDIR="${OUTDIR:-/path/to/host_deconv_out_v2}"
HUMAN_REF_FASTA="${HUMAN_REF_FASTA:-/path/to/GRCh38.primary_assembly.fa.gz}"

HBV_MIN_MAPQ="${HBV_MIN_MAPQ:-20}"
HBV_MIN_MATCH_LEN="${HBV_MIN_MATCH_LEN:-100}"
HBV_MERGE_DISTANCE="${HBV_MERGE_DISTANCE:-10}"

HUMAN_MIN_MAPQ="${HUMAN_MIN_MAPQ:-20}"
HUMAN_MIN_MATCH_LEN="${HUMAN_MIN_MATCH_LEN:-100}"
HUMAN_MIN_FRAC="${HUMAN_MIN_FRAC:-0.05}"
HUMAN_MIN_PID="${HUMAN_MIN_PID:-0.80}"
HUMAN_MIN_TOTAL_ALN="${HUMAN_MIN_TOTAL_ALN:-200}"
HUMAN_MIN_TOTAL_FRAC="${HUMAN_MIN_TOTAL_FRAC:-0.05}"

TAIL_NEAR_END="${TAIL_NEAR_END:-80}"
TAIL_MAX_LEN="${TAIL_MAX_LEN:-300}"
TAIL_MAX_FRAC="${TAIL_MAX_FRAC:-0.20}"
TAIL_MERGE_DISTANCE="${TAIL_MERGE_DISTANCE:-10}"
MIN_KEEP_LEN="${MIN_KEEP_LEN:-800}"
DISCARD_INTERNAL="${DISCARD_INTERNAL:-1}"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH" >&2; exit 1; }; }
need_cmd minimap2
need_cmd samtools
need_cmd pigz
need_cmd python3
need_cmd seqkit

if [[ -z "${THREADS}" || "${THREADS}" -le 0 ]]; then
  if command -v nproc >/dev/null 2>&1; then THREADS=$(nproc)
  elif [[ "$OSTYPE" == "darwin"* ]]; then THREADS=$(sysctl -n hw.ncpu); else THREADS=8; fi
fi

mkdir -p "${OUTDIR}"/{ref,logs,tmp}
LOG="${OUTDIR}/logs/run_$(date +%Y%m%d_%H%M%S).log"
echo "=== Host deconvolution v2 ===" | tee -a "$LOG"

HBV_MMI="${OUTDIR}/ref/hbv_panel.mmi"
HUMAN_MMI="${OUTDIR}/ref/human.mmi"
[[ -s "${HBV_MMI}" ]] || { echo "[build] mm2 index HBV" | tee -a "$LOG"; minimap2 -d "${HBV_MMI}" "${HBV_PANEL}"; }
[[ -s "${HUMAN_MMI}" ]] || { echo "[build] mm2 index HUMAN" | tee -a "$LOG"; minimap2 -d "${HUMAN_MMI}" "${HUMAN_REF_FASTA}"; }

SUMMARY="${OUTDIR}/summary.tsv"
echo -e "sample\ttotal_reads\thbv_hit_reads\tmasked_reads\tmasked_bases\thost_reads\tviral_enriched_reads\tviral_unmasked_reads\thost_unmasked_reads\ttrimmed_reads\ttrimmed_bp" > "${SUMMARY}"

shopt -s nullglob
for fq in "${RAW_DIR}"/*.fastq.gz; do
  sample=$(basename "$fq"); sample=${sample%.fastq.gz}
  echo ">>> Processing ${sample}" | tee -a "$LOG"
  SAMPDIR="${OUTDIR}/${sample}"; mkdir -p "${SAMPDIR}"
  total_reads=$( (gzip -cd "$fq" | awk 'NR%4==1{c++} END{print c+0}') )

  PAF_HBV="${SAMPDIR}/${sample}.vsHBV.paf.gz"
  if [[ ! -s "${PAF_HBV}" ]]; then
    minimap2 -x map-ont -t "${THREADS}" --secondary=no -c "${HBV_MMI}" "$fq" | pigz -c > "${PAF_HBV}"
  fi

  MASKED_FQ="${SAMPDIR}/${sample}.maskHBV.fastq.gz"
  MASK_LOG="${SAMPDIR}/${sample}.maskHBV.log"
  python3 "$(dirname "$0")/6.1_mask_hbv_regions.py"     --paf "${PAF_HBV}"     --fastq "$fq"     --out "${MASKED_FQ}"     --min-mapq "${HBV_MIN_MAPQ}"     --merge-distance "${HBV_MERGE_DISTANCE}"     --min-match-len "${HBV_MIN_MATCH_LEN}"     2>&1 | tee -a "${MASK_LOG}" | tee -a "$LOG" >/dev/null

  [[ -s "${MASKED_FQ}" ]] || { echo "[${sample}] ERROR: masked FASTQ missing." | tee -a "$LOG"; continue; }
  masked_reads=$(grep -a "\[mask_hbv_regions\]" "${MASK_LOG}" | tail -n1 | sed -E 's/.*masked_reads=([0-9]+).*/\1/')
  masked_bases=$(grep -a "\[mask_hbv_regions\]" "${MASK_LOG}" | tail -n1 | sed -E 's/.*masked_bases=([0-9]+).*/\1/')

  PAF_MASKED_HUMAN="${SAMPDIR}/${sample}.maskHBV.vsHuman.paf.gz"
  minimap2 -x map-ont -t "${THREADS}" --secondary=no -c "${HUMAN_MMI}" "${MASKED_FQ}" | pigz -c > "${PAF_MASKED_HUMAN}"

  HOST_IDS="${SAMPDIR}/host.ids.txt"
  python3 "$(dirname "$0")/6.4_paf_partition_ids.py"     --paf "${PAF_MASKED_HUMAN}"     --out-host-ids "${HOST_IDS}"     --min-mapq "${HUMAN_MIN_MAPQ}"     --min-match-len "${HUMAN_MIN_MATCH_LEN}"     --min-frac "${HUMAN_MIN_FRAC}"     --min-pid "${HUMAN_MIN_PID}"     --min-total-aln "${HUMAN_MIN_TOTAL_ALN}"     --min-total-frac "${HUMAN_MIN_TOTAL_FRAC}"     2>&1 | tee -a "$LOG" >/dev/null

  HOST_FQ="${SAMPDIR}/${sample}.host_only.fastq.gz"
  VIRAL_FQ="${SAMPDIR}/${sample}.viral_enriched.fastq.gz"
  if [[ -s "${HOST_IDS}" ]]; then
    seqkit grep -n -f "${HOST_IDS}" "${MASKED_FQ}" | pigz -c > "${HOST_FQ}"
    seqkit grep -n -v -f "${HOST_IDS}" "${MASKED_FQ}" | pigz -c > "${VIRAL_FQ}"
    host_reads=$(seqkit stats -T "${HOST_FQ}" | awk 'NR==2{print $4}')
    viral_reads=$(seqkit stats -T "${VIRAL_FQ}" | awk 'NR==2{print $4}')
  else
    pigz -c "${MASKED_FQ}" > "${VIRAL_FQ}"
    : > "${HOST_FQ}"
    host_reads=0
    viral_reads="${total_reads}"
  fi

  PAF_UNMASKED_HUMAN="${SAMPDIR}/${sample}.orig.vsHuman.paf.gz"
  if [[ ! -s "${PAF_UNMASKED_HUMAN}" ]]; then
    minimap2 -x map-ont -t "${THREADS}" --secondary=no -c "${HUMAN_MMI}" "$fq" | pigz -c > "${PAF_UNMASKED_HUMAN}"
  fi

  VIRAL_UNMASKED="${SAMPDIR}/${sample}.viral_enriched.unmasked.fastq.gz"
  HOST_UNMASKED="${SAMPDIR}/${sample}.host_only.unmasked.fastq.gz"
  TAIL_LOG="${SAMPDIR}/${sample}.trim_tails.log"
  python3 "$(dirname "$0")/6.3_trim_host_tails.py"     --paf "${PAF_UNMASKED_HUMAN}"     --fastq "$fq"     --out-viral "${VIRAL_UNMASKED}"     --out-host "${HOST_UNMASKED}"     --min-mapq "${HUMAN_MIN_MAPQ}"     --min-match-len "${HUMAN_MIN_MATCH_LEN}"     --min-frac "${HUMAN_MIN_FRAC}"     --min-pid "${HUMAN_MIN_PID}"     --tail-near-end "${TAIL_NEAR_END}"     --tail-max-len "${TAIL_MAX_LEN}"     --tail-max-frac "${TAIL_MAX_FRAC}"     --merge-distance "${TAIL_MERGE_DISTANCE}"     --min-keep-len "${MIN_KEEP_LEN}"     --discard-internal "${DISCARD_INTERNAL}"     2> "${TAIL_LOG}"

  viral_unmasked_reads=$(seqkit stats -T "${VIRAL_UNMASKED}" | awk 'NR==2{print $4}')
  host_unmasked_reads=$(seqkit stats -T "${HOST_UNMASKED}" | awk 'NR==2{print $4}')
  trimmed_reads=$(grep -a 'trimmed_reads=' "${TAIL_LOG}" | sed -E 's/.*trimmed_reads=([0-9]+).*/\1/')
  trimmed_bp=$(grep -a 'trimmed_bp=' "${TAIL_LOG}" | sed -E 's/.*trimmed_bp=([0-9]+).*/\1/')

  hbv_hit_reads=$(gzip -cd "${PAF_HBV}" | awk '{print $1}' | LC_ALL=C sort -u | wc -l | awk '{print $1}')

  echo -e "${sample}\t${total_reads}\t${hbv_hit_reads}\t${masked_reads:-NA}\t${masked_bases:-NA}\t${host_reads}\t${viral_reads}\t${viral_unmasked_reads}\t${host_unmasked_reads}\t${trimmed_reads:-0}\t${trimmed_bp:-0}" >> "${SUMMARY}"
  echo "[${sample}] DONE." | tee -a "$LOG"
done

echo "All samples finished. Summary at: ${SUMMARY}"