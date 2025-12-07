#!/usr/bin/env bash
# 7_TA1_map V2 - Updated for unmasked reads
# Align ONT (Nanopore) reads to a reference using minimap2 with secondary on/off.
# Outputs: sorted BAMs + BAI, alignment stats, and soft-clip/secondary hotspots per sample & mode.
# V2: Process unmasked viral-enriched reads from host_deconv_out_5.v2

set -euo pipefail

############################
# User paths (from prompt) #
############################
SAMPLES_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/host_deconv_out_5"
OUT_DIR="/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/TA1_map_6_V2"
REF="/Users/jck/Desktop/ref/gold/reordered_1824_3221_reference_A2763.fasta"

###################################
# Performance / resource settings #
###################################
THREADS_TOTAL="${THREADS_TOTAL:-16}"     # total CPU threads available
MM2_THREADS="${MM2_THREADS:-4}"          # threads per mapping job
SORT_MEM="${SORT_MEM:-1G}"               # samtools sort memory per thread
MAX_JOBS=$(( THREADS_TOTAL / MM2_THREADS ))
(( MAX_JOBS >= 1 )) || MAX_JOBS=1

###################################
# Tools check                     #
###################################
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not found in PATH."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH."; exit 1; }

mkdir -p "${OUT_DIR}"/{logs,hotspots,stats,bam}

###################################
# Build reference index if needed #
###################################
REF_MMI="${REF%.*}.mmi"
if [[ ! -s "${REF_MMI}" ]]; then
  echo "[`date`] Building minimap2 index: ${REF} -> ${REF_MMI}"
  minimap2 -d "${REF_MMI}" "${REF}"
fi

###################################
# Discover FASTQ files            #
###################################
mapfile -t FASTQS < <(find "${SAMPLES_DIR}" -type f -name "*.viral_enriched.unmasked.fastq.gz" | sort)
if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "ERROR: No *.viral_enriched.unmasked.fastq.gz found under ${SAMPLES_DIR}"; exit 1;
fi

echo "[`date`] Found ${#FASTQS[@]} FASTQ files."
echo "[`date`] Processing unmasked viral-enriched reads (V2)"

###################################
# Helpers                         #
###################################
# Extract sample ID (before first dot) from file name like:
#   10090.porechop.viral_enriched.unmasked.fastq.gz  -> 10090
get_sample_id() {
  local fq="$1"
  local fn bn sid
  fn="$(basename "$fq")"
  bn="${fn%.viral_enriched.unmasked.fastq.gz}"    # e.g. 10090.porechop
  sid="${bn%%.*}"                                  # e.g. 10090
  printf "%s" "$sid"
}

# Produce hotspots for a BAM:
# - soft-clip hotspots (primary only; exclude unmapped/secondary/supplementary)
#   5' clip counted at leftmost POS; 3' clip counted at reference end pos
# - secondary start hotspots (flag 0x100) at POS
# Output:
#   hotspots/${prefix}.hotspots_all.tsv        (all positions with counts)
#   hotspots/${prefix}.hotspots_top100.tsv     (header + top 100 by total)
make_hotspots() {
  local bam="$1"
  local prefix="$2"
  local sc_tsv="${OUT_DIR}/hotspots/${prefix}.softclip_counts.tsv"
  local sec_tsv="${OUT_DIR}/hotspots/${prefix}.secondary_counts.tsv"
  local all_tsv="${OUT_DIR}/hotspots/${prefix}.hotspots_all.tsv"
  local top_tsv="${OUT_DIR}/hotspots/${prefix}.hotspots_top100.tsv"

  # Soft-clip counts (primary alignments only: -F 0x904)
  samtools view -@ "${MM2_THREADS}" -F 0x904 "${bam}" \
  | awk '
    BEGIN{OFS="\t";}
    function parse_cigar(c, op_lens, op_codes,   i,num,ch,n) {
      n=0; num="";
      for (i=1; i<=length(c); i++) {
        ch=substr(c,i,1);
        if (ch ~ /[0-9]/) { num = num ch; }
        else { n++; op_lens[n]=num+0; op_codes[n]=ch; num=""; }
      }
      return n;
    }
    function ref_aln_len(c,   op_lens,op_codes,n,i,len) {
      len=0; n=parse_cigar(c,op_lens,op_codes);
      for (i=1;i<=n;i++) if (op_codes[i]=="M" || op_codes[i]=="D" || op_codes[i]=="N" || op_codes[i]=="=" || op_codes[i]=="X") len+=op_lens[i];
      return len;
    }
    function first_op_S_len(c, op_lens,op_codes,n){ n=parse_cigar(c,op_lens,op_codes); return (n>=1 && op_codes[1]=="S")?op_lens[1]:0; }
    function last_op_S_len(c,  op_lens,op_codes,n){ n=parse_cigar(c,op_lens,op_codes); return (n>=1 && op_codes[n]=="S")?op_lens[n]:0; }
    {
      r=$3; p=$4+0; cg=$6;
      s5=first_op_S_len(cg);
      s3=last_op_S_len(cg);
      ral=ref_aln_len(cg);
      pe=p+ral-1;
      if (s5>0) { k=r SUBSEP p; sc5[k]++; }
      if (s3>0) { k=r SUBSEP pe; sc3[k]++; }
    }
    END{
      for (k in sc5){ split(k,a,SUBSEP); print a[1], a[2], "softclip5", sc5[k]; }
      for (k in sc3){ split(k,a,SUBSEP); print a[1], a[2], "softclip3", sc3[k]; }
    }' > "${sc_tsv}"

  # Secondary start hotspots (flag 0x100)
  samtools view -@ "${MM2_THREADS}" -f 0x100 "${bam}" \
  | awk 'BEGIN{OFS="\t"} {k=$3 SUBSEP $4; c[k]++} END{for(k in c){split(k,a,SUBSEP); print a[1],a[2],"secondary_start",c[k];}}' > "${sec_tsv}"

  # Merge and rank
  awk 'BEGIN{OFS="\t"}
       FNR==NR { key=$1 SUBSEP $2; if($3=="softclip5") sc5[key]=$4; else if($3=="softclip3") sc3[key]=$4; all[key]=1; next }
       { key=$1 SUBSEP $2; sec[key]=$4; all[key]=1 }
       END{
         print "contig","pos","softclip5","softclip3","secondary","total";
         for (k in all) {
           s5=(k in sc5)?sc5[k]:0; s3=(k in sc3)?sc3[k]:0; s=(k in sec)?sec[k]:0;
           split(k,a,SUBSEP); print a[1],a[2],s5,s3,s,s5+s3+s;
         }
       }' "${sc_tsv}" "${sec_tsv}" \
  | sort -k6,6nr -k1,1 -k2,2n > "${all_tsv}"

  # header + top100
  { head -n 1 "${all_tsv}"; tail -n +2 "${all_tsv}" | head -n 100; } > "${top_tsv}"
}

# One mapping job: fq + mode (sec_on|sec_off)
run_one() {
  local fq="$1"; local mode="$2"
  local sid; sid="$(get_sample_id "$fq")"
  local sec_flag log prefix bam

  case "$mode" in
    sec_on)  sec_flag="yes" ;;
    sec_off) sec_flag="no"  ;;
    *) echo "Unknown mode: $mode"; return 2 ;;
  esac

  prefix="${sid}.${mode}"
  log="${OUT_DIR}/logs/${prefix}.log"
  bam="${OUT_DIR}/bam/${prefix}.sorted.bam"

  echo "[`date`] ${sid} ${mode} : mapping ..." | tee -a "${log}"

  # Mapping + coordinate sort + index
  # Note: pipe straight into sort to avoid temporary BAM
  set -o pipefail
  minimap2 -t "${MM2_THREADS}" -x map-ont -a -Y --secondary="${sec_flag}" "${REF_MMI}" "${fq}" \
    2>> "${log}" \
  | samtools sort -@ "${MM2_THREADS}" -m "${SORT_MEM}" -o "${bam}" - \
    2>> "${log}"

  samtools index -@ "${MM2_THREADS}" "${bam}" 2>> "${log}"

  # Stats
  samtools flagstat -@ "${MM2_THREADS}" "${bam}" > "${OUT_DIR}/stats/${prefix}.flagstat.txt"
  samtools stats    -@ "${MM2_THREADS}" "${bam}" > "${OUT_DIR}/stats/${prefix}.stats.txt"
  samtools idxstats "${bam}" > "${OUT_DIR}/stats/${prefix}.idxstats.tsv"

  # Hotspots
  make_hotspots "${bam}" "${prefix}"

  echo "[`date`] ${sid} ${mode} : done." | tee -a "${log}"
}

###################################
# Submit jobs with concurrency    #
###################################
echo "[`date`] Using THREADS_TOTAL=${THREADS_TOTAL}, MM2_THREADS=${MM2_THREADS}, MAX_JOBS=${MAX_JOBS}"

job_guard() {
  # Limit background jobs to MAX_JOBS
  while (( $(jobs -pr | wc -l | tr -d ' ') >= MAX_JOBS )); do sleep 0.5; done
}

for fq in "${FASTQS[@]}"; do
  for mode in sec_on sec_off; do
    job_guard
    run_one "${fq}" "${mode}" &
  done
done

wait
echo "[`date`] All jobs completed. Outputs under: ${OUT_DIR}"
