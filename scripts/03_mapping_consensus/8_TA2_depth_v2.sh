#!/usr/bin/env bash
# ============================================================================
# 8_TA2_depth_v2.sh - Coverage Analysis & QC Masks
# ============================================================================
# Features:
#   - per-base coverage (mosdepth)
#   - sliding-window coverage
#   - low-coverage & high-variance windows (BED)
#   - homopolymer masks from reference (BED)
#   - per-sample and combined coverage summary (TSV)
# ============================================================================

set -euo pipefail

########################
# Path Configuration
########################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"

BAMS_DIR="${BAMS_DIR:-${PROJECT_DIR}/TA1_map_6_V2/bam}"
REF="${REF:-${HBV_REF:-/path/to/reference.fasta}}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/TA2_depth_7}"

THREADS="${THREADS:-4}"

# Sliding window: HBV ~3.2kb, window=200bp, step=50bp
WIN="${WIN:-200}"
STEP="${STEP:-50}"

# Coverage thresholds
LOW_FRAC="${LOW_FRAC:-0.20}"   # Low = below 20% of median
LOW_ABS="${LOW_ABS:-100}"      # Absolute minimum threshold

# High variance: MAD multiplier
MAD_MULT="${MAD_MULT:-4}"

# Homopolymer masks: >=4bp and >=6bp runs, with padding
HPOLY_MIN4="${HPOLY_MIN4:-4}"
HPOLY_MIN6="${HPOLY_MIN6:-6}"
HPOLY_PAD="${HPOLY_PAD:-5}"

# mosdepth filter flags (exclude unmapped/secondary/supplementary = 0x904)
MOS_FFLAGS="${MOS_FFLAGS:-2308}"

########################
# Tools check
########################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found"; exit 1; }; }
need mosdepth
need samtools
need python3
need gzip
need awk
need sort

mkdir -p "${OUT_DIR}"/{mosdepth,windows,perbase,flags,homopolymers,summary,logs,tmp}

########################
# Reference index + windows
########################
if [[ ! -s "${REF}.fai" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing reference"
  samtools faidx "${REF}"
fi

WIN_BED="${OUT_DIR}/tmp/ref.win${WIN}_step${STEP}.bed"
if [[ ! -s "${WIN_BED}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating sliding windows: ${WIN_BED}"
  awk -v OFS="\t" -v W="${WIN}" -v S="${STEP}" '
    FNR==NR {len[$1]=$2; chr[NR]=$1; nchr=NR; next}
    END{
      for(i=1;i<=nchr;i++){
        c=chr[i]; L=len[c]+0;
        for (start=0; start<L; start+=S) {
          end = start+W; if (end>L) end=L;
          print c, start, end;
          if (end==L) break;
        }
      }
    }' "${REF}.fai" "${REF}.fai" > "${WIN_BED}"
fi

########################
# Homopolymer masks
########################
HP4="${OUT_DIR}/homopolymers/ref.hpoly_ge${HPOLY_MIN4}.bed"
HP6="${OUT_DIR}/homopolymers/ref.hpoly_ge${HPOLY_MIN6}.bed"
HP6_PAD="${OUT_DIR}/homopolymers/ref.hpoly_ge${HPOLY_MIN6}.pad${HPOLY_PAD}.bed"

if [[ ! -s "${HP4}" || ! -s "${HP6}" || ! -s "${HP6_PAD}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building homopolymer masks"
  python3 - "$REF" "$HPOLY_MIN4" "$HPOLY_MIN6" "$HPOLY_PAD" "$HP4" "$HP6" "$HP6_PAD" <<'PY'
import sys
from itertools import groupby
ref = sys.argv[1]
min4 = int(sys.argv[2])
min6 = int(sys.argv[3])
pad  = int(sys.argv[4])
hp4  = sys.argv[5]
hp6  = sys.argv[6]
hp6p = sys.argv[7]

def write_bed(items, path):
    with open(path, 'w') as out:
        for ctg, start, end in items:
            out.write(f"{ctg}\t{start}\t{end}\n")

def padded(items, L, pad):
    for ctg, s, e in items:
        ps = max(0, s - pad)
        pe = min(L[ctg], e + pad)
        yield (ctg, ps, pe)

seqs = {}
lens = {}
ctg=None; buf=[]
with open(ref) as fh:
    for line in fh:
        if line.startswith('>'):
            if ctg:
                seq = ''.join(buf).replace('\n','').upper()
                seqs[ctg]=seq; lens[ctg]=len(seq)
            ctg = line.strip().split()[0][1:]
            buf=[]
        else:
            buf.append(line.strip())
    if ctg:
        seq = ''.join(buf).replace('\n','').upper()
        seqs[ctg]=seq; lens[ctg]=len(seq)

hp4_list=[]; hp6_list=[]
for ctg, seq in seqs.items():
    i=0; n=len(seq)
    while i<n:
        j=i+1
        while j<n and seq[j]==seq[i]:
            j+=1
        runlen = j-i
        if runlen>=min4:
            hp4_list.append((ctg, i, j))
        if runlen>=min6:
            hp6_list.append((ctg, i, j))
        i=j

write_bed(hp4_list, hp4)
write_bed(hp6_list, hp6)
write_bed(list(padded(hp6_list, lens, pad)), hp6p)
PY
fi

########################
# Helper functions
########################
merge_bed() {
  awk -v OFS="\t" '
    NR==1{c=$1;s=$2;e=$3; next}
    ($1==c && $2<=e){ if($3>e)e=$3; next }
    { print c,s,e; c=$1;s=$2;e=$3 }
    END{ if(NR>0) print c,s,e }'
}

compute_coverage_stats() {
  local sample="$1"
  local regions_gz="$2"
  local low_bed="$3"
  local hv_bed="$4"
  local summary_tsv="$5"
  local low_frac="$6"
  local low_abs="$7"
  local mad_mult="$8"

  local tmp_cov
  tmp_cov="$(mktemp)"
  gzip -cd "${regions_gz}" | awk '{print $4}' | sort -n > "${tmp_cov}"
  local n
  n=$(wc -l < "${tmp_cov}" | tr -d ' ')

  if [[ "$n" -eq 0 ]]; then
    > "${low_bed}"
    > "${hv_bed}"
    echo -e "${sample}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0" > "${summary_tsv}"
    rm -f "${tmp_cov}"
    return
  fi

  # Calculate median
  local median
  if (( n % 2 == 1 )); then
    local k=$(( (n+1)/2 ))
    median=$(awk -v k="$k" 'NR==k{print $1}' "${tmp_cov}")
  else
    local k=$(( n/2 ))
    median=$(awk -v k="$k" 'NR==k{a=$1} NR==k+1{printf "%.6f", (a+$1)/2; exit}' "${tmp_cov}")
  fi

  local mean stdev cv mn mx low_thr mad mad_thr low_cnt hv_cnt low_frac_val hv_frac_val
  mean=$(awk '{s+=$1} END{printf "%.6f", s/NR}' "${tmp_cov}")
  stdev=$(awk -v m="$mean" '{d=($1-m); s+=d*d} END{printf "%.6f", sqrt(s/NR)}' "${tmp_cov}")
  cv=$(awk -v m="$mean" -v sd="$stdev" 'BEGIN{printf "%.6f", (m>0? sd/m : 0)}')
  mn=$(head -n1 "${tmp_cov}")
  mx=$(tail -n1 "${tmp_cov}")
  low_thr=$(awk -v m="$median" -v lf="$low_frac" -v la="$low_abs" 'BEGIN{t=lf*m; if (t<la) t=la; printf "%.6f", t}')
  
  mad=$(awk -v med="$median" '{x=$1-med; if(x<0)x=-x; print x}' "${tmp_cov}" | sort -n \
        | awk -v n="$n" '{a[NR]=$1} END{if(n%2){k=(n+1)/2; print a[k]} else{k=n/2; printf "%.6f", (a[k]+a[k+1])/2}}')
  mad_thr=$(awk -v mm="$mad_mult" -v mad="$mad" 'BEGIN{printf "%.6f", mm*mad}')

  low_cnt=$(awk -v t="$low_thr" '$1 < t{c++} END{print 0+c}' "${tmp_cov}")
  hv_cnt=$(awk -v med="$median" -v mt="$mad_thr" '{x=$1-med; if(x<0)x=-x; if(x>mt)c++} END{print 0+c}' "${tmp_cov}")
  low_frac_val=$(awk -v a="$low_cnt" -v n="$n" 'BEGIN{printf "%.6f", a/n}')
  hv_frac_val=$(awk -v a="$hv_cnt" -v n="$n" 'BEGIN{printf "%.6f", a/n}')

  echo -e "${sample}\t${n}\t${mean}\t${median}\t${stdev}\t${cv}\t${mn}\t${mx}\t${low_thr}\t${mad_thr}\t${low_frac_val}\t${hv_frac_val}" > "${summary_tsv}"

  gzip -cd "${regions_gz}" | awk -v t="$low_thr" '$4 < t {print $1"\t"$2"\t"$3}' > "${low_bed}"
  gzip -cd "${regions_gz}" | awk -v med="$median" -v mt="$mad_thr" '{x=$4-med; if(x<0)x=-x; if(x>mt) print $1"\t"$2"\t"$3}' > "${hv_bed}"

  rm -f "${tmp_cov}"
}

########################
# Main loop
########################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TA2 Coverage Analysis starting..."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] BAMs: ${BAMS_DIR}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reference: ${REF}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output: ${OUT_DIR}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================"

find "${BAMS_DIR}" -type f -name "*.sorted.bam" | sort | while read -r BAM; do
  BN="$(basename "${BAM}")"
  ID="${BN%.sorted.bam}"
  PREFIX_BASE="${OUT_DIR}/mosdepth/${ID}.base"
  PREFIX_WIN="${OUT_DIR}/mosdepth/${ID}.w${WIN}s${STEP}"

  echo ""
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing: ${ID}"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')]   mosdepth per-base"
  mosdepth -t "${THREADS}" -F "${MOS_FFLAGS}" "${PREFIX_BASE}" "${BAM}" \
    2> "${OUT_DIR}/logs/${ID}.mosdepth.base.log"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')]   mosdepth windows"
  mosdepth -t "${THREADS}" -F "${MOS_FFLAGS}" -b "${WIN_BED}" "${PREFIX_WIN}" "${BAM}" \
    2> "${OUT_DIR}/logs/${ID}.mosdepth.win.log"

  REG_BED_GZ="${PREFIX_WIN}.regions.bed.gz"
  if [[ ! -s "${REG_BED_GZ}" ]]; then
    echo "WARN: missing ${REG_BED_GZ}" >&2
    continue
  fi

  LOW_RAW="${OUT_DIR}/flags/${ID}.lowcov.raw.bed"
  HV_RAW="${OUT_DIR}/flags/${ID}.highvar.raw.bed"
  SUM_TSV="${OUT_DIR}/summary/${ID}.coverage_summary.tsv"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')]   computing stats"
  compute_coverage_stats "${ID}" "${REG_BED_GZ}" "${LOW_RAW}" "${HV_RAW}" "${SUM_TSV}" "${LOW_FRAC}" "${LOW_ABS}" "${MAD_MULT}"

  sort -k1,1 -k2,2n "${LOW_RAW}" | merge_bed > "${OUT_DIR}/flags/${ID}.lowcov.merged.bed"
  sort -k1,1 -k2,2n "${HV_RAW}"  | merge_bed > "${OUT_DIR}/flags/${ID}.highvar.merged.bed"

  ln -sf "${PREFIX_BASE}.per-base.bed.gz"  "${OUT_DIR}/perbase/${ID}.per-base.bed.gz" 2>/dev/null || true
  ln -sf "${PREFIX_WIN}.regions.bed.gz"    "${OUT_DIR}/windows/${ID}.w${WIN}s${STEP}.bed.gz" 2>/dev/null || true
done

########################
# Aggregate summaries
########################
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building combined summary..."

ALL_SUMMARY="${OUT_DIR}/summary/ALL.coverage_summary.tsv"
echo -e "sample\twindows\tmean\tmedian\tstdev\tcv\tmin\tmax\tlow_thr\tmad_thr\tlow_frac\thighvar_frac" > "${ALL_SUMMARY}"

shopt -s nullglob
for tsv in "${OUT_DIR}/summary"/*.coverage_summary.tsv; do
  bn="$(basename "${tsv}")"
  [[ "${bn}" == "ALL.coverage_summary.tsv" ]] && continue
  cat "${tsv}" >> "${ALL_SUMMARY}"
done

{
  head -n1 "${ALL_SUMMARY}"
  tail -n +2 "${ALL_SUMMARY}" | sort -t$'\t' -k1,1
} > "${ALL_SUMMARY}.tmp" && mv "${ALL_SUMMARY}.tmp" "${ALL_SUMMARY}"

TOTAL_SAMPLES=$(tail -n +2 "${ALL_SUMMARY}" | wc -l | tr -d ' ')

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Complete! Samples: ${TOTAL_SAMPLES}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Summary: ${ALL_SUMMARY}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================"
