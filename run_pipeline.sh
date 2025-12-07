#!/usr/bin/env bash
# ============================================================================
# HBV Nanopore WGS Pipeline - Main Driver Script
# ============================================================================
# 
# Complete pipeline from raw FASTQ to variant calls.
# Runs all modules sequentially for single sample or batch processing.
#
# Usage:
#   ./run_pipeline.sh CONFIG [OPTIONS]
#
# Examples:
#   ./run_pipeline.sh config/config.yaml --sample 10090
#   ./run_pipeline.sh config/config.yaml --batch config/samplesheet.csv
#   ./run_pipeline.sh config/config.yaml --sample 10090 --from 3
#   ./run_pipeline.sh config/config.yaml --sample 10090 --step 4
#
# ============================================================================

set -Eeuo pipefail

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
CONFIG_PARSER="${SCRIPTS_DIR}/utils/config_parser.py"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# ============================================================================
# Helper Functions
# ============================================================================

usage() {
    cat << EOF
HBV Nanopore WGS Pipeline

Usage:
    $(basename "$0") CONFIG [OPTIONS]

Arguments:
    CONFIG          Path to configuration YAML file (required)

Options:
    --sample ID     Process single sample with given ID
    --batch FILE    Process samples from samplesheet CSV
    --step N        Run only step N (1-11)
    --from N        Run from step N to end (default: 1)
    --to N          Run up to step N
    --list          List available steps and exit
    --dry-run       Show what would be run without executing
    -y, --yes       Skip all confirmations
    -h, --help      Show this help message

Steps:
    1   Dorado adapter trimming
    2   Length and quality filtering
    3   Porechop mid-adapter removal
    4   Multi-tool QC (FastQC, NanoPlot, etc.)
    5   Kraken2 contamination screening
    6   Host decontamination
    7   Mapping to unified reference
    8   Depth analysis
    9   Medaka consensus (R1 + R2)
    10  Consensus comparison
    11  Variant calling pipeline

Examples:
    # Single sample, full pipeline
    $(basename "$0") config/config.yaml --sample 10090

    # Batch processing
    $(basename "$0") config/config.yaml --batch samplesheet.csv

    # Start from step 6 (host decontamination)
    $(basename "$0") config/config.yaml --sample 10090 --from 6

    # Run only QC steps (1-5)
    $(basename "$0") config/config.yaml --sample 10090 --from 1 --to 5

    # Dry run to see what would execute
    $(basename "$0") config/config.yaml --sample 10090 --dry-run

EOF
    exit 0
}

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo -e "${GREEN}${msg}${NC}"
    [[ -n "${LOG_FILE:-}" ]] && echo "$msg" >> "$LOG_FILE" || true
}

log_step() {
    local msg="$1"
    echo ""
    echo -e "${PURPLE}========================================${NC}"
    echo -e "${PURPLE}${msg}${NC}"
    echo -e "${PURPLE}========================================${NC}"
    echo ""
    if [[ -n "${LOG_FILE:-}" ]]; then
        echo "" >> "$LOG_FILE"
        echo "========================================" >> "$LOG_FILE"
        echo "$msg" >> "$LOG_FILE"
        echo "========================================" >> "$LOG_FILE"
    fi
}

log_error() {
    local msg="[ERROR] $1"
    echo -e "${RED}${msg}${NC}" >&2
    [[ -n "${LOG_FILE:-}" ]] && echo "$msg" >> "$LOG_FILE" || true
}

log_warn() {
    local msg="[WARN] $1"
    echo -e "${YELLOW}${msg}${NC}"
    [[ -n "${LOG_FILE:-}" ]] && echo "$msg" >> "$LOG_FILE" || true
}

log_info() {
    local msg="[INFO] $1"
    echo -e "${CYAN}${msg}${NC}"
    [[ -n "${LOG_FILE:-}" ]] && echo "$msg" >> "$LOG_FILE" || true
}

check_config() {
    local config_file="$1"
    if [[ ! -f "$config_file" ]]; then
        log_error "Config file not found: $config_file"
        log_info "Copy the example config: cp config/config.example.yaml config/config.yaml"
        exit 1
    fi
}

# Parse config using simple grep-based YAML parsing
# Note: For complex YAML, install PyYAML and use scripts/utils/config_parser.py
get_config() {
    local key="$1"
    local default="${2:-}"
    local leaf_key="${key##*.}"
    local result
    
    # Simple grep-based YAML parsing (handles basic flat key: value cases)
    result=$(grep -E "^[[:space:]]*${leaf_key}:" "$CONFIG_FILE" 2>/dev/null | head -1 | cut -d: -f2- | sed 's/^ *//; s/ *$//; s/^"//; s/"$//' 2>/dev/null) || result=""
    
    # Return result or default
    if [[ -n "$result" ]]; then
        printf '%s' "$result"
    else
        printf '%s' "$default"
    fi
}

check_conda_env() {
    local env_name="$1"
    if ! conda info --envs 2>/dev/null | grep -q "^${env_name}\s"; then
        log_warn "Conda environment not found: $env_name"
        log_info "Create it with: mamba env create -f environment/env_*.yml"
        return 1
    fi
    return 0
}

activate_env() {
    local env_name="$1"
    log_info "Activating conda environment: $env_name"
    
    # Source conda initialization if needed
    if [[ -f "${HOME}/miniforge3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/miniforge3/etc/profile.d/conda.sh"
    elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/anaconda3/etc/profile.d/conda.sh"
    fi
    
    conda activate "$env_name" 2>/dev/null || true
}

list_steps() {
    cat << EOF
Available Pipeline Steps:

  Step   Module              Description
  ────   ──────              ───────────
  1      01_qc               Dorado adapter trimming
  2      01_qc               Length and quality filtering
  3      01_qc               Porechop mid-adapter removal
  4      01_qc               Multi-tool QC (FastQC, NanoPlot, etc.)
  5      01_qc               Kraken2 contamination screening
  6      02_host_deconv      Host decontamination
  7      03_mapping          Mapping to unified reference
  8      03_mapping          Depth analysis
  9      03_consensus        Medaka consensus (R1 + R2)
  10     03_consensus        Consensus comparison
  11     04_variants         Variant calling pipeline

Environment Requirements:
  Steps 1:       dorado
  Steps 2-5:     hbv_qc
  Steps 6-8:     hbv_base / hbv_kraken2
  Steps 9-10:    hbv_medaka
  Step 11:       hbv_variants + hbv_clair3

EOF
    exit 0
}

setup_logging() {
    local sample_id="${1:-batch}"
    LOG_DIR="${WORK_DIR}/logs"
    mkdir -p "$LOG_DIR"
    LOG_FILE="${LOG_DIR}/pipeline_${sample_id}_$(date '+%Y%m%d_%H%M%S').log"
    log_info "Logging to: $LOG_FILE"
}

# ============================================================================
# Step Functions
# ============================================================================

run_step_1() {
    log_step "Step 1: Dorado Adapter Trimming"
    local step_log="${LOG_DIR}/step01_dorado_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 1: Dorado Adapter Trimming
  Command: bash ${SCRIPTS_DIR}/01_qc/1_trim_adapters_dorado_fastq.sh --run
  Input:   ${RAW_FASTQ_DIR}
  Output:  ${WORK_DIR}/fastq_dorado
  Kit:     ${DORADO_KIT}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running Dorado adapter trimming..."
    SRC_DIR="${RAW_FASTQ_DIR}" \
    OUT_DIR="${WORK_DIR}/fastq_dorado" \
    KIT="${DORADO_KIT}" \
        bash "${SCRIPTS_DIR}/01_qc/1_trim_adapters_dorado_fastq.sh" --run 2>&1 | tee "$step_log"
}

run_step_2() {
    log_step "Step 2: Length and Quality Filtering"
    local step_log="${LOG_DIR}/step02_filter_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 2: Length and Quality Filtering
  Command: python ${SCRIPTS_DIR}/01_qc/2_filter_fastq_gz_2-4k_Q17.py
  Input:   ${WORK_DIR}/fastq_dorado
  Output:  ${WORK_DIR}/fastq_filter_2
  Filter:  Length ${MIN_READ_LENGTH}-${MAX_READ_LENGTH}, Q>${MIN_QUALITY}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running length/quality filtering..."
    python "${SCRIPTS_DIR}/01_qc/2_filter_fastq_gz_2-4k_Q17.py" \
        --input "${WORK_DIR}/fastq_dorado" \
        --output "${WORK_DIR}/fastq_filter_2" \
        --min-len "${MIN_READ_LENGTH}" \
        --max-len "${MAX_READ_LENGTH}" \
        --min-quality "${MIN_QUALITY}" 2>&1 | tee "$step_log"
}

run_step_3() {
    log_step "Step 3: Porechop Mid-Adapter Removal"
    local step_log="${LOG_DIR}/step03_porechop_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 3: Porechop Mid-Adapter Removal
  Command: bash ${SCRIPTS_DIR}/01_qc/3_run_porechop_midcheck.sh
  Input:   ${WORK_DIR}/fastq_filter_2
  Output:  ${WORK_DIR}/fastq_porechop_3
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running Porechop..."
    INPUT_DIR="${WORK_DIR}/fastq_filter_2" \
    OUTPUT_DIR="${WORK_DIR}/fastq_porechop_3" \
    THREADS="${THREADS}" \
        bash "${SCRIPTS_DIR}/01_qc/3_run_porechop_midcheck.sh" 2>&1 | tee "$step_log"
}

run_step_4() {
    log_step "Step 4: Multi-Tool QC"
    local step_log="${LOG_DIR}/step04_multiqc_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 4: Multi-Tool QC (FastQC, NanoPlot, nanoQC, SeqKit)
  Command: bash ${SCRIPTS_DIR}/01_qc/4_ont_multi_qc.sh
  Input:   ${WORK_DIR}/fastq_porechop_3
  Output:  ${WORK_DIR}/multi_tool_qc_4
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running multi-tool QC..."
    IN_DIR="${WORK_DIR}/fastq_porechop_3" \
    OUT_DIR="${WORK_DIR}/multi_tool_qc_4" \
    TOTAL_CORES="${THREADS}" \
        bash "${SCRIPTS_DIR}/01_qc/4_ont_multi_qc.sh" 2>&1 | tee "$step_log"
}

run_step_5() {
    log_step "Step 5: Kraken2 Contamination Screening"
    local step_log="${LOG_DIR}/step05_kraken2_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 5: Kraken2 Contamination Screening
  Command: bash ${SCRIPTS_DIR}/01_qc/5_kraken2_nanopore_qc.sh
  Input:   ${WORK_DIR}/fastq_porechop_3
  Output:  ${WORK_DIR}/multi_tool_qc_4
  Database: ${KRAKEN_DB}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running Kraken2 QC..."
    INDIR="${WORK_DIR}/fastq_porechop_3" \
    OUTDIR="${WORK_DIR}/multi_tool_qc_4" \
    KRAKEN_DB="${KRAKEN_DB}" \
    THREADS="${THREADS}" \
        bash "${SCRIPTS_DIR}/01_qc/5_kraken2_nanopore_qc.sh" 2>&1 | tee "$step_log"
}

run_step_6() {
    log_step "Step 6: Host Decontamination"
    local step_log="${LOG_DIR}/step06_host_deconv_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 6: Host Decontamination
  Command: bash ${SCRIPTS_DIR}/02_host_deconv/6.2_host_deconv_mask_hbv.v2.sh
  Input:   ${WORK_DIR}/fastq_porechop_3
  Output:  ${WORK_DIR}/host_deconv_out_5
  HBV Panel: ${HBV_PANEL}
  Human Ref: ${HUMAN_REF}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running host decontamination..."
    RAW_DIR="${WORK_DIR}/fastq_porechop_3" \
    OUTDIR="${WORK_DIR}/host_deconv_out_5" \
    HBV_PANEL="${HBV_PANEL}" \
    HUMAN_REF_FASTA="${HUMAN_REF}" \
    THREADS="${THREADS}" \
        bash "${SCRIPTS_DIR}/02_host_deconv/6.2_host_deconv_mask_hbv.v2.sh" 2>&1 | tee "$step_log"
}

run_step_7() {
    log_step "Step 7: Mapping to Unified Reference"
    local step_log="${LOG_DIR}/step07_mapping_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 7: Mapping to Unified Reference
  Command: bash ${SCRIPTS_DIR}/03_mapping_consensus/7_TA1_map_V2_unmaskreads.sh
  Input:   ${WORK_DIR}/host_deconv_out_5
  Output:  ${WORK_DIR}/TA1_map_6_V2/bam
  Reference: ${HBV_UNIFIED}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running mapping..."
    FASTQ_DIR="${WORK_DIR}/host_deconv_out_5" \
    BAMS_DIR="${WORK_DIR}/TA1_map_6_V2/bam" \
    REF="${HBV_UNIFIED}" \
    THREADS="${THREADS}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/7_TA1_map_V2_unmaskreads.sh" 2>&1 | tee "$step_log"
}

run_step_8() {
    log_step "Step 8: Depth Analysis"
    local step_log="${LOG_DIR}/step08_depth_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 8: Depth Analysis
  Command: bash ${SCRIPTS_DIR}/03_mapping_consensus/8_TA2_depth_v2.sh
  Input:   ${WORK_DIR}/TA1_map_6_V2/bam
  Output:  ${WORK_DIR}/TA2_depth_7
  Reference: ${HBV_UNIFIED}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running depth analysis..."
    BAMS_DIR="${WORK_DIR}/TA1_map_6_V2/bam" \
    REF="${HBV_UNIFIED}" \
    OUT_DIR="${WORK_DIR}/TA2_depth_7" \
    THREADS="${THREADS}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/8_TA2_depth_v2.sh" 2>&1 | tee "$step_log"
}

run_step_9() {
    log_step "Step 9: Medaka Consensus (R1 + R2)"
    local step_log="${LOG_DIR}/step09_medaka_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 9: Medaka Consensus (R1 + R2)
  Commands:
    9.1: ${SCRIPTS_DIR}/03_mapping_consensus/9.1_run_r1_minimap2_and_qc_v2.sh
    9.2: ${SCRIPTS_DIR}/03_mapping_consensus/9.2_run_medaka_from_bam_v2.sh
    9.3: ${SCRIPTS_DIR}/03_mapping_consensus/9.3_run_medaka_round2_v2.sh
  Input:   ${WORK_DIR}/host_deconv_out_5
  Output:  ${WORK_DIR}/Medaka_consensus_8/r2
  Model:   ${MEDAKA_MODEL}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running Round 1 alignment..."
    FASTQ_ROOT="${WORK_DIR}/host_deconv_out_5" \
    REF_FASTA="${HBV_UNIFIED}" \
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r1" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.1_run_r1_minimap2_and_qc_v2.sh" 2>&1 | tee "$step_log"
    
    log_info "Running Round 1 Medaka consensus..."
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r1" \
    REF_FASTA="${HBV_UNIFIED}" \
    MEDAKA_MODEL="${MEDAKA_MODEL}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.2_run_medaka_from_bam_v2.sh" 2>&1 | tee -a "$step_log"
    
    log_info "Running Round 2 polish..."
    R1_CONS_DIR="${WORK_DIR}/Medaka_consensus_8/r1/medaka_R1_consensus" \
    READS_ROOT="${WORK_DIR}/host_deconv_out_5" \
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r2" \
    MEDAKA_MODEL="${MEDAKA_MODEL}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.3_run_medaka_round2_v2.sh" 2>&1 | tee -a "$step_log"
}

run_step_10() {
    log_step "Step 10: Consensus Comparison"
    local step_log="${LOG_DIR}/step10_comparison_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 10: Consensus Comparison with Coordinate Transform
  Command: python ${SCRIPTS_DIR}/03_mapping_consensus/10_consensus_comparison_with_transform.py
  Input:   ${WORK_DIR}/Medaka_consensus_8/r2
  Output:  ${WORK_DIR}/consensus_ref_viewa_9
  Offset:  ${HBV_OFFSET}
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running consensus comparison..."
    HBV_CONSENSUS_DIR="${WORK_DIR}/Medaka_consensus_8/r2" \
    HBV_REF="${HBV_UNIFIED}" \
    HBV_OUTPUT_DIR="${WORK_DIR}/consensus_ref_viewa_9" \
        python "${SCRIPTS_DIR}/03_mapping_consensus/10_consensus_comparison_with_transform.py" \
        --offset "${HBV_OFFSET}" \
        --no-confirm 2>&1 | tee "$step_log"
}

run_step_11() {
    log_step "Step 11: Variant Calling Pipeline"
    local step_log="${LOG_DIR}/step11_variants_${SAMPLE_ID:-batch}.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        cat << EOF
[DRY-RUN] Step 11: Variant Calling Pipeline (iVar + Clair3)
  Command: bash ${SCRIPTS_DIR}/04_variants/11_run_variants_call_10.sh
  Input:   ${WORK_DIR}/Medaka_consensus_8/r2 (consensus)
           ${WORK_DIR}/host_deconv_out_5 (reads)
  Output:  ${WORK_DIR}/variants_call_10
  Log:     ${step_log}
EOF
        return 0
    fi
    
    log_info "Running variant calling pipeline..."
    cd "${SCRIPTS_DIR}/04_variants"
    
    HBV_WORK_DIR="${WORK_DIR}/variants_call_10" \
    CONSENSUS_DIR="${WORK_DIR}/Medaka_consensus_8/r2" \
    FASTQ_DIR="${WORK_DIR}/host_deconv_out_5" \
        bash 11_run_variants_call_10.sh ${AUTO_YES:+-y} 2>&1 | tee "$step_log"
    
    cd "${SCRIPT_DIR}"
}

# ============================================================================
# Parse Arguments
# ============================================================================

CONFIG_FILE=""
SAMPLE_ID=""
BATCH_FILE=""
STEP_ONLY=""
START_FROM=1
END_AT=11
DRY_RUN=false
AUTO_YES=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --batch)
            BATCH_FILE="$2"
            shift 2
            ;;
        --step)
            STEP_ONLY="$2"
            shift 2
            ;;
        --from)
            START_FROM="$2"
            shift 2
            ;;
        --to)
            END_AT="$2"
            shift 2
            ;;
        --list)
            list_steps
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -y|--yes)
            AUTO_YES=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        -*)
            log_error "Unknown option: $1"
            usage
            ;;
        *)
            if [[ -z "$CONFIG_FILE" ]]; then
                CONFIG_FILE="$1"
            else
                log_error "Unexpected argument: $1"
                usage
            fi
            shift
            ;;
    esac
done

# ============================================================================
# Validate Arguments and Load Config
# ============================================================================

if [[ -z "$CONFIG_FILE" ]]; then
    log_error "Configuration file is required"
    usage
fi

check_config "$CONFIG_FILE"

# Load configuration values using simple grep-based YAML parsing
log_info "Loading configuration from: $CONFIG_FILE"

# Disable strict mode temporarily for config loading (grep may fail on missing keys)
set +e

WORK_DIR=$(get_config "project.work_dir" "")
RAW_FASTQ_DIR=$(get_config "project.raw_fastq_dir" "")
RESULTS_DIR=$(get_config "project.results_dir" "${WORK_DIR}/results")

# Reference paths
HBV_PANEL=$(get_config "reference.hbv_panel" "")
HBV_UNIFIED=$(get_config "reference.hbv_unified" "")
HUMAN_REF=$(get_config "reference.human_ref" "")
HBV_GENOME_LENGTH=$(get_config "reference.hbv_genome_length" "3221")
HBV_OFFSET=$(get_config "reference.hbv_rotation_offset" "1823")

# Database paths
KRAKEN_DB=$(get_config "databases.kraken2_db" "")
CLAIR3_MODELS=$(get_config "databases.clair3_models" "")
CLAIR3_MODEL=$(get_config "databases.clair3_model_name" "r1041_e82_400bps_sup_v500")

# Resources
THREADS=$(get_config "resources.threads" "16")

# QC parameters
DORADO_KIT=$(get_config "qc.dorado_kit" "SQK-NBD114-24")
MIN_READ_LENGTH=$(get_config "qc.min_read_length" "2000")
MAX_READ_LENGTH=$(get_config "qc.max_read_length" "4000")
MIN_QUALITY=$(get_config "qc.min_quality" "17")

# Medaka parameters
MEDAKA_MODEL=$(get_config "medaka.model" "r1041_e82_400bps_sup_v500")

# Re-enable strict mode
set -e

# Validate required paths
if [[ -z "$WORK_DIR" ]]; then
    log_error "project.work_dir not set in config"
    exit 1
fi

# Create directories
mkdir -p "$WORK_DIR"
mkdir -p "${WORK_DIR}/logs"

# Setup logging
setup_logging "${SAMPLE_ID:-batch}"

# ============================================================================
# Main Execution
# ============================================================================

log "============================================================"
log "HBV Nanopore WGS Pipeline v1.0"
log "============================================================"
log "Configuration: $CONFIG_FILE"
log "Working directory: $WORK_DIR"
log "Threads: $THREADS"
[[ -n "$SAMPLE_ID" ]] && log "Sample: $SAMPLE_ID"
[[ -n "$BATCH_FILE" ]] && log "Batch file: $BATCH_FILE"

if [[ "$DRY_RUN" == "true" ]]; then
    log_warn "DRY-RUN MODE - No actual commands will be executed"
    log_info "Commands and paths will be shown for review"
fi

# Determine which steps to run
if [[ -n "$STEP_ONLY" ]]; then
    START_FROM="$STEP_ONLY"
    END_AT="$STEP_ONLY"
fi

log "Running steps ${START_FROM} to ${END_AT}"
log ""

# Execute steps
for step in $(seq "$START_FROM" "$END_AT"); do
    case "$step" in
        1)  run_step_1 ;;
        2)  run_step_2 ;;
        3)  run_step_3 ;;
        4)  run_step_4 ;;
        5)  run_step_5 ;;
        6)  run_step_6 ;;
        7)  run_step_7 ;;
        8)  run_step_8 ;;
        9)  run_step_9 ;;
        10) run_step_10 ;;
        11) run_step_11 ;;
        *)
            log_error "Unknown step: $step"
            exit 1
            ;;
    esac
done

log_step "Pipeline Complete!"
log "Results are in: $WORK_DIR"
log "Final variants: ${WORK_DIR}/variants_call_10/variants/unified/"
log "Log file: $LOG_FILE"
