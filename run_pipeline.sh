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

EOF
    exit 0
}

log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

log_step() {
    echo ""
    echo -e "${PURPLE}========================================${NC}"
    echo -e "${PURPLE}$1${NC}"
    echo -e "${PURPLE}========================================${NC}"
    echo ""
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_info() {
    echo -e "${CYAN}[INFO]${NC} $1"
}

check_config() {
    local config_file="$1"
    if [[ ! -f "$config_file" ]]; then
        log_error "Config file not found: $config_file"
        log_info "Copy the example config: cp config/config.example.yaml config/config.yaml"
        exit 1
    fi
}

check_conda_env() {
    local env_name="$1"
    if ! conda info --envs | grep -q "^${env_name}\s"; then
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
    
    conda activate "$env_name"
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
  Steps 6-8:     hbv_base
  Steps 9-10:    hbv_medaka
  Step 11:       hbv_variants + hbv_clair3

EOF
    exit 0
}

# ============================================================================
# Step Functions
# ============================================================================

run_step_1() {
    log_step "Step 1: Dorado Adapter Trimming"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/01_qc/1_trim_adapters_dorado_fastq.sh --run"
        return 0
    fi
    
    # Note: This step requires the dorado conda environment
    log_info "Running Dorado adapter trimming..."
    bash "${SCRIPTS_DIR}/01_qc/1_trim_adapters_dorado_fastq.sh" --run \
        --input "${RAW_FASTQ_DIR}" \
        --output "${WORK_DIR}/fastq_dorado" \
        --kit "${DORADO_KIT:-SQK-NBD114-24}"
}

run_step_2() {
    log_step "Step 2: Length and Quality Filtering"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: python ${SCRIPTS_DIR}/01_qc/2_filter_fastq_gz_2-4k_Q17.py"
        return 0
    fi
    
    log_info "Running length/quality filtering..."
    python "${SCRIPTS_DIR}/01_qc/2_filter_fastq_gz_2-4k_Q17.py" \
        --input "${WORK_DIR}/fastq_dorado" \
        --output "${WORK_DIR}/fastq_filter_2"
}

run_step_3() {
    log_step "Step 3: Porechop Mid-Adapter Removal"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/01_qc/3_run_porechop_midcheck.sh"
        return 0
    fi
    
    log_info "Running Porechop..."
    INPUT_DIR="${WORK_DIR}/fastq_filter_2" \
    OUTPUT_DIR="${WORK_DIR}/fastq_porechop_3" \
        bash "${SCRIPTS_DIR}/01_qc/3_run_porechop_midcheck.sh"
}

run_step_4() {
    log_step "Step 4: Multi-Tool QC"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/01_qc/4_ont_multi_qc.sh"
        return 0
    fi
    
    log_info "Running multi-tool QC..."
    IN_DIR="${WORK_DIR}/fastq_porechop_3" \
    OUT_DIR="${WORK_DIR}/multi_tool_qc_4" \
        bash "${SCRIPTS_DIR}/01_qc/4_ont_multi_qc.sh"
}

run_step_5() {
    log_step "Step 5: Kraken2 Contamination Screening"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/01_qc/5_kraken2_nanopore_qc.sh"
        return 0
    fi
    
    log_info "Running Kraken2 QC..."
    INDIR="${WORK_DIR}/fastq_porechop_3" \
    OUTDIR="${WORK_DIR}/multi_tool_qc_4" \
        bash "${SCRIPTS_DIR}/01_qc/5_kraken2_nanopore_qc.sh"
}

run_step_6() {
    log_step "Step 6: Host Decontamination"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/02_host_deconv/6.2_host_deconv_mask_hbv.v2.sh"
        return 0
    fi
    
    log_info "Running host decontamination..."
    RAW_DIR="${WORK_DIR}/fastq_porechop_3" \
    OUTDIR="${WORK_DIR}/host_deconv_out_5" \
    HBV_PANEL="${HBV_PANEL}" \
    HUMAN_REF_FASTA="${HUMAN_REF}" \
    THREADS="${THREADS:-16}" \
        bash "${SCRIPTS_DIR}/02_host_deconv/6.2_host_deconv_mask_hbv.v2.sh"
}

run_step_7() {
    log_step "Step 7: Mapping to Unified Reference"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/03_mapping_consensus/7_TA1_map_V2_unmaskreads.sh"
        return 0
    fi
    
    log_info "Running mapping..."
    FASTQ_DIR="${WORK_DIR}/host_deconv_out_5" \
    BAMS_DIR="${WORK_DIR}/TA1_map_6_V2/bam" \
    REF="${HBV_UNIFIED}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/7_TA1_map_V2_unmaskreads.sh"
}

run_step_8() {
    log_step "Step 8: Depth Analysis"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/03_mapping_consensus/8_TA2_depth_v2.sh"
        return 0
    fi
    
    log_info "Running depth analysis..."
    BAMS_DIR="${WORK_DIR}/TA1_map_6_V2/bam" \
    REF="${HBV_UNIFIED}" \
    OUT_DIR="${WORK_DIR}/TA2_depth_7" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/8_TA2_depth_v2.sh"
}

run_step_9() {
    log_step "Step 9: Medaka Consensus (R1 + R2)"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run medaka consensus pipeline (9.1, 9.2, 9.3)"
        return 0
    fi
    
    log_info "Running Round 1 alignment..."
    FASTQ_ROOT="${WORK_DIR}/host_deconv_out_5" \
    REF_FASTA="${HBV_UNIFIED}" \
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r1" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.1_run_r1_minimap2_and_qc_v2.sh"
    
    log_info "Running Round 1 Medaka consensus..."
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r1" \
    REF_FASTA="${HBV_UNIFIED}" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.2_run_medaka_from_bam_v2.sh"
    
    log_info "Running Round 2 polish..."
    R1_CONS_DIR="${WORK_DIR}/Medaka_consensus_8/r1/medaka_R1_consensus" \
    READS_ROOT="${WORK_DIR}/host_deconv_out_5" \
    OUT_ROOT="${WORK_DIR}/Medaka_consensus_8/r2" \
        bash "${SCRIPTS_DIR}/03_mapping_consensus/9.3_run_medaka_round2_v2.sh"
}

run_step_10() {
    log_step "Step 10: Consensus Comparison"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: python ${SCRIPTS_DIR}/03_mapping_consensus/10_consensus_comparison_with_transform.py"
        return 0
    fi
    
    log_info "Running consensus comparison..."
    HBV_CONSENSUS_DIR="${WORK_DIR}/Medaka_consensus_8/r2" \
    HBV_REF="${HBV_UNIFIED}" \
    HBV_OUTPUT_DIR="${WORK_DIR}/consensus_ref_viewa_9" \
        python "${SCRIPTS_DIR}/03_mapping_consensus/10_consensus_comparison_with_transform.py" \
        --offset 1823 \
        --no-confirm
}

run_step_11() {
    log_step "Step 11: Variant Calling Pipeline"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] Would run: bash ${SCRIPTS_DIR}/04_variants/11_run_variants_call_10.sh"
        return 0
    fi
    
    log_info "Running variant calling pipeline..."
    cd "${SCRIPTS_DIR}/04_variants"
    
    if [[ "$AUTO_YES" == "true" ]]; then
        bash 11_run_variants_call_10.sh -y
    else
        bash 11_run_variants_call_10.sh
    fi
    
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
# Validate Arguments
# ============================================================================

if [[ -z "$CONFIG_FILE" ]]; then
    log_error "Configuration file is required"
    usage
fi

check_config "$CONFIG_FILE"

# Simple YAML parsing (for key paths)
# In production, consider using yq or a proper YAML parser
parse_yaml_value() {
    local key="$1"
    grep -E "^\s*${key}:" "$CONFIG_FILE" | head -1 | sed 's/.*:\s*//' | tr -d '"' | tr -d "'"
}

# Load configuration values
WORK_DIR=$(parse_yaml_value "work_dir")
RAW_FASTQ_DIR=$(parse_yaml_value "raw_fastq_dir")
HBV_PANEL=$(parse_yaml_value "hbv_panel")
HBV_UNIFIED=$(parse_yaml_value "hbv_unified")
HUMAN_REF=$(parse_yaml_value "human_ref")
THREADS=$(parse_yaml_value "threads")
DORADO_KIT=$(parse_yaml_value "dorado_kit")

# Validate required paths
if [[ -z "$WORK_DIR" ]]; then
    log_error "work_dir not set in config"
    exit 1
fi

# Create work directory if needed
mkdir -p "$WORK_DIR"

# ============================================================================
# Main Execution
# ============================================================================

log "HBV Nanopore WGS Pipeline"
log "Configuration: $CONFIG_FILE"
log "Working directory: $WORK_DIR"

if [[ "$DRY_RUN" == "true" ]]; then
    log_warn "DRY-RUN MODE - No actual commands will be executed"
fi

# Determine which steps to run
if [[ -n "$STEP_ONLY" ]]; then
    START_FROM="$STEP_ONLY"
    END_AT="$STEP_ONLY"
fi

log "Running steps ${START_FROM} to ${END_AT}"

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

