#!/bin/bash

################################################################################
# HBV Variant Calling Pipeline - variants_call_10
# 
# Purpose:
#   Complete pipeline for HBV variant calling from nanopore sequencing data.
#   Includes mapping, filtering, variant calling (iVar + Clair3), and unified
#   reference coordinate transformation.
#
# Input Data:
#   - FASTQ: host_deconv_out_5/{sample_id}_subsampled.trimmed_filtered.porechop/
#   - Consensus: Medaka_consensus_8/r2/{sample_id}/consensus.fasta
#   - View A: consensus_ref_viewa_9/*.csv
#
# Output:
#   - variants_call_10/mapping/
#   - variants_call_10/filtering/
#   - variants_call_10/variants/{ivar,clair3,filtered,transformed,combined,unified}/
#
# Usage:
#   ./11_run_variants_call_10.sh [--step N] [--from N] [--test SAMPLE]
#
# Options:
#   --step N     : Run only step N (1-8)
#   --from N     : Run from step N to end
#   --test SAMPLE: Test mode with single sample
#   --dry-run    : Show what would be run without executing
#   -y           : Skip all confirmations
#
# Version: 1.0
################################################################################

set -e
set -u

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# ========== PATH CONFIGURATION ==========
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Output directory
WORK_DIR="${PROJECT_DIR}/variants_call_10"

# Log directory
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Conda initialization
CONDA_INIT="${HOME}/miniforge3/etc/profile.d/conda.sh"
# ========================================

# Functions
timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

log() {
    echo -e "${GREEN}[$(timestamp)]${NC} $1"
}

log_error() {
    echo -e "${RED}[$(timestamp)] ERROR:${NC} $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[$(timestamp)] WARN:${NC} $1"
}

log_step() {
    echo ""
    echo -e "${PURPLE}========================================${NC}"
    echo -e "${PURPLE}$1${NC}"
    echo -e "${PURPLE}========================================${NC}"
    echo ""
}

start_timer() {
    STEP_START_TIME=$(date +%s)
}

end_timer() {
    local step_end_time=$(date +%s)
    local duration=$((step_end_time - STEP_START_TIME))
    local minutes=$((duration / 60))
    local seconds=$((duration % 60))
    echo -e "${GREEN}Step completed in ${minutes}m ${seconds}s${NC}"
}

# Default parameters
STEP_ONLY=""
START_FROM=1
TEST_SAMPLE=""
DRY_RUN=false
AUTO_YES=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --step)
            STEP_ONLY="$2"
            shift 2
            ;;
        --from)
            START_FROM="$2"
            shift 2
            ;;
        --test)
            TEST_SAMPLE="$2"
            shift 2
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
            head -50 "$0" | grep -E "^#" | sed 's/^# //'
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Title
echo -e "${CYAN}"
echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║           HBV Variant Calling Pipeline v1.0                   ║"
echo "║                   variants_call_10                            ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

log "Working directory: ${WORK_DIR}"
log "Script directory: ${SCRIPT_DIR}"
log "Log directory: ${LOG_DIR}"

if [[ -n "$TEST_SAMPLE" ]]; then
    log "TEST MODE: Single sample = ${TEST_SAMPLE}"
fi

if [[ "$DRY_RUN" == true ]]; then
    log_warn "DRY RUN: Commands will be shown but not executed"
fi

# Create output directory
mkdir -p "${WORK_DIR}"

# Record start time
OVERALL_START=$(date +%s)

# Initialize conda
source "${CONDA_INIT}"

################################################################################
# STEP 1: Mapping
################################################################################
run_step1() {
    log_step "STEP 1: Mapping (minimap2)"
    start_timer

    log "Activating base environment..."
    conda activate base

    local cmd="python3 ${SCRIPT_DIR}/11.1_mapping_batch_individual.py -t 2 -p 4"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step1_mapping.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 1 completed successfully"
        else
            log_error "Step 1 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 2: Filtering
################################################################################
run_step2() {
    log_step "STEP 2: BAM Filtering"
    start_timer

    log "Activating base environment..."
    conda activate base

    local cmd="python3 ${SCRIPT_DIR}/11.2combined_bam_filter_batch.py --all -j 4"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step2_filtering.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 2 completed successfully"
        else
            log_error "Step 2 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 3a: iVar Variant Calling
################################################################################
run_step3a() {
    log_step "STEP 3a: iVar Variant Calling"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local yes_flag=""
    [[ "$AUTO_YES" == true ]] && yes_flag="-y"
    
    local cmd="python3 ${SCRIPT_DIR}/11.3a_variant_call_ivar_batch.py --all -j 4 ${yes_flag}"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step3a_ivar.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 3a completed successfully"
        else
            log_error "Step 3a failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 3b: Clair3 Variant Calling
################################################################################
run_step3b() {
    log_step "STEP 3b: Clair3 Variant Calling"
    start_timer

    log "Activating Clair3 environment..."
    conda activate Clair3

    local cmd="bash ${SCRIPT_DIR}/11.3b_batch_clair3_variants.sh --all"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step3b_clair3.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 3b completed successfully"
        else
            log_error "Step 3b failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 4: Filter & Tier Variants
################################################################################
run_step4() {
    log_step "STEP 4: Filter & Tier Variants (iVar + Clair3)"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local yes_flag=""
    [[ "$AUTO_YES" == true ]] && yes_flag="-y"
    
    local cmd="python3 ${SCRIPT_DIR}/11.4_ivar_clair3_filter_batch.py --all --profile lenient ${yes_flag}"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step4_filter.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 4 completed successfully"
        else
            log_error "Step 4 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 5: Summarize Mutations
################################################################################
run_step5() {
    log_step "STEP 5: Summarize Mutations"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local yes_flag=""
    [[ "$AUTO_YES" == true ]] && yes_flag="-y"
    
    local cmd="python3 ${SCRIPT_DIR}/11.5_summarize_mutations_batch.py --all ${yes_flag}"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step5_summarize.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 5 completed successfully"
        else
            log_error "Step 5 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 6: Coordinate Transformation
################################################################################
run_step6() {
    log_step "STEP 6: Coordinate Transformation"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local cmd="python3 ${SCRIPT_DIR}/11.6_coordinate_transform_batch.py --jobs 4"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step6_transform.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 6 completed successfully"
        else
            log_error "Step 6 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 7: Merge Variants
################################################################################
run_step7() {
    log_step "STEP 7: Merge Variants"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local cmd="python3 ${SCRIPT_DIR}/11.7_merge_variants_batch.py --jobs 4"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step7_merge.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 7 completed successfully"
        else
            log_error "Step 7 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# STEP 8: Flip to Unified Reference
################################################################################
run_step8() {
    log_step "STEP 8: Flip to Unified Reference"
    start_timer

    log "Activating ivar environment..."
    conda activate ivar

    local yes_flag=""
    [[ "$AUTO_YES" == true ]] && yes_flag="-y"
    
    local cmd="python3 ${SCRIPT_DIR}/11.8_flip_to_unified_batch.py --all ${yes_flag}"
    log "Command: $cmd"
    
    if [[ "$DRY_RUN" == false ]]; then
        $cmd 2>&1 | tee "${LOG_DIR}/step8_flip.log"
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log "✓ Step 8 completed successfully"
        else
            log_error "Step 8 failed!"
            exit 1
        fi
    fi
    
    end_timer
}

################################################################################
# Execute steps
################################################################################

if [[ -n "$STEP_ONLY" ]]; then
    # Run single step
    case $STEP_ONLY in
        1) run_step1 ;;
        2) run_step2 ;;
        3a) run_step3a ;;
        3b) run_step3b ;;
        4) run_step4 ;;
        5) run_step5 ;;
        6) run_step6 ;;
        7) run_step7 ;;
        8) run_step8 ;;
        *)
            log_error "Invalid step: $STEP_ONLY. Use 1, 2, 3a, 3b, 4, 5, 6, 7, or 8."
            exit 1
            ;;
    esac
else
    # Run from START_FROM to end
    [[ $START_FROM -le 1 ]] && run_step1
    [[ $START_FROM -le 2 ]] && run_step2
    [[ $START_FROM -le 3 ]] && { run_step3a; run_step3b; }
    [[ $START_FROM -le 4 ]] && run_step4
    [[ $START_FROM -le 5 ]] && run_step5
    [[ $START_FROM -le 6 ]] && run_step6
    [[ $START_FROM -le 7 ]] && run_step7
    [[ $START_FROM -le 8 ]] && run_step8
fi

################################################################################
# Summary
################################################################################

OVERALL_END=$(date +%s)
OVERALL_DURATION=$((OVERALL_END - OVERALL_START))
OVERALL_MINUTES=$((OVERALL_DURATION / 60))
OVERALL_SECONDS=$((OVERALL_DURATION % 60))

log_step "Pipeline Execution Complete!"
echo -e "${GREEN}Total execution time: ${OVERALL_MINUTES}m ${OVERALL_SECONDS}s${NC}"
echo ""

# Quick summary
log_step "Output Summary"
echo "Working directory: ${WORK_DIR}"
echo ""

if [[ -d "${WORK_DIR}/mapping" ]]; then
    echo "Mapping outputs:     $(ls -d ${WORK_DIR}/mapping/*/ 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/filtering" ]]; then
    echo "Filtering outputs:   $(ls -d ${WORK_DIR}/filtering/*/ 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/ivar" ]]; then
    echo "iVar outputs:        $(ls -d ${WORK_DIR}/variants/ivar/*/ 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/clair3" ]]; then
    echo "Clair3 outputs:      $(ls -d ${WORK_DIR}/variants/clair3/*_clair3/ 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/filtered" ]]; then
    echo "Filtered outputs:    $(ls -d ${WORK_DIR}/variants/filtered/*/ 2>/dev/null | grep -v cohort_summary | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/transformed" ]]; then
    echo "Transformed outputs: $(ls -d ${WORK_DIR}/variants/transformed/*/ 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/combined" ]]; then
    echo "Combined outputs:    $(ls ${WORK_DIR}/variants/combined/*_combined.tsv 2>/dev/null | wc -l | tr -d ' ') samples"
fi
if [[ -d "${WORK_DIR}/variants/unified" ]]; then
    echo "Unified outputs:     $(ls ${WORK_DIR}/variants/unified/*_combined.unified.tsv 2>/dev/null | wc -l | tr -d ' ') samples"
fi

echo ""
log "Pipeline complete!"

