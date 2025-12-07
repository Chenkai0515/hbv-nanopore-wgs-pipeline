#!/usr/bin/env bash
set -euo pipefail

# ========= HBV Sample Batch Variant Detection Script =========
# Function: Perform Clair3 variant detection on all HBV samples, generating detailed reports
# Support: Single sample testing, batch processing, 16-core parallel, detailed reporting
# New in v1.1: Haploid sensitive calling mode, sequence head/tail variant detection, optimized quality thresholds
# Version: v2.0 (Adapted for variants_call_10 pipeline)
#
# Required Environment:
#   - Clair3 conda environment with all dependencies
#     Create via:
#       conda create -n Clair3 -c bioconda clair3 samtools bcftools parallel
#       conda activate Clair3
#
# Required Model:
#   - Download from GitHub: https://github.com/HKU-BAL/Clair3
#   - Model: r1041_e82_400bps_sup_v500
#   - Installation:
#       cd /path/to/conda/envs/Clair3/bin/
#       mkdir -p models
#       cd models
#       wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v500.tar.gz
#       tar -xzf r1041_e82_400bps_sup_v500.tar.gz
#       rm r1041_e82_400bps_sup_v500.tar.gz

# ========= PATH CONFIGURATION =========
# Get script directory for relative paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Working directory (can be set via HBV_WORK_DIR environment variable)
WORK_DIR="${HBV_WORK_DIR:-${PROJECT_DIR}/variants_call_10}"

# Consensus sequences directory (can be set via CONSENSUS_DIR environment variable)
CONSENSUS_DIR="${CONSENSUS_DIR:-${PROJECT_DIR}/Medaka_consensus_8/r2}"

# Filtering directory (BAM files)
FILTERING_DIR="${WORK_DIR}/filtering"

# Output directories
OUTPUT_DIR="${WORK_DIR}/variants/clair3"
REPORT_DIR="${OUTPUT_DIR}/reports"
# ========================================

# Clair3 Configuration (set CLAIR3_PATH if Clair3 is not in PATH)
CLAIR3_PATH="${CLAIR3_PATH:-}"
# Note: Model path should be set according to your conda environment
# Typical path: /path/to/conda/envs/Clair3/bin/models/r1041_e82_400bps_sup_v500
# Please update MODEL_DIR to match your installation
MODEL_DIR="${CONDA_PREFIX}/bin/models/r1041_e82_400bps_sup_v500"
THREADS_PER_SAMPLE=4  # Number of threads per sample
MAX_PARALLEL_SAMPLES=4  # Maximum parallel samples (4 samples x 4 threads = 16 cores)

# Detection Parameters
SNP_MIN_AF=0.01
INDEL_MIN_AF=0.05
MIN_MQ=10
MIN_COV=200

# New experimental parameters (v1.1)
# --var_pct_full=1.0: Include all 0/1 and 1/1 variants in full-alignment mode calling (default 0.3)
# --ref_pct_full=1.0: Include all 0/0 variants in full-alignment mode calling (default 0.1 for ont)
# --haploid_sensitive: Enable haploid sensitive mode, both 0/1 and 1/1 considered as variants
# --enable_variant_calling_at_sequence_head_and_tail: Enable variant detection in 16bp window at sequence head/tail

# ========= Color Output Definitions =========
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# ========= Helper Functions =========
log_info() {
    echo -e "${GREEN}[INFO $(date +'%H:%M:%S')]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN $(date +'%H:%M:%S')]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR $(date +'%H:%M:%S')]${NC} $1" >&2
}

log_progress() {
    echo -e "${BLUE}[PROGRESS $(date +'%H:%M:%S')]${NC} $1"
}

# Display help information
show_help() {
    cat << EOF
HBV Sample Batch Variant Detection Script v2.0 (variants_call_10)

Usage:
    $0 [options] [sample_IDs...]

Options:
    -h, --help          Display this help information
    -t, --test SAMPLE   Single sample test mode (e.g.: -t 10090)
    -a, --all           Batch process all samples
    -l, --list          List all available samples
    -r, --report        Generate report only (do not run detection)
    -c, --continue      Continue interrupted batch processing
    --threads N         Threads per sample (default: 4)
    --parallel N        Maximum parallel samples (default: 4)
    --check-model       Check if Clair3 model is properly installed

Examples:
    $0 --test 10090                 # Test single sample
    $0 --all                        # Batch process all samples
    $0 10090 10893                  # Process specified samples
    $0 --list                       # List all samples
    $0 --report                     # Generate summary report
    $0 --check-model                # Check model installation

EOF
}

# Check environment and dependencies
check_environment() {
    log_info "Checking runtime environment..."
    
    # Check Clair3 environment
    if ! conda info --envs | grep -q "Clair3"; then
        log_error "Clair3 environment not found, please create environment first"
        log_info "Create via: conda create -n Clair3 -c bioconda clair3 samtools bcftools parallel"
        exit 1
    fi
    
    # Check necessary directories
    for dir in "$FILTERING_DIR" "$CLAIR3_PATH"; do
        if [[ ! -d "$dir" ]]; then
            log_error "Directory does not exist: $dir"
            exit 1
        fi
    done
    
    # Check model (warning only, not fatal)
    if [[ ! -d "$MODEL_DIR" ]]; then
        log_warn "Clair3 model directory may not exist: $MODEL_DIR"
        log_warn "Please ensure you have downloaded the model and updated MODEL_DIR in this script"
        log_warn "Download instructions: https://github.com/HKU-BAL/Clair3"
    fi
    
    # Create output directories
    mkdir -p "$OUTPUT_DIR" "$REPORT_DIR"
    
    log_info "Environment check complete"
}

# Check model installation
check_model() {
    log_info "Checking Clair3 model installation..."
    
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        log_error "Not in conda environment. Please activate Clair3 environment first:"
        log_info "  conda activate Clair3"
        exit 1
    fi
    
    log_info "Current conda environment: ${CONDA_PREFIX}"
    log_info "Expected model path: ${MODEL_DIR}"
    
    if [[ -d "$MODEL_DIR" ]]; then
        log_info "Model directory found!"
        log_info "Model contents:"
        ls -lh "$MODEL_DIR"
    else
        log_error "Model directory not found: $MODEL_DIR"
        log_info ""
        log_info "To download and install the model:"
        log_info "  cd ${CONDA_PREFIX}/bin/"
        log_info "  mkdir -p models"
        log_info "  cd models"
        log_info "  wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v500.tar.gz"
        log_info "  tar -xzf r1041_e82_400bps_sup_v500.tar.gz"
        log_info "  rm r1041_e82_400bps_sup_v500.tar.gz"
        log_info ""
        log_info "For more models, visit: https://github.com/HKU-BAL/Clair3"
        exit 1
    fi
}

# Get all available samples
get_all_samples() {
    find "$FILTERING_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort
}

# List all samples
list_samples() {
    log_info "Available sample list:"
    local samples=($(get_all_samples))
    local total=${#samples[@]}
    
    echo -e "${CYAN}Total samples: $total${NC}"
    echo "----------------------------------------"
    
    local count=0
    for sample in "${samples[@]}"; do
        count=$((count + 1))
        printf "%3d. %s\n" $count "$sample"
        
        # New line every 10 samples
        if [[ $((count % 10)) -eq 0 ]]; then
            echo ""
        fi
    done
    echo ""
}

# Check if sample files exist
check_sample_files() {
    local sample_id="$1"
    local bam_file="${FILTERING_DIR}/${sample_id}/${sample_id}.final.filtered.sorted.bam"
    local consensus_file="${CONSENSUS_DIR}/${sample_id}/consensus.fasta"
    
    if [[ ! -f "$bam_file" ]]; then
        log_error "Sample $sample_id: BAM file does not exist - $bam_file"
        return 1
    fi
    
    if [[ ! -f "$consensus_file" ]]; then
        log_error "Sample $sample_id: consensus file does not exist - $consensus_file"
        return 1
    fi
    
    return 0
}

# Process single sample variant detection
process_single_sample() {
    local sample_id="$1"
    local output_dir="${OUTPUT_DIR}/${sample_id}_clair3"
    local log_file="${output_dir}.log"
    
    log_info "Starting to process sample: $sample_id"
    
    # Check files
    if ! check_sample_files "$sample_id"; then
        return 1
    fi
    
    # Check if already completed
    if [[ -f "${output_dir}/merge_output.vcf.gz" ]] && [[ -f "${output_dir}/COMPLETED" ]]; then
        log_warn "Sample $sample_id already completed, skipping"
        return 0
    fi
    
    # Clean up previous incomplete output
    [[ -d "$output_dir" ]] && rm -rf "$output_dir"
    
    local bam_file="${FILTERING_DIR}/${sample_id}/${sample_id}.final.filtered.sorted.bam"
    local consensus_file="${CONSENSUS_DIR}/${sample_id}/consensus.fasta"
    
    # Record start time
    local start_time=$(date +%s)
    
    {
        echo "========================================="
        echo "Sample: $sample_id"
        echo "Start time: $(date)"
        echo "BAM file: $bam_file"
        echo "Reference sequence: $consensus_file"
        echo "Output directory: $output_dir"
        echo "========================================="
        
        # Activate environment and run Clair3
        source ~/miniforge3/etc/profile.d/conda.sh
        conda activate Clair3
        # Ensure GNU getopt is used (required for ARM Mac)
        export PATH="/opt/homebrew/opt/gnu-getopt/bin:$CONDA_PREFIX/bin:$CLAIR3_PATH:$PATH"
        cd "$CLAIR3_PATH"
        
        ./run_clair3.sh \
            --platform=ont \
            --threads="$THREADS_PER_SAMPLE" \
            --model_path="$MODEL_DIR" \
            --bam_fn="$bam_file" \
            --ref_fn="$consensus_file" \
            --output="$output_dir" \
            --include_all_ctgs \
            --snp_min_af="$SNP_MIN_AF" \
            --indel_min_af="$INDEL_MIN_AF" \
            --min_mq="$MIN_MQ" \
            --min_coverage="$MIN_COV" \
            --gvcf \
            --sample_name="$sample_id" \
            --pypy=python3 \
            --var_pct_full=1.0 \
            --ref_pct_full=1.0 \
            --haploid_sensitive \
            --no_phasing_for_fa \
            --enable_variant_calling_at_sequence_head_and_tail
        
        # Calculate runtime
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        echo "========================================="
        echo "Completion time: $(date)"
        echo "Runtime: ${duration} seconds"
        echo "========================================="
        
        # Create completion marker
        echo "COMPLETED at $(date)" > "${output_dir}/COMPLETED"
        echo "Duration: ${duration} seconds" >> "${output_dir}/COMPLETED"
        
    } > "$log_file" 2>&1
    
    if [[ $? -eq 0 ]]; then
        log_info "Sample $sample_id processing complete (${duration} seconds)"
        return 0
    else
        log_error "Sample $sample_id processing failed, check log: $log_file"
        return 1
    fi
}

# Generate detailed report for single sample
generate_sample_report() {
    local sample_id="$1"
    local output_dir="${OUTPUT_DIR}/${sample_id}_clair3"
    local report_file="${REPORT_DIR}/${sample_id}_report.md"
    
    if [[ ! -f "${output_dir}/merge_output.vcf.gz" ]]; then
        log_warn "Sample $sample_id has no valid output, skipping report generation"
        return 0
    fi
    
    log_info "Generating detailed report for sample $sample_id..."
    
    # Activate environment for VCF analysis
    cd "$output_dir"
    
    # Activate conda environment
    source ~/miniforge3/etc/profile.d/conda.sh
    conda activate Clair3 || {
        log_warn "Cannot activate Clair3 environment, using basic report"
        return 0
    }
    
    # Basic statistics
    local total_variants=$(bcftools view -H merge_output.vcf.gz | wc -l 2>/dev/null || echo "0")
    local pass_variants=$(bcftools view -f PASS -H merge_output.vcf.gz | wc -l 2>/dev/null || echo "0")  
    local real_variants=$(bcftools view -H merge_output.vcf.gz | grep -v 'RefCall' | wc -l 2>/dev/null || echo "0")
    
    # Classification statistics
    local snp_count=$(bcftools view -H merge_output.vcf.gz | grep -v 'RefCall' | awk 'length($4)==1 && length($5)==1' | wc -l 2>/dev/null || echo "0")
    local indel_count=$(bcftools view -H merge_output.vcf.gz | grep -v 'RefCall' | awk 'length($4)!=1 || length($5)!=1' | wc -l 2>/dev/null || echo "0")
    
    # Runtime information
    local duration="unknown"
    if [[ -f COMPLETED ]]; then
        duration=$(grep 'Duration:' COMPLETED | cut -d' ' -f2 2>/dev/null || echo "unknown")
    fi
    
    # Generate Markdown report
    cat > "$report_file" << EOF
# Sample $sample_id Variant Detection Report

## Basic Information
- **Sample ID**: $sample_id
- **Analysis time**: $(date)
- **Runtime**: ${duration} seconds
- **Output directory**: $PWD

## Detection Statistics
| Metric | Count |
|--------|-------|
| Total sites | $total_variants |
| PASS sites | $pass_variants |
| **Real variants** | **$real_variants** |
| SNP count | $snp_count |
| Indel count | $indel_count |

## File Outputs
- Main VCF file: \`merge_output.vcf.gz\`
- Pileup VCF: \`pileup.vcf.gz\`  
- Full-alignment VCF: \`full_alignment.vcf.gz\`
- gVCF file: \`merge_output.gvcf.gz\`

EOF

    # If there are real variants, add detailed information
    if [[ $real_variants -gt 0 ]]; then
        cat >> "$report_file" << 'DETAIL_EOF'

## Variant Details

| Position | Reference | Alternative | Frequency | Quality | Filter | Type |
|----------|-----------|-------------|-----------|---------|--------|------|
DETAIL_EOF
        
        bcftools query -f '%POS\t%REF\t%ALT\t[%AF]\t%QUAL\t%FILTER\t%INFO/P%INFO/F\n' merge_output.vcf.gz | grep -v 'RefCall' | while read pos ref alt af qual filter info; do
            # Determine variant type
            if [[ ${#ref} -eq 1 ]] && [[ ${#alt} -eq 1 ]]; then
                var_type='SNP'
            else
                var_type='Indel'
            fi
            
            # Format frequency
            af_percent=$(echo "$af" | awk '{printf "%.2f%%", $1*100}')
            
            echo "| $pos | $ref | $alt | $af_percent | $qual | $filter | $var_type |" >> "$report_file"
        done
    
    else
        cat >> "$report_file" << 'NOVAR_EOF'

## Variant Details
> **No real variants detected**
> 
> This sample is highly consistent with its own consensus sequence, no significant variants detected.
> This usually indicates high sample quality or a single clone.

NOVAR_EOF
    fi
    
    # Add quality assessment
    cat >> "$report_file" << 'ASSESS_EOF'

## Quality Assessment

ASSESS_EOF
    
    if [[ $real_variants -eq 0 ]]; then
        cat >> "$report_file" << 'HIGHQ_EOF'
- High-quality sample: No variants detected, strong sequence consistency
- Processing successful: All output files generated normally
- Data characteristics: Single clone or highly purified sample
HIGHQ_EOF
    else
        local max_af=$(bcftools query -f '[%AF]\n' merge_output.vcf.gz | grep -v 'RefCall' | sort -rn | head -1 | awk '{printf "%.1f%%", $1*100}' 2>/dev/null || echo "unknown")
        cat >> "$report_file" << VAR_EOF
- Successfully detected: Found $real_variants real variants
- Diverse sample: Quasispecies diversity or mixed infection present
- Variant frequency: Maximum $max_af
VAR_EOF
    fi
    
    log_info "Sample $sample_id report generation complete: $report_file"
}

# Batch process all samples
batch_process_all() {
    local samples=($(get_all_samples))
    local total=${#samples[@]}
    local success_count=0
    local failed_count=0
    local start_time=$(date +%s)
    
    log_info "Starting batch processing of $total samples (parallelism: $MAX_PARALLEL_SAMPLES)"
    
    # Create temporary task list
    local task_list="${OUTPUT_DIR}/batch_tasks.txt"
    printf '%s\n' "${samples[@]}" > "$task_list"
    
    # Use parallel for parallel processing
    if command -v parallel >/dev/null 2>&1; then
        log_info "Using GNU parallel for parallel processing..."
        
        export -f process_single_sample generate_sample_report log_info log_error log_warn check_sample_files
        export FILTERING_DIR CONSENSUS_DIR OUTPUT_DIR REPORT_DIR CLAIR3_PATH MODEL_DIR
        export THREADS_PER_SAMPLE SNP_MIN_AF INDEL_MIN_AF MIN_MQ MIN_COV
        export RED GREEN YELLOW BLUE PURPLE CYAN NC
        
        parallel -j "$MAX_PARALLEL_SAMPLES" --joblog "${OUTPUT_DIR}/parallel.log" \
            'process_single_sample {} && generate_sample_report {}' :::: "$task_list"
            
    else
        # Fallback to basic parallel processing
        log_warn "GNU parallel not found, using basic parallel processing"
        
        local running_jobs=0
        for sample in "${samples[@]}"; do
            # Control parallelism
            while [[ $running_jobs -ge $MAX_PARALLEL_SAMPLES ]]; do
                wait -n  # Wait for any background task to complete
                running_jobs=$((running_jobs - 1))
            done
            
            # Start new task
            {
                if process_single_sample "$sample"; then
                    generate_sample_report "$sample"
                fi
            } &
            
            running_jobs=$((running_jobs + 1))
            log_progress "Started: $sample (active tasks: $running_jobs)"
        done
        
        # Wait for all tasks to complete
        wait
    fi
    
    # Summarize results
    for sample in "${samples[@]}"; do
        if [[ -f "${OUTPUT_DIR}/${sample}_clair3/COMPLETED" ]]; then
            success_count=$((success_count + 1))
        else
            failed_count=$((failed_count + 1))
        fi
    done
    
    local end_time=$(date +%s)
    local total_duration=$((end_time - start_time))
    
    log_info "Batch processing complete!"
    log_info "Total samples: $total"
    log_info "Success: $success_count"
    log_info "Failed: $failed_count"
    log_info "Total time: $total_duration seconds"
    
    # Generate summary report
    generate_summary_report
}

# Generate summary report
generate_summary_report() {
    set +u  # Temporarily disable unbound variable check
    local summary_file="${REPORT_DIR}/SUMMARY_REPORT.md"
    
    log_info "Generating summary report..."
    
    conda run -n Clair3 bash -c "
        set +u  # Disable unbound variable check for this subshell
        # Count all sample results
        total_samples=0
        completed_samples=0
        total_variants=0
        samples_with_variants=0
        
        echo '# HBV Sample Variant Detection Summary Report' > '$summary_file'
        echo '' >> '$summary_file'
        echo 'Generated: $(date)' >> '$summary_file'
        echo '' >> '$summary_file'
        
        echo '## Processing Statistics' >> '$summary_file'
        
        for sample_dir in '${OUTPUT_DIR}'/*_clair3; do
            if [[ -d \"\$sample_dir\" ]]; then
                total_samples=\$((total_samples + 1))
                sample_id=\$(basename \"\$sample_dir\" | sed 's/_clair3\$//')
                
                if [[ -f \"\$sample_dir/COMPLETED\" ]] && [[ -f \"\$sample_dir/merge_output.vcf.gz\" ]]; then
                    completed_samples=\$((completed_samples + 1))
                    
                    # Count variants
                    variants=\$(bcftools view -H \"\$sample_dir/merge_output.vcf.gz\" | grep -v 'RefCall' | wc -l)
                    total_variants=\$((total_variants + variants))
                    
                    if [[ \$variants -gt 0 ]]; then
                        samples_with_variants=\$((samples_with_variants + 1))
                    fi
                fi
            fi
        done
        
        cat >> '$summary_file' << EOF

| Metric | Count | Percentage |
|--------|-------|------------|
| Total samples | \$total_samples | 100% |
| Completed samples | \$completed_samples | \$((\$completed_samples * 100 / \$total_samples))% |
| Samples with variants | \$samples_with_variants | \$((\$samples_with_variants * 100 / \$completed_samples))% |
| Total variants | \$total_variants | - |
| Average variants | \$((\$total_variants / \$completed_samples)) | - |

## Sample Details

| Sample ID | Status | Variants | Runtime | Report Link |
|-----------|--------|----------|---------|-------------|
EOF
        
        # Add detailed information for each sample
        for sample_dir in '${OUTPUT_DIR}'/*_clair3; do
            if [[ -d \"\$sample_dir\" ]]; then
                sample_id=\$(basename \"\$sample_dir\" | sed 's/_clair3\$//')
                
                if [[ -f \"\$sample_dir/COMPLETED\" ]] && [[ -f \"\$sample_dir/merge_output.vcf.gz\" ]]; then
                    variants=\$(bcftools view -H \"\$sample_dir/merge_output.vcf.gz\" | grep -v 'RefCall' | wc -l)
                    duration=\$(grep 'Duration:' \"\$sample_dir/COMPLETED\" | cut -d' ' -f2 2>/dev/null || echo 'unknown')
                    status='Completed'
                    report_link=\"[\${sample_id}_report.md](./${sample_id}_report.md)\"
                else
                    variants='N/A'
                    duration='N/A'
                    status='Failed'
                    report_link='N/A'
                fi
                
                echo \"| \$sample_id | \$status | \$variants | \${duration}s | \$report_link |\" >> '$summary_file'
            fi
        done
    "
    
    log_info "Summary report generation complete: $summary_file"
    set -u  # Re-enable unbound variable check
}

# Main function
main() {
    # Parse command line arguments
    local test_mode=false
    local batch_mode=false
    local list_mode=false
    local report_only=false
    local check_model_mode=false
    local test_sample=""
    local specific_samples=()
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -t|--test)
                test_mode=true
                test_sample="$2"
                shift 2
                ;;
            -a|--all)
                batch_mode=true
                shift
                ;;
            -l|--list)
                list_mode=true
                shift
                ;;
            -r|--report)
                report_only=true
                shift
                ;;
            --check-model)
                check_model_mode=true
                shift
                ;;
            --threads)
                THREADS_PER_SAMPLE="$2"
                shift 2
                ;;
            --parallel)
                MAX_PARALLEL_SAMPLES="$2"
                shift 2
                ;;
            -*)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
            *)
                specific_samples+=("$1")
                shift
                ;;
        esac
    done
    
    # Display title
    echo -e "${PURPLE}"
    echo "========================================="
    echo "    HBV Sample Batch Variant Detection System v2.0"
    echo "    (variants_call_10 pipeline)"
    echo "========================================="
    echo -e "${NC}"
    
    # Model check mode
    if [[ "$check_model_mode" == true ]]; then
        check_model
        exit 0
    fi
    
    check_environment
    
    # Execute corresponding function
    if [[ "$list_mode" == true ]]; then
        list_samples
        
    elif [[ "$report_only" == true ]]; then
        log_info "Report-only mode"
        generate_summary_report
        
    elif [[ "$test_mode" == true ]]; then
        log_info "Single sample test mode: $test_sample"
        if check_sample_files "$test_sample"; then
            if process_single_sample "$test_sample"; then
                generate_sample_report "$test_sample"
                log_info "Test complete! Check report: ${REPORT_DIR}/${test_sample}_report.md"
            fi
        fi
        
    elif [[ "$batch_mode" == true ]]; then
        log_info "Batch process all samples mode"
        batch_process_all
        
    elif [[ ${#specific_samples[@]} -gt 0 ]]; then
        log_info "Processing specified samples: ${specific_samples[*]}"
        for sample in "${specific_samples[@]}"; do
            if check_sample_files "$sample"; then
                if process_single_sample "$sample"; then
                    generate_sample_report "$sample"
                fi
            fi
        done
        generate_summary_report
        
    else
        log_error "Please specify operation mode"
        show_help
        exit 1
    fi
    
    echo -e "${GREEN}"
    echo "========================================="
    echo "           Task Execution Complete!"
    echo "========================================="
    echo -e "${NC}"
}

# Run main function
main "$@"
