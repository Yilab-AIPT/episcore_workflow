#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY_RUN="${DRY_RUN:-false}"

NF_MAIN="/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf"
# SAMPLESHEET="/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/samplesheet.csv"
# RESULTS_BASE="/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search"
SAMPLESHEET="/lustre1/cqyi/syfan/nf_methylbert_infer/data/220k_data_mb_vs_mq/infered_samplesheet.csv"
RESULTS_BASE="/lustre1/cqyi/syfan/nf_methylbert_infer/old_model_downsampled_episcore"

LOG_DIR="${RESULTS_BASE}/logs"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="${LOG_DIR}/grid_search_${TIMESTAMP}.log"

THRESHOLDS=(0.1 0.33 0.5 0.67 0.9)

mkdir -p "$LOG_DIR"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG_FILE"; }

log "=== Grid search started ==="
log "Profile: grid_search,singularity"
log "Thresholds: ${THRESHOLDS[*]}"
log "Samplesheet: ${SAMPLESHEET}"
log "Working directory: ${SCRIPT_DIR}"

total=${#THRESHOLDS[@]}
passed=0
failed=0
failed_runs=()

for threshold in "${THRESHOLDS[@]}"; do
    outdir="${RESULTS_BASE}/threshold_${threshold}"
    run_label="threshold=${threshold}"

    log "--- [$(( passed + failed + 1 ))/${total}] Running: ${run_label} ---"
    log "  outdir : ${outdir}"

    mkdir -p "$outdir"

    nf_cmd=(
        nextflow run "$NF_MAIN"
        -profile grid_search,singularity
        --input "$SAMPLESHEET"
        --outdir "$outdir"
        --threshold "$threshold"
    )

    if [[ "$DRY_RUN" == "true" ]]; then
        log "  [DRY RUN] ${nf_cmd[*]}"
        passed=$(( passed + 1 ))
        continue
    fi

    run_log="${outdir}/nextflow_run.log"
    if "${nf_cmd[@]}" 2>&1 | tee "$run_log"; then
        log "  PASSED: ${run_label}"
        passed=$(( passed + 1 ))
    else
        log "  FAILED: ${run_label} (see ${run_log})"
        failed=$(( failed + 1 ))
        failed_runs+=("${run_label}")
    fi
done

log "=== Grid search complete ==="
log "Total: ${total} | Passed: ${passed} | Failed: ${failed}"

if (( ${#failed_runs[@]} > 0 )); then
    log "Failed runs:"
    for label in "${failed_runs[@]}"; do
        log "  - ${label}"
    done
    exit 1
fi

exit 0
