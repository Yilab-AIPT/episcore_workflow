#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/results/grid_search/logs"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="${LOG_DIR}/grid_search_${TIMESTAMP}.log"

PROFILE="${PROFILE:-at_ref,alioth_slurm}"
NF_OPTS="${NF_OPTS:-}"
DRY_RUN="${DRY_RUN:-false}"

THRESHOLDS=(0.1 0.33 0.5 0.7 0.9)
CPG_RECALLS=(0.25 0.45 0.65 0.85)

mkdir -p "$LOG_DIR"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG_FILE"; }

log "=== Grid search started ==="
log "Profile: ${PROFILE}"
log "Thresholds: ${THRESHOLDS[*]}"
log "CpG recalls: ${CPG_RECALLS[*]}"
log "Working directory: ${SCRIPT_DIR}"

total=$(( ${#THRESHOLDS[@]} * ${#CPG_RECALLS[@]} ))
passed=0
failed=0
skipped=0
failed_combos=()

for recall in "${CPG_RECALLS[@]}"; do
    cpg_file="${SCRIPT_DIR}/assets/grid_search/CpG_recall${recall}.txt"
    if [[ ! -f "$cpg_file" ]]; then
        log "ERROR: CpG list not found: ${cpg_file}"
        log "Skipping all runs with CpG_recall=${recall}"
        skipped=$(( skipped + ${#THRESHOLDS[@]} ))
        continue
    fi
done

for threshold in "${THRESHOLDS[@]}"; do
    for recall in "${CPG_RECALLS[@]}"; do
        cpg_file="${SCRIPT_DIR}/assets/grid_search/CpG_recall${recall}.txt"
        outdir="/lustre1/cqyi/syfan/nf_autogluon_infer/results/220k_panel/grid_search/thres_${threshold}_cpg_${recall}_ref"
        run_label="threshold=${threshold}, CpG_recall=${recall}"

        if [[ ! -f "$cpg_file" ]]; then
            skipped=$(( skipped + 1 ))
            continue
        fi

        log "--- [$(( passed + failed + skipped + 1 ))/${total}] Running: ${run_label} ---"
        log "  cpg_list : ${cpg_file}"
        log "  outdir   : ${outdir}"

        mkdir -p "$outdir"

        nf_cmd=(
            nextflow run "${SCRIPT_DIR}/main.nf"
            -profile "$PROFILE"
            -work-dir "${outdir}/work"
            --outdir "$outdir"
            --threshold "$threshold"
            --cpg_list "$cpg_file"
            $NF_OPTS
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
            log "  Cleaning up work dir: ${outdir}/work"
            rm -rf "${outdir}/work"
        else
            log "  FAILED: ${run_label} (see ${run_log})"
            failed=$(( failed + 1 ))
            failed_combos+=("${run_label}")
        fi
    done
done

log "=== Grid search complete ==="
log "Total: ${total} | Passed: ${passed} | Failed: ${failed} | Skipped: ${skipped}"

if (( ${#failed_combos[@]} > 0 )); then
    log "Failed combinations:"
    for combo in "${failed_combos[@]}"; do
        log "  - ${combo}"
    done
    exit 1
fi

exit 0
