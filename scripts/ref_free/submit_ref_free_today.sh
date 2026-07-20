#!/usr/bin/bash
# Submit today's clean 48+48 ref_free runs (fixed + filtered combos) and plots.
#
# Layout:
#   ${TODAY_BASE}/fixed_combo/{ref_free_ezscore,plots}
#   ${TODAY_BASE}/filtered_combos/{ref_free_ezscore,plots}

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
TODAY_BASE=${TODAY_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-ref_free}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
MAX_ARRAY_JOBS=${MAX_ARRAY_JOBS:-50}
REF_N=${REF_N:-48}
SEED=${SEED:-42}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
DRY_RUN=${DRY_RUN:-0}

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
if [ "$N_JOBS" -gt "$MAX_ARRAY_JOBS" ]; then
    REPEATS_PER_JOB=$(( (TOTAL_REPEATS + MAX_ARRAY_JOBS - 1) / MAX_ARRAY_JOBS ))
    N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
fi
ARRAY_LAST=$((N_JOBS - 1))

FIXED_BASE="${TODAY_BASE}/fixed_combo"
FILTERED_BASE="${TODAY_BASE}/filtered_combos"
mkdir -p "$FIXED_BASE" "$FILTERED_BASE"

echo "Today base     : $TODAY_BASE"
echo "Array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB} repeats/job)"

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] would submit fixed + filtered arrays and plot jobs"
    exit 0
fi

export TOTAL_REPEATS REPEATS_PER_JOB REF_N SEED SIF
export EZ_CUTOFF_MIN=3.0 EZ_CUTOFF_MAX=4.5 EZ_CUTOFF_STEP=0.1
export CUTOFF=3.0 MIN_FF=0

# --- fixed combo ---
export COMBO_MODE=fixed
export EP_THRESHOLD=0.5 EP_RECALL=0.65
export Z_THRESHOLD=0.85 Z_RECALL=0.95
unset EP_THRESHOLD_MIN EP_THRESHOLD_MAX EP_RECALL_MIN EP_RECALL_MAX
unset Z_THRESHOLD_MIN Z_THRESHOLD_MAX Z_RECALL_MIN Z_RECALL_MAX

fixed_job=$(sbatch --parsable --job-name=ref_free_fixed \
    --array="0-${ARRAY_LAST}" \
    run_ref_free_ezscore.slurm "$INPUT_DIR" "$FIXED_BASE")
echo "Submitted fixed array job_id=${fixed_job}"

fixed_plot=$(sbatch --parsable --job-name=plot_fixed \
    --dependency="afterok:${fixed_job}" \
    run_aggregate_and_plot.slurm \
    "$FIXED_BASE" \
    "48+48 fixed combo (ep 0.5/0.65, z 0.85/0.95)")
echo "Submitted fixed aggregate+plot job_id=${fixed_plot}"

# --- filtered combos ---
export COMBO_MODE=all
unset EP_THRESHOLD EP_RECALL Z_THRESHOLD Z_RECALL
export EP_THRESHOLD_MIN=0.33 EP_THRESHOLD_MAX=0.67
export EP_RECALL_MIN=0.5 EP_RECALL_MAX=0.75
export Z_THRESHOLD_MIN=0.8 Z_THRESHOLD_MAX=0.95
export Z_RECALL_MIN=0.9 Z_RECALL_MAX=0.99

filtered_job=$(sbatch --parsable --job-name=ref_free_filt \
    --array="0-${ARRAY_LAST}" \
    run_ref_free_ezscore.slurm "$INPUT_DIR" "$FILTERED_BASE")
echo "Submitted filtered array job_id=${filtered_job}"

filtered_plot=$(sbatch --parsable --job-name=plot_filt \
    --dependency="afterok:${filtered_job}" \
    run_aggregate_and_plot.slurm \
    "$FILTERED_BASE" \
    "48+48 filtered combos")
echo "Submitted filtered aggregate+plot job_id=${filtered_plot}"

echo "Outputs: $TODAY_BASE"
