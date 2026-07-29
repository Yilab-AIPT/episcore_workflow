#!/usr/bin/bash
# Fast AUC-vs-repeats curve for fixed combo (dev/test only).
#
# Design: one 1e6-repeat run stores per-repeat flags; bootstrap over flags
# gives mean AUC + percentile band at each repeat_n — no re-simulation.
#
# Output:
#   ${OUT_BASE}/fixed_flags/flags_*.npz
#   ${OUT_BASE}/auc_learning_curve/auc_vs_repeats.html

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUT_BASE=${OUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260721-ref_free/fixed_auc_curve}
TOTAL_REPEATS=${TOTAL_REPEATS:-1000000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-20000}
MAX_ARRAY_JOBS=${MAX_ARRAY_JOBS:-50}
REF_N=${REF_N:-40}
SEED=${SEED:-42}
N_BOOT=${N_BOOT:-200}
FF_MIN=${FF_MIN:-0.01}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
if [ "$N_JOBS" -gt "$MAX_ARRAY_JOBS" ]; then
    REPEATS_PER_JOB=$(( (TOTAL_REPEATS + MAX_ARRAY_JOBS - 1) / MAX_ARRAY_JOBS ))
    N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
fi
ARRAY_LAST=$((N_JOBS - 1))

mkdir -p "$OUT_BASE"
echo "OUT_BASE       : $OUT_BASE"
echo "total repeats  : $TOTAL_REPEATS"
echo "array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB}/job)"
echo "bootstrap      : n_boot=$N_BOOT ff_min=$FF_MIN"

export TOTAL_REPEATS REPEATS_PER_JOB REF_N SEED SIF N_BOOT FF_MIN
export EP_THRESHOLD=0.5 EP_RECALL=0.65 Z_THRESHOLD=0.85 Z_RECALL=0.95
export CUTOFF=3.0 EZ_CUTOFF=3.0

arr=$(sbatch --parsable --job-name=fixed_flags \
    --array="0-${ARRAY_LAST}" \
    run_fixed_flags.slurm "$INPUT_DIR" "$OUT_BASE")
echo "Submitted flags array job_id=${arr}"

plot=$(sbatch --parsable --job-name=plot_auc_curve \
    --dependency="afterok:${arr}" \
    run_plot_auc_curve.slurm \
    "${OUT_BASE}/fixed_flags" \
    "${OUT_BASE}/auc_learning_curve")
echo "Submitted plot job_id=${plot}"
echo "HTML will be: ${OUT_BASE}/auc_learning_curve/auc_vs_repeats.html"
