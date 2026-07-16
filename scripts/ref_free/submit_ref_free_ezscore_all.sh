#!/usr/bin/bash
# Submit 48+48 dev Normal split ref_free_ezscore sweep (all combos).

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${OUTPUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260715-ref_free_ezscore_all}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
MAX_ARRAY_JOBS=${MAX_ARRAY_JOBS:-50}
REF_N=${REF_N:-48}
COMBO_MODE=${COMBO_MODE:-all}
EP_THRESHOLD=${EP_THRESHOLD:-}
EP_RECALL=${EP_RECALL:-}
Z_THRESHOLD=${Z_THRESHOLD:-}
Z_RECALL=${Z_RECALL:-}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF=${EZ_CUTOFF:-}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0}
DRY_RUN=${DRY_RUN:-0}
AGGREGATE=${AGGREGATE:-1}

cd "$SCRIPT_DIR"
mkdir -p logs

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
if [ "$N_JOBS" -gt "$MAX_ARRAY_JOBS" ]; then
    REPEATS_PER_JOB=$(( (TOTAL_REPEATS + MAX_ARRAY_JOBS - 1) / MAX_ARRAY_JOBS ))
    N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
fi
ARRAY_LAST=$((N_JOBS - 1))

echo "Input dir      : $INPUT_DIR"
echo "Output base    : $OUTPUT_BASE"
echo "Total repeats  : $TOTAL_REPEATS"
echo "Array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB} repeats each)"
echo "Ref split      : ${REF_N}+${REF_N} dev Normal"
echo "Combo mode     : $COMBO_MODE"
if [ "$COMBO_MODE" = "fixed" ]; then
    echo "Episcore combo : threshold=$EP_THRESHOLD recall=$EP_RECALL"
    echo "Zscore combo   : threshold=$Z_THRESHOLD recall=$Z_RECALL"
fi

echo "Ep/z cutoff    : $CUTOFF"
if [ -n "$EZ_CUTOFF" ]; then
    echo "Ez cutoff      : $EZ_CUTOFF"
else
    echo "Ez cutoff      : $CUTOFF (same as ep/z)"
fi

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --array=0-${ARRAY_LAST} run_ref_free_ezscore.slurm '$INPUT_DIR' '$OUTPUT_BASE'"
    exit 0
fi

array_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" REPEATS_PER_JOB="$REPEATS_PER_JOB" \
    REF_N="$REF_N" COMBO_MODE="$COMBO_MODE" \
    EP_THRESHOLD="$EP_THRESHOLD" EP_RECALL="$EP_RECALL" \
    Z_THRESHOLD="$Z_THRESHOLD" Z_RECALL="$Z_RECALL" \
    CUTOFF="$CUTOFF" EZ_CUTOFF="$EZ_CUTOFF" SEED="$SEED" MIN_FF="$MIN_FF" \
    sbatch --parsable --job-name=ref_free_ezscore \
    --array="0-${ARRAY_LAST}" run_ref_free_ezscore.slurm "$INPUT_DIR" "$OUTPUT_BASE")
echo "Submitted ref_free_ezscore array job_id=${array_jobid}"

if [ "$AGGREGATE" = 1 ]; then
    agg_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" sbatch --parsable \
        --job-name=aggregate_ref_free_ez \
        --dependency="afterok:${array_jobid}" \
        run_aggregate_ref_free_ezscore.slurm "$OUTPUT_BASE")
    echo "Submitted aggregate_ref_free_ez job_id=${agg_jobid}"
fi
