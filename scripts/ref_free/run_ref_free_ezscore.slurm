#!/usr/bin/bash
#SBATCH -J ref_free_ezscore
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${2:?OUTPUT_BASE required}
TOTAL_REPEATS=${TOTAL_REPEATS:-100000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-2000}
REF_N=${REF_N:-40}
COMBO_MODE=${COMBO_MODE:-all}
EP_THRESHOLD=${EP_THRESHOLD:-}
EP_RECALL=${EP_RECALL:-}
Z_THRESHOLD=${Z_THRESHOLD:-}
Z_RECALL=${Z_RECALL:-}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF_MIN=${EZ_CUTOFF_MIN:-3.0}
EZ_CUTOFF_MAX=${EZ_CUTOFF_MAX:-4.5}
EZ_CUTOFF_STEP=${EZ_CUTOFF_STEP:-0.1}
EP_THRESHOLD_MIN=${EP_THRESHOLD_MIN:-}
EP_THRESHOLD_MAX=${EP_THRESHOLD_MAX:-}
EP_RECALL_MIN=${EP_RECALL_MIN:-}
EP_RECALL_MAX=${EP_RECALL_MAX:-}
Z_THRESHOLD_MIN=${Z_THRESHOLD_MIN:-}
Z_THRESHOLD_MAX=${Z_THRESHOLD_MAX:-}
Z_RECALL_MIN=${Z_RECALL_MIN:-}
Z_RECALL_MAX=${Z_RECALL_MAX:-}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0}
STORE_PAIR_COUNTS=${STORE_PAIR_COUNTS:-0}
COMPRESS=${COMPRESS:-1}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

REPEAT_START=$((SLURM_ARRAY_TASK_ID * REPEATS_PER_JOB))
REPEAT_END=$((REPEAT_START + REPEATS_PER_JOB))
if [ "$REPEAT_END" -gt "$TOTAL_REPEATS" ]; then
    REPEAT_END=$TOTAL_REPEATS
fi
if [ "$REPEAT_START" -ge "$TOTAL_REPEATS" ]; then
    echo "[$(date)] array=${SLURM_ARRAY_TASK_ID}: nothing to do"
    exit 0
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

echo "[$(date)] ref_free_ezscore array=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "  output base : $OUTPUT_BASE"
echo "  repeats     : [$REPEAT_START, $REPEAT_END) of $TOTAL_REPEATS"
echo "  combo-mode  : $COMBO_MODE store_pair_counts=$STORE_PAIR_COUNTS compress=$COMPRESS"

PY_ARGS=(
    --input-dir "$INPUT_DIR"
    --output-base "$OUTPUT_BASE"
    --total-repeats "$TOTAL_REPEATS"
    --repeat-start "$REPEAT_START"
    --repeat-end "$REPEAT_END"
    --ref-n "$REF_N"
    --seed "$SEED"
    --cutoff "$CUTOFF"
    --min-ff "$MIN_FF"
    --combo-mode "$COMBO_MODE"
    --ez-cutoff-min "$EZ_CUTOFF_MIN"
    --ez-cutoff-max "$EZ_CUTOFF_MAX"
    --ez-cutoff-step "$EZ_CUTOFF_STEP"
)
if [ "$COMPRESS" = "1" ]; then
    PY_ARGS+=(--compress)
else
    PY_ARGS+=(--no-compress)
fi
if [ "$STORE_PAIR_COUNTS" = "1" ]; then
    PY_ARGS+=(--store-pair-counts)
else
    PY_ARGS+=(--no-store-pair-counts)
fi
if [ "$COMBO_MODE" = "fixed" ]; then
    PY_ARGS+=(
        --ep-threshold "$EP_THRESHOLD"
        --ep-recall "$EP_RECALL"
        --z-threshold "$Z_THRESHOLD"
        --z-recall "$Z_RECALL"
    )
else
    [[ -n "$EP_THRESHOLD_MIN" ]] && PY_ARGS+=(--ep-threshold-min "$EP_THRESHOLD_MIN")
    [[ -n "$EP_THRESHOLD_MAX" ]] && PY_ARGS+=(--ep-threshold-max "$EP_THRESHOLD_MAX")
    [[ -n "$EP_RECALL_MIN" ]] && PY_ARGS+=(--ep-recall-min "$EP_RECALL_MIN")
    [[ -n "$EP_RECALL_MAX" ]] && PY_ARGS+=(--ep-recall-max "$EP_RECALL_MAX")
    [[ -n "$Z_THRESHOLD_MIN" ]] && PY_ARGS+=(--z-threshold-min "$Z_THRESHOLD_MIN")
    [[ -n "$Z_THRESHOLD_MAX" ]] && PY_ARGS+=(--z-threshold-max "$Z_THRESHOLD_MAX")
    [[ -n "$Z_RECALL_MIN" ]] && PY_ARGS+=(--z-recall-min "$Z_RECALL_MIN")
    [[ -n "$Z_RECALL_MAX" ]] && PY_ARGS+=(--z-recall-max "$Z_RECALL_MAX")
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 ref_free_ezscore.py "${PY_ARGS[@]}"

echo "[$(date)] Finished ref_free_ezscore array=${SLURM_ARRAY_TASK_ID}"
