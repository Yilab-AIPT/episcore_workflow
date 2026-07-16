#!/usr/bin/bash
#SBATCH -J ref_free
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

# Reference-free episcore/zscore abnormality sweep (all parquet combos).
# One array task processes a chunk of repeats (default 200 per task, 50 jobs).
#
# Usage:
#     sbatch run_ref_free.slurm [<input_dir>] [<output_base>]
#
# Environment (optional):
#     TOTAL_REPEATS    default 10000
#     REPEATS_PER_JOB  default 200
#     REF_N            default 50
#     COMBO_MODE       default all  (all | fixed)
#     EP_THRESHOLD     required when COMBO_MODE=fixed
#     EP_RECALL        required when COMBO_MODE=fixed
#     Z_THRESHOLD      required when COMBO_MODE=fixed
#     Z_RECALL         required when COMBO_MODE=fixed
#     CUTOFF           default 3.0
#     SEED             default 42
#     MIN_FF           default 0

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260713-ref_free}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
REF_N=${REF_N:-50}
COMBO_MODE=${COMBO_MODE:-all}
EP_THRESHOLD=${EP_THRESHOLD:-}
EP_RECALL=${EP_RECALL:-}
Z_THRESHOLD=${Z_THRESHOLD:-}
Z_RECALL=${Z_RECALL:-}
CUTOFF=${CUTOFF:-3.0}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0}
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
    echo "[$(date)] array=${SLURM_ARRAY_TASK_ID}: repeat_start=${REPEAT_START} >= total=${TOTAL_REPEATS}, nothing to do"
    exit 0
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] ref_free array=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "  input dir       : $INPUT_DIR"
echo "  output base     : $OUTPUT_BASE"
echo "  repeats / chunk : $REPEATS_PER_JOB"
echo "  repeat range    : [$REPEAT_START, $REPEAT_END) of $TOTAL_REPEATS"
echo "  ref-n           : $REF_N"
echo "  combo-mode      : $COMBO_MODE"
if [ "$COMBO_MODE" = "fixed" ]; then
    echo "  episcore combo  : threshold=$EP_THRESHOLD recall=$EP_RECALL"
    echo "  zscore combo    : threshold=$Z_THRESHOLD recall=$Z_RECALL"
else
    echo "  combos          : all available in parquets"
fi
echo "  cutoff          : $CUTOFF"

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
)
if [ "$COMBO_MODE" = "fixed" ]; then
    PY_ARGS+=(
        --ep-threshold "$EP_THRESHOLD"
        --ep-recall "$EP_RECALL"
        --z-threshold "$Z_THRESHOLD"
        --z-recall "$Z_RECALL"
    )
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 ref_free.py "${PY_ARGS[@]}"

echo "[$(date)] Finished ref_free array=${SLURM_ARRAY_TASK_ID}"
