#!/usr/bin/bash
#SBATCH -J grid_search_ref40
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

# Random ref-40 episcore/zscore/ezscore grid search.
# One array task processes a chunk of repeats (default 100 per task).
#
# Usage:
#     sbatch run_grid_search_ref40.slurm [<input_dir>] [<output_base>]
#
# Environment (optional):
#     TOTAL_REPEATS    default 100
#     REPEATS_PER_JOB  default 100  (repeats per array task)
#     SPLIT_MODE       default dev_test_split  (dev_test_split | even_split | all)
#     COMBO_MODE       default flexible        (flexible | fix)
#     MIN_FF           default 0
#
# Defaults:
#     input_dir   : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
#     output_base : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng}
TOTAL_REPEATS=${TOTAL_REPEATS:-100}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-100}
SPLIT_MODE=${SPLIT_MODE:-dev_test_split}
COMBO_MODE=${COMBO_MODE:-flexible}
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

echo "[$(date)] grid_search_ref40 array=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "  input dir       : $INPUT_DIR"
echo "  output base     : $OUTPUT_BASE"
echo "  repeats / chunk : $REPEATS_PER_JOB"
echo "  repeat range    : [$REPEAT_START, $REPEAT_END) of $TOTAL_REPEATS"
echo "  split-mode      : $SPLIT_MODE"
echo "  combo-mode      : $COMBO_MODE"
echo "  min-ff          : $MIN_FF"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 grid_search_ref40.py \
        --input-dir "$INPUT_DIR" \
        --output-base "$OUTPUT_BASE" \
        --total-repeats "$TOTAL_REPEATS" \
        --repeat-start "$REPEAT_START" \
        --repeat-end "$REPEAT_END" \
        --split-mode "$SPLIT_MODE" \
        --combo-mode "$COMBO_MODE" \
        --min-ff "$MIN_FF"

echo "[$(date)] Finished grid_search_ref40 array=${SLURM_ARRAY_TASK_ID}"
