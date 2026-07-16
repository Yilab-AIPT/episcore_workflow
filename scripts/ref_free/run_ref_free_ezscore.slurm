#!/usr/bin/bash
#SBATCH -J ref_free_ezscore
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260715-ref_free_ezscore_all}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
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

echo "[$(date)] ref_free_ezscore array=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "  input dir       : $INPUT_DIR"
echo "  output base     : $OUTPUT_BASE"
echo "  repeat range    : [$REPEAT_START, $REPEAT_END) of $TOTAL_REPEATS"
echo "  ref split       : ${REF_N}+${REF_N} dev Normal"
echo "  combo-mode      : $COMBO_MODE"

echo "  ep/z cutoff     : $CUTOFF"
if [ -n "$EZ_CUTOFF" ]; then
    echo "  ez cutoff       : $EZ_CUTOFF"
else
    echo "  ez cutoff       : $CUTOFF (same as ep/z)"
fi

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
if [ -n "$EZ_CUTOFF" ]; then
    PY_ARGS+=(--ez-cutoff "$EZ_CUTOFF")
fi
if [ "$COMBO_MODE" = "fixed" ]; then
    PY_ARGS+=(
        --ep-threshold "$EP_THRESHOLD"
        --ep-recall "$EP_RECALL"
        --z-threshold "$Z_THRESHOLD"
        --z-recall "$Z_RECALL"
    )
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 ref_free_ezscore.py "${PY_ARGS[@]}"

echo "[$(date)] Finished ref_free_ezscore array=${SLURM_ARRAY_TASK_ID}"
