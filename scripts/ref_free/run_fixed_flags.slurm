#!/usr/bin/bash
#SBATCH -J fixed_flags
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

set -euo pipefail

INPUT_DIR=${1:?input dir}
OUTPUT_BASE=${2:?output base}
TOTAL_REPEATS=${TOTAL_REPEATS:-1000000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-20000}
REF_N=${REF_N:-40}
SEED=${SEED:-42}
EP_THRESHOLD=${EP_THRESHOLD:-0.5}
EP_RECALL=${EP_RECALL:-0.65}
Z_THRESHOLD=${Z_THRESHOLD:-0.85}
Z_RECALL=${Z_RECALL:-0.95}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF=${EZ_CUTOFF:-3.0}
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
    echo "nothing to do"
    exit 0
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

echo "[$(date)] fixed_flags array=${SLURM_ARRAY_TASK_ID} [$REPEAT_START, $REPEAT_END)"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 ref_free_fixed_flags.py \
        --input-dir "$INPUT_DIR" \
        --output-base "$OUTPUT_BASE" \
        --total-repeats "$TOTAL_REPEATS" \
        --repeat-start "$REPEAT_START" \
        --repeat-end "$REPEAT_END" \
        --ref-n "$REF_N" \
        --seed "$SEED" \
        --cutoff "$CUTOFF" \
        --ez-cutoff "$EZ_CUTOFF" \
        --ep-threshold "$EP_THRESHOLD" \
        --ep-recall "$EP_RECALL" \
        --z-threshold "$Z_THRESHOLD" \
        --z-recall "$Z_RECALL"

echo "[$(date)] Finished array=${SLURM_ARRAY_TASK_ID}"
