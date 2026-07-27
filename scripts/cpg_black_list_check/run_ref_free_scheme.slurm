#!/usr/bin/bash
#SBATCH -J bl_ref_free
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log

# Array task: one repeat slice for one scheme.
# Env: SCHEME, TOTAL_REPEATS, REPEATS_PER_JOB, REF_N, SEED, MIN_FF, EZ_*, CUTOFF
#
# Usage (via sbatch --export):
#   sbatch --array=0-49 --export=ALL,SCHEME=baseline run_ref_free_scheme.slurm \
#       <meta_csv> <pipeline_root> <scheme_output_base>

set -euo pipefail

META_CSV=${1:?meta.csv}
PIPELINE_ROOT=${2:?pipeline root}
OUTPUT_BASE=${3:?scheme output base}

SCHEME=${SCHEME:?SCHEME required}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
REF_N=${REF_N:-40}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0.01}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF_MIN=${EZ_CUTOFF_MIN:-3.0}
EZ_CUTOFF_MAX=${EZ_CUTOFF_MAX:-4.5}
EZ_CUTOFF_STEP=${EZ_CUTOFF_STEP:-0.1}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}

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

echo "[$(date)] scheme=${SCHEME} array=${SLURM_ARRAY_TASK_ID} repeats=[${REPEAT_START},${REPEAT_END})"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 "${SCRIPT_DIR}/ref_free_scheme.py" \
        --meta-csv "$META_CSV" \
        --pipeline-root "$PIPELINE_ROOT" \
        --scheme "$SCHEME" \
        --output-base "$OUTPUT_BASE" \
        --total-repeats "$TOTAL_REPEATS" \
        --repeat-start "$REPEAT_START" \
        --repeat-end "$REPEAT_END" \
        --ref-n "$REF_N" \
        --seed "$SEED" \
        --cutoff "$CUTOFF" \
        --min-ff "$MIN_FF" \
        --ez-cutoff-min "$EZ_CUTOFF_MIN" \
        --ez-cutoff-max "$EZ_CUTOFF_MAX" \
        --ez-cutoff-step "$EZ_CUTOFF_STEP"

echo "[$(date)] Finished scheme=${SCHEME} array=${SLURM_ARRAY_TASK_ID}"
