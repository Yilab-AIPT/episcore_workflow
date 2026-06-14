#!/usr/bin/bash
#SBATCH -J random_ref40_scores
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x_%a.log
#SBATCH -e logs/%x_%a.log
#SBATCH --array=0-19

# Randomly select 40 dev+Normal reference samples (excluding PTAY0586P8S1),
# compute zscore/episcore/ezscore for 1000 unique draws.
#
# Array task k handles repeats [k*50, (k+1)*50).
#
# Usage:
#     sbatch run_random_ref40_scores.slurm [<input_dir>] [<output_base>]
#
# Defaults:
#     input_dir   : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40
#     output_base : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40}
OUTPUT_BASE=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

REPEAT_START=$((SLURM_ARRAY_TASK_ID * 50))
REPEAT_END=$((REPEAT_START + 50))
THREADS=${SLURM_CPUS_PER_TASK:-8}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] random_ref40_scores array=${SLURM_ARRAY_TASK_ID} host=$(hostname)"
echo "  input dir     : $INPUT_DIR"
echo "  output base   : $OUTPUT_BASE"
echo "  repeat range  : [$REPEAT_START, $REPEAT_END)"
echo "  threads       : $THREADS"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 calc_random_ref40_scores.py \
        --input-dir "$INPUT_DIR" \
        --output-base "$OUTPUT_BASE" \
        --total-repeats 1000 \
        --repeat-start "$REPEAT_START" \
        --repeat-end "$REPEAT_END" \
        --threads "$THREADS"

echo "[$(date)] Finished random_ref40_scores array=${SLURM_ARRAY_TASK_ID}"
