#!/usr/bin/bash
#SBATCH -J aggregate_ref_free
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

set -euo pipefail

OUTPUT_BASE=${1:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260713-ref_free}
TOTAL_REPEATS=${TOTAL_REPEATS:-}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] aggregate_ref_free host=$(hostname)"
echo "  output base : $OUTPUT_BASE"

AGG_ARGS=(--output-base "$OUTPUT_BASE")
if [ -n "$TOTAL_REPEATS" ]; then
    AGG_ARGS+=(--total-repeats "$TOTAL_REPEATS")
    echo "  total reps  : $TOTAL_REPEATS"
else
    echo "  total reps  : from run_config.json"
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 aggregate_ref_free.py "${AGG_ARGS[@]}"

echo "[$(date)] Finished aggregate_ref_free"
