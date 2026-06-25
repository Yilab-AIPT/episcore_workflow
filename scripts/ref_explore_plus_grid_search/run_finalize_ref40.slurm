#!/usr/bin/bash
#SBATCH -J finalize_ref40
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Post-aggregate ref-40 finalization: annotate scores, confusion PNGs, ezscore ref search.
#
# Usage:
#     sbatch run_finalize_ref40.slurm <output_base> [<input_dir>]
#
# Environment (optional):
#     N_EZSCORE_REF   default 20
#     N_REPEATS       default 5000
#     SEED            default 42
#     SKIP_EZSCORE    set to 1 to skip ezscore reference search

set -euo pipefail

OUTPUT_BASE=${1:?output_base required}
INPUT_DIR=${2:-}
N_EZSCORE_REF=${N_EZSCORE_REF:-20}
N_REPEATS=${N_REPEATS:-5000}
SEED=${SEED:-42}
SKIP_EZSCORE=${SKIP_EZSCORE:-0}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] finalize_ref40 host=$(hostname)"
echo "  output base   : $OUTPUT_BASE"
echo "  input dir     : ${INPUT_DIR:-<none>}"
echo "  n-ezscore-ref : $N_EZSCORE_REF"
echo "  n-repeats     : $N_REPEATS"
echo "  seed          : $SEED"
echo "  skip ezscore  : $SKIP_EZSCORE"

args=(
    --output-base "$OUTPUT_BASE"
    --n-ezscore-ref "$N_EZSCORE_REF"
    --n-repeats "$N_REPEATS"
    --seed "$SEED"
)
if [ -n "$INPUT_DIR" ]; then
    args+=(--input-dir "$INPUT_DIR")
fi
if [ "$SKIP_EZSCORE" = 1 ]; then
    args+=(--skip-ezscore-search)
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 finalize_ref40.py "${args[@]}"

echo "[$(date)] Finished finalize_ref40"
