#!/usr/bin/bash
#SBATCH -J plot_auc_curve
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

set -euo pipefail

FLAG_DIR=${1:?flag dir …/fixed_flags}
OUT_DIR=${2:-}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
N_BOOT=${N_BOOT:-200}
FF_MIN=${FF_MIN:-0.01}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

ARGS=(--flag-dir "$FLAG_DIR" --ff-min "$FF_MIN" --n-boot "$N_BOOT")
if [[ -n "$OUT_DIR" ]]; then
    ARGS+=(--output-dir "$OUT_DIR")
fi

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 plot_auc_vs_repeats.py "${ARGS[@]}"

echo "[$(date)] Done"
