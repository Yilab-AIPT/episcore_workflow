#!/usr/bin/bash
#SBATCH -J bl_rf_rec_agg
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

# Aggregate baseline-0.6 / baseline-0.65 ref_free slices and plot.
#
# Usage:
#   sbatch run_baseline_recall_aggregate.slurm <ref_free_root>

set -euo pipefail

REF_FREE_ROOT=${1:?ref_free root}
EZ_CUTOFF=${EZ_CUTOFF:-4.5}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

SCHEMES=(baseline-0.6 baseline-0.65 baseline-deeper-0.65)
PLOT_DIR="${REF_FREE_ROOT}/plots_baseline_recall"

echo "[$(date)] Aggregate + plot under $REF_FREE_ROOT"
for scheme in "${SCHEMES[@]}"; do
    base="${REF_FREE_ROOT}/${scheme}"
    if [ ! -d "${base}/ref_free_ezscore" ]; then
        echo "SKIP missing ${base}/ref_free_ezscore"
        continue
    fi
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
        python3 "${SCRIPT_DIR}/aggregate_ref_free_scheme.py" \
            --output-base "$base"
done

mkdir -p "$PLOT_DIR"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 "${SCRIPT_DIR}/plot_baseline_recall_scatter.py" \
        --ref-free-root "$REF_FREE_ROOT" \
        --output-dir "$PLOT_DIR" \
        --ez-cutoff "$EZ_CUTOFF"

echo "[$(date)] Done → $PLOT_DIR"
