#!/usr/bin/bash
#SBATCH -J bl_rf_agg
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

# Aggregate all schemes then draw ff vs signal-ratio scatters.
#
# Usage:
#   sbatch run_ref_free_aggregate_plot.slurm <ref_free_root>

set -euo pipefail

REF_FREE_ROOT=${1:?ref_free root containing <scheme>/}
EZ_CUTOFF=${EZ_CUTOFF:-4.5}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

SCHEMES=(baseline 15-site1-J 15-site1-M 20-site1-M 15-site2-J 15-site2-M 20-site2-M)

echo "[$(date)] Aggregate schemes under $REF_FREE_ROOT"
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

PLOT_DIR="${REF_FREE_ROOT}/plots"
mkdir -p "$PLOT_DIR"
echo "[$(date)] Scatter plots → $PLOT_DIR (ez_cutoff=${EZ_CUTOFF})"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 "${SCRIPT_DIR}/plot_ref_free_scatter.py" \
        --ref-free-root "$REF_FREE_ROOT" \
        --output-dir "$PLOT_DIR" \
        --ez-cutoff "$EZ_CUTOFF"

echo "[$(date)] Done"
