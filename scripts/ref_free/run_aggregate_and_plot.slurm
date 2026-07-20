#!/usr/bin/bash
#SBATCH -J agg_plot_ref_free
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

set -euo pipefail

OUTPUT_BASE=${1:?OUTPUT_BASE required}
TITLE=${2:-48+48 reference-free signal ratio}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] Aggregate + plot for $OUTPUT_BASE"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 aggregate_ref_free_ezscore.py --output-base "$OUTPUT_BASE"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 plot_ref_free_interactive.py \
        --result-dir "$OUTPUT_BASE" \
        --title "$TITLE"

echo "[$(date)] Done"
