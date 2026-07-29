#!/usr/bin/bash
#SBATCH -J agg_plot_ref_free
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log

set -euo pipefail

OUTPUT_BASE=${1:?OUTPUT_BASE required}
TITLE=${2:-40+40 reference-free signal ratio}
FF_MIN=${3:-0.01}
DO_PAIR_SEARCH=${4:-0}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] Aggregate + plot for $OUTPUT_BASE (ff-min=$FF_MIN pair_search=$DO_PAIR_SEARCH)"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 aggregate_ref_free_ezscore.py --output-base "$OUTPUT_BASE" --ff-min "$FF_MIN"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 plot_ref_free_interactive.py \
        --result-dir "$OUTPUT_BASE" \
        --title "$TITLE" \
        --ff-min "$FF_MIN"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 plot_ref_free_dual_ezscore.py \
        --result-dir "$OUTPUT_BASE" \
        --title "$TITLE" \
        --ff-min "$FF_MIN"

if [ "$DO_PAIR_SEARCH" = "1" ]; then
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
        python3 search_ez_pair_subsets.py --output-base "$OUTPUT_BASE" --ff-min "$FF_MIN"
    # dual + 3-panel for best subset
    BEST_TSV="$OUTPUT_BASE/ref_free_ezscore/abnormality_signal_ratio_best_subset.tsv"
    if [ -f "$BEST_TSV" ]; then
        singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
            python3 plot_ref_free_interactive.py \
                --result-dir "$OUTPUT_BASE" \
                --scores-tsv "$BEST_TSV" \
                --output-html "$OUTPUT_BASE/plots/signal_ratio_best_subset.html" \
                --title "$TITLE (best ez-pair subset)" \
                --ff-min "$FF_MIN"
        singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
            python3 plot_ref_free_dual_ezscore.py \
                --result-dir "$OUTPUT_BASE" \
                --scores-tsv "$BEST_TSV" \
                --output-html "$OUTPUT_BASE/plots/ezscore_eval_vs_val_best_subset.html" \
                --title "$TITLE (best ez-pair subset)" \
                --ff-min "$FF_MIN"
    fi
fi

echo "[$(date)] Done"
