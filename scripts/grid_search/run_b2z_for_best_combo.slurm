#!/usr/bin/bash
#SBATCH -J b2z_best_combo
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=8
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Recompute z_intra / z_inter / s_inter using the best (threshold, recall) per chr.
#
# Usage:
#     sbatch run_b2z_for_best_combo.slurm \
#            [<best_combo_csv>] [<output_base>]
#
# Defaults match the 20260508-grid_search run:
#     best_combo_csv : <RESULT_BASE>/best_combo_per_chr.csv
#     output_base    : <RESULT_BASE>/output
#
# Inputs:
#     best_combo_csv : CSV with columns chr,threshold,recall (extra columns ignored)
#     output_base    : Directory containing threshold_<t>_recall_<r>/ subdirs
#                      with _analyze_zscore.tsv.gz and _reference_zscore.tsv.gz
#                      (i.e. the output dir from beta_to_zscore.py).
#
# Outputs (under <output_base>/best_combo/):
#     _reference_zscore.tsv.gz
#     _analyze_zscore.tsv.gz

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search

BEST_COMBO_CSV=${1:-"${RESULT_BASE}/best_combo_per_chr.csv"}
OUTPUT_BASE=${2:-"${RESULT_BASE}/output"}

# Run from the directory containing b2z_for_best_combo.py (where sbatch was called).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

if [ ! -f "$BEST_COMBO_CSV" ]; then
    echo "ERROR: best-combo CSV not found: $BEST_COMBO_CSV" >&2
    exit 1
fi
if [ ! -d "$OUTPUT_BASE" ]; then
    echo "ERROR: output base dir not found: $OUTPUT_BASE" >&2
    exit 1
fi

echo "[$(date)] Starting b2z_for_best_combo on host $(hostname)"
echo "  best-combo CSV : $BEST_COMBO_CSV"
echo "  output base    : $OUTPUT_BASE"
echo "  output dir     : $OUTPUT_BASE/best_combo"
echo "  cpus           : ${SLURM_CPUS_PER_TASK:-8}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 b2z_for_best_combo.py \
        --best-combo-csv "$BEST_COMBO_CSV" \
        --output-base    "$OUTPUT_BASE"

echo "[$(date)] Finished b2z_for_best_combo"
