#!/usr/bin/bash
#SBATCH -J b2z_global_recall
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Stitch per-chr z-score columns across combos under a single shared recall.
#
# Usage:
#     sbatch run_b2z_for_global_recall.slurm \
#            [<best_combo_csv>] [<output_base>]
#
# Defaults match the 20260508-grid_search run:
#     best_combo_csv : <RESULT_BASE>/best_combo_per_chr_global.csv
#     output_base    : <RESULT_BASE>/output
#
# Inputs:
#     best_combo_csv : CSV with columns chr,threshold,recall (extra columns ignored).
#                      Rows normally share the same recall.
#     output_base    : Directory containing threshold_<t>_recall_<r>/ subdirs
#                      with _analyze_zscore.tsv.gz and _reference_zscore.tsv.gz
#                      (i.e. the output dir from beta_to_episcore.py).
#
# Outputs (under <output_base>/global_recall/):
#     _reference_zscore.tsv.gz
#     _analyze_zscore.tsv.gz

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search

BEST_COMBO_CSV=${1:-"${RESULT_BASE}/best_combo_per_chr_global.csv"}
OUTPUT_BASE=${2:-"${RESULT_BASE}/output"}

# Run from the directory containing b2z_for_global_recall.py (where sbatch was called).
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

echo "[$(date)] Starting b2z_for_global_recall on host $(hostname)"
echo "  best-combo CSV : $BEST_COMBO_CSV"
echo "  output base    : $OUTPUT_BASE"
echo "  output dir     : $OUTPUT_BASE/global_recall"
echo "  cpus           : ${SLURM_CPUS_PER_TASK:-4}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 b2z_for_global_recall.py \
        --best-combo-csv "$BEST_COMBO_CSV" \
        --output-base    "$OUTPUT_BASE"

echo "[$(date)] Finished b2z_for_global_recall"
