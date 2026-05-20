#!/usr/bin/bash
#SBATCH -J b2z_enlarged_ref
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=128
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Sweep over enlarged reference sizes and report per-run MCC.
#
# For each N from --n-min (default 10) up to the size of the dev+Normal pool,
# up to --runs *unique* random sub-sets of N reference candidates are drawn
# (default 10000 per N). For every draw the script:
#   1. recomputes per-chromosome reference statistics from the chosen N
#      samples;
#   2. recomputes s_inter for a *fixed* evaluation set (test set + dev samples
#      with a non-Normal label) using the new reference statistics;
#   3. records TP/TN/FP/FN/MCC in report.csv.
#
# No per-run _reference_zscore.tsv.gz / _analyze_zscore.tsv.gz are produced;
# only <output_base>/<output_subdir>/report.csv is written.
#
# Usage:
#     sbatch run_b2z_for_enlarged_ref.slurm \
#            [<combo_dir>] [<meta_csv>] [<output_base>] [<output_subdir>]
#
# Defaults match the 20260508-grid_search run:
#     combo_dir      : <RESULT_BASE>/output/threshold_0.5_recall_0.65
#     meta_csv       : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv
#     output_base    : <RESULT_BASE>/output
#     output_subdir  : enlarged_reference

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
META_CSV_DEFAULT=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv

COMBO_DIR=${1:-"${RESULT_BASE}/output/threshold_0.5_recall_0.65"}
META_CSV=${2:-"$META_CSV_DEFAULT"}
OUTPUT_BASE=${3:-"${RESULT_BASE}"}
OUTPUT_SUBDIR=${4:-"enlarged_reference"}

# Run from the directory containing b2z_for_enlarged_ref.py (where sbatch was called).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

if [ ! -d "$COMBO_DIR" ]; then
    echo "ERROR: combo dir not found: $COMBO_DIR" >&2
    exit 1
fi
if [ ! -f "$META_CSV" ]; then
    echo "ERROR: meta CSV not found: $META_CSV" >&2
    exit 1
fi
mkdir -p "$OUTPUT_BASE"

THREADS=${SLURM_CPUS_PER_TASK:-8}

# Keep BLAS/OpenMP from oversubscribing while our own ThreadPoolExecutor
# manages parallelism; each NumPy call should run single-threaded inside a
# worker thread.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] Starting b2z_for_enlarged_ref on host $(hostname)"
echo "  combo dir      : $COMBO_DIR"
echo "  meta csv       : $META_CSV"
echo "  output base    : $OUTPUT_BASE"
echo "  output subdir  : $OUTPUT_SUBDIR"
echo "  output dir     : $OUTPUT_BASE/$OUTPUT_SUBDIR"
echo "  cpus / threads : ${THREADS}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 b2z_for_enlarged_ref.py \
        --combo-dir     "$COMBO_DIR" \
        --meta-csv      "$META_CSV" \
        --output-base   "$OUTPUT_BASE" \
        --output-subdir "$OUTPUT_SUBDIR" \
        --threads       "$THREADS"

echo "[$(date)] Finished b2z_for_enlarged_ref"
