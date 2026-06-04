#!/usr/bin/bash
#SBATCH -J b2z_enlarged_ref
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=128
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Sweep over enlarged reference sizes and report per-run MCC.
#
# Modes (--mode, passed through to b2z_for_enlarged_ref.py):
#   all (default)
#       For each N from --n-min up to the dev+Normal pool size, up to --runs
#       unique random reference sub-sets; single MCC on test + dev non-Normal.
#   isolated
#       100 random 50-sample reference pools from dev Normal; per pool, 50
#       unique reference draws per N in [5, 45]; separate MCC_dev and MCC_test.
#       Also writes top_100_reference_list.csv for N in [35, 45].
#
# Usage:
#     sbatch run_b2z_for_enlarged_ref.slurm \
#            [<combo_dir>] [<meta_csv>] [<output_dir>] [<mode>]
#
# Defaults match the 20260508-grid_search run:
#     combo_dir   : <RESULT_BASE>/output/threshold_0.5_recall_0.65
#     meta_csv    : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv
#     output_dir  : <RESULT_BASE>/output/enlarged_reference  (all mode)
#                   <RESULT_BASE>/output/enlarged_reference_isolated  (isolated)
#     mode        : all

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
META_CSV_DEFAULT=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv

COMBO_DIR=${1:-"${RESULT_BASE}/output/threshold_0.5_recall_0.65"}
META_CSV=${2:-"$META_CSV_DEFAULT"}
MODE=${4:-all}

case "$MODE" in
    all)
        OUTPUT_DIR=${3:-"${RESULT_BASE}/output/enlarged_reference"}
        ;;
    isolated)
        OUTPUT_DIR=${3:-"${RESULT_BASE}/output/enlarged_reference_isolated"}
        ;;
    *)
        echo "ERROR: mode must be 'all' or 'isolated', got: $MODE" >&2
        exit 1
        ;;
esac

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
mkdir -p "$(dirname "$OUTPUT_DIR")"

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
echo "  output dir     : $OUTPUT_DIR"
echo "  mode           : $MODE"
echo "  cpus / threads : ${THREADS}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 b2z_for_enlarged_ref.py \
        --combo-dir   "$COMBO_DIR" \
        --meta-csv    "$META_CSV" \
        --output-dir  "$OUTPUT_DIR" \
        --mode        "$MODE" \
        --threads     "$THREADS"

echo "[$(date)] Finished b2z_for_enlarged_ref"
