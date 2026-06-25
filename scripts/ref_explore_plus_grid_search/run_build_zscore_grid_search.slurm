#!/usr/bin/bash
#SBATCH -J build_zscore_grid
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32

# Parse zscore grid-search CSVs into one long-format parquet table.
#
# Usage:
#     sbatch run_build_zscore_grid_search.slurm [<zscore_grid_search_dir>] [<output_dir>]
#
# Defaults:
#     zscore_grid_search_dir : /lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search/zscore_results
#     output_dir             : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng

set -euo pipefail

ZSCORE_GRID_SEARCH_DIR=${1:-/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search/zscore_results}
OUTPUT_DIR=${2:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs "$OUTPUT_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] Starting build_zscore_grid_search on host $(hostname)"
echo "  zscore grid dir : $ZSCORE_GRID_SEARCH_DIR"
echo "  output dir      : $OUTPUT_DIR"
echo "  cpus            : ${SLURM_CPUS_PER_TASK:-32}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 build_zscore_grid_search.py \
        --zscore-grid-search-dir "$ZSCORE_GRID_SEARCH_DIR" \
        --output-dir "$OUTPUT_DIR" \
        --threads "${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] Finished build_zscore_grid_search"
