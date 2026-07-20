#!/usr/bin/bash
#SBATCH -J patch_zscore_pq
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/patch_zscore_parquet_%j.log
#SBATCH -e logs/patch_zscore_parquet_%j.log

set -euo pipefail

MISSING_TSV=${1:?missing_coverage.tsv required}
PARQUET_PATH=${2:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng/zscore_grid_search.parquet}
ZSCORE_ROOT=${3:-/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search/zscore_results}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] Patching zscore parquet"
echo "  missing tsv : $MISSING_TSV"
echo "  parquet     : $PARQUET_PATH"
echo "  zscore root : $ZSCORE_ROOT"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 patch_zscore_grid_parquet.py \
        --missing-tsv "$MISSING_TSV" \
        --zscore-root "$ZSCORE_ROOT" \
        --parquet-path "$PARQUET_PATH"

echo "[$(date)] Re-checking coverage"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 check_grid_coverage.py \
        --input-dir "$INPUT_DIR" \
        --output-dir "$(dirname "$MISSING_TSV")" \
        --fail-on-missing

echo "[$(date)] Patch + coverage check complete"
