#!/usr/bin/bash
#SBATCH -J snp_excl
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Build exclude-overlap SNP subset for one CpG recall level.
#
# Usage:
#     sbatch --job-name=snp_excl_r<recall> run_build_exclude_snp_subset.slurm <recall>
#
# Example:
#     sbatch --job-name=snp_excl_r0.65 run_build_exclude_snp_subset.slurm 0.65

set -euo pipefail

recall=${1:?Usage: sbatch run_build_exclude_snp_subset.slurm <recall>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

RECALL_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list
OUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260604-snp_subset/snp_exclude_subset

echo "[$(date)] Starting exclude recall=${recall} on host $(hostname)"
echo "  recall dir : $RECALL_DIR"
echo "  output dir : $OUT_DIR"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 build_exclude_snp_subset.py \
        --recall "$recall" \
        --recall-dir "$RECALL_DIR" \
        --out-dir "$OUT_DIR"

echo "[$(date)] Finished exclude recall=${recall}"
