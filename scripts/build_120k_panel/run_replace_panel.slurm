#!/usr/bin/bash
#SBATCH -J cpg_replace
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o logs/cpg_replace.log
#SBATCH -e logs/cpg_replace.log

# Merge per-sample depths and build replaced recall-0.65 site list.
#
# Usage:
#     sbatch run_replace_panel.slurm
#
# Typically submitted with dependency after depth array completes.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

OUT_DIR=${OUT_DIR:-/lustre1/cqyi/AIPT_2.0/results/small_panel}
RECALL_DIR=${RECALL_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list}
SAMPLESHEET=${SAMPLESHEET:-/lustre1/cqyi/syfan/nipt_article_plot/dev_and_test_mqres_samplesheet.csv}
CONTAINER=${CONTAINER:-/lustre1/cqyi/syfan/images/common_tools.sif}
N_SAMPLES=${N_SAMPLES:-50}
SEED=${SEED:-42}

mkdir -p logs "$OUT_DIR"

echo "[$(date)] merge + replace on host $(hostname)"
echo "  out_dir : $OUT_DIR"
echo "  n_samples: $N_SAMPLES seed: $SEED"

singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" \
    python3 "${SCRIPT_DIR}/replace_deeper_sites.py" \
        --samplesheet "$SAMPLESHEET" \
        --recall-dir "$RECALL_DIR" \
        --n-samples "$N_SAMPLES" \
        --seed "$SEED" \
        --out-dir "$OUT_DIR" \
        --selected-samples "${OUT_DIR}/selected_samples.csv" \
        --depth-dir "${OUT_DIR}/depth_per_sample"

echo "[$(date)] Finished replace panel build"
