#!/usr/bin/bash
#SBATCH -J cpg_depth
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o logs/cpg_depth_%a.log
#SBATCH -e logs/cpg_depth_%a.log

# Compute per-site CpG depth for one sample (SLURM array task).
#
# Usage:
#     sbatch --array=1-N run_depth_array.slurm
#
# Expects:
#     selected_samples.csv  in OUT_DIR (written by submit.sh step 1)
#     OUT_DIR/depth_per_sample/ for outputs

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

OUT_DIR=${OUT_DIR:-/lustre1/cqyi/AIPT_2.0/results/small_panel}
SITES=${SITES:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list/220k_cpg_recall_0.6.txt}
CONTAINER=${CONTAINER:-/lustre1/cqyi/syfan/images/common_tools.sif}
SELECTED_SAMPLES="${OUT_DIR}/selected_samples.csv"
DEPTH_DIR="${OUT_DIR}/depth_per_sample"

mkdir -p logs "$DEPTH_DIR"

if [ ! -f "$SELECTED_SAMPLES" ]; then
    echo "ERROR: missing $SELECTED_SAMPLES — run submit.sh first" >&2
    exit 1
fi

mapfile -t SAMPLE_IDS < <(tail -n +2 "$SELECTED_SAMPLES" | cut -d, -f1)
mapfile -t BAM_PATHS < <(tail -n +2 "$SELECTED_SAMPLES" | cut -d, -f2)

task_idx=$((SLURM_ARRAY_TASK_ID - 1))
sample_id="${SAMPLE_IDS[$task_idx]}"
bam_path="${BAM_PATHS[$task_idx]}"

echo "[$(date)] array=${SLURM_ARRAY_TASK_ID} sample=${sample_id}"
echo "  BAM: ${bam_path}"
echo "  OUT: ${DEPTH_DIR}"

singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" \
    python3 "${SCRIPT_DIR}/compute_sample_cpg_depth.py" \
        --bam "$bam_path" \
        --sample-id "$sample_id" \
        --sites "$SITES" \
        --output-dir "$DEPTH_DIR"

echo "[$(date)] Finished ${sample_id}"
