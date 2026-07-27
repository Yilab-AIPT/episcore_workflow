#!/usr/bin/bash
#SBATCH -J allo_epi
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Episcore (beta → chrX s_inter) for one recall; BAM-split threshold was 0.5.

set -euo pipefail

recall=${1:?Usage: sbatch run_episcore_recall.slurm <recall>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples
OUTPUT_DIR=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_240k
SIF=/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif

SAMPLES_META="${INPUT_DIR}/episcore_samples_meta.csv"
CPG_LIST="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
OUT_DIR="${OUTPUT_DIR}/episcore_recall/recall_${recall}"

if [ ! -f "$SAMPLES_META" ]; then
    echo "ERROR: missing $SAMPLES_META" >&2
    exit 1
fi
if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: missing $CPG_LIST" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
echo "[$(date)] episcore recall=${recall} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-16}"

# Skip if already complete
if [ -f "${OUT_DIR}/_analyze_zscore.tsv.gz" ]; then
    echo "[$(date)] recall=${recall} already done, skipping"
    exit 0
fi

singularity exec -B /lustre1,/lustre2,/appsnew \
    "$SIF" \
    python3 run_episcore_recall.py \
        --samples-meta "$SAMPLES_META" \
        --cpg-list "$CPG_LIST" \
        --output-prefix "${OUT_DIR}/" \
        --chr-list 1-22,X \
        --ncpus "${SLURM_CPUS_PER_TASK:-16}"

echo "[$(date)] Finished episcore recall=${recall}"
