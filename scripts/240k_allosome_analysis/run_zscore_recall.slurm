#!/usr/bin/bash
#SBATCH -J allo_z
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Read-count zscore (chrX) for one recall; prob_class_1 cutoff=0.85.

set -euo pipefail

recall=${1:?Usage: sbatch run_zscore_recall.slurm <recall>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples
OUTPUT_DIR=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_240k
SIF=/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif

SAMPLES_META="${INPUT_DIR}/zscore_samples_meta.csv"
CPG_LIST="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
OUT_DIR="${OUTPUT_DIR}/zscore_recall/recall_${recall}"

if [ ! -f "$SAMPLES_META" ]; then
    echo "ERROR: missing $SAMPLES_META" >&2
    exit 1
fi
if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: missing $CPG_LIST" >&2
    exit 1
fi

# Skip only if output exists and chrX_percentage is non-zero (legacy bug: 23→chr23).
if [ -f "${OUT_DIR}/_analyze_zscore.tsv.gz" ]; then
    if singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3 - "${OUT_DIR}/_analyze_zscore.tsv.gz" <<'PY'
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t")
raise SystemExit(0 if "chrX_percentage" in df.columns and float(df["chrX_percentage"].max()) > 0 else 1)
PY
    then
        echo "[$(date)] recall=${recall} already valid, skipping"
        exit 0
    fi
    echo "[$(date)] recall=${recall} invalid (chrX all-zero); regenerating"
    rm -f "${OUT_DIR}/_analyze_zscore.tsv.gz" "${OUT_DIR}/_reference_percentage.tsv.gz"
fi

mkdir -p "$OUT_DIR"
echo "[$(date)] zscore recall=${recall} cutoff=0.85 host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-4}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    "$SIF" \
    python3 run_zscore_recall.py \
        --samples-meta "$SAMPLES_META" \
        --cpg-list "$CPG_LIST" \
        --output-prefix "${OUT_DIR}/" \
        --cutoff 0.85 \
        --ncpus "${SLURM_CPUS_PER_TASK:-4}"

echo "[$(date)] Finished zscore recall=${recall}"
