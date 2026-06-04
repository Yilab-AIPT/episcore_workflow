#!/usr/bin/bash
#SBATCH -J xo_b2z
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run run_xo_episcore.py for one recall level.
#
# Usage:
#     sbatch --job-name=xo_b2z_r<recall> run_xo_episcore.slurm <recall>
#
# Inputs:
#     samples meta : <META_DIR>/samples_meta.csv
#     cpg list     : <CPG_DIR>/240k_cpg_recall_<recall>.txt
#
# Outputs:
#     <RESULT_BASE>/recall_<recall>/
#         _reference_zscore.tsv.gz
#         _analyze_zscore.tsv.gz

set -euo pipefail

recall=${1:?Usage: sbatch run_xo_episcore.slurm <recall>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

META_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260526-grid_search_XO
RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260526-grid_search_XO
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_240k

SAMPLES_META="${META_DIR}/samples_meta.csv"
CPG_LIST="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
OUT_DIR="${RESULT_BASE}/recall_${recall}"

if [ ! -f "$SAMPLES_META" ]; then
    echo "ERROR: samples meta not found: $SAMPLES_META" >&2
    exit 1
fi
if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: CpG list not found: $CPG_LIST" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "[$(date)] Starting recall=${recall} on host $(hostname)"
echo "  samples meta  : $SAMPLES_META"
echo "  cpg list      : $CPG_LIST"
echo "  output dir    : $OUT_DIR"
echo "  cpus          : ${SLURM_CPUS_PER_TASK:-32}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 run_xo_episcore.py \
        --samples-meta "$SAMPLES_META" \
        --cpg-list       "$CPG_LIST" \
        --output-prefix  "${OUT_DIR}/" \
        --chr-list       1-22,X \
        --ncpus          "${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] Finished recall=${recall}"
