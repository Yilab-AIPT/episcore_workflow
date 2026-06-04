#!/usr/bin/bash
#SBATCH -J b2z
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run beta_to_zscore.py for one (threshold, recall) combination.
#
# Usage:
#     sbatch --job-name=b2z_t<thres>_r<recall> \
#            run_beta_to_zscore.slurm <threshold> <recall>
#
# Inputs:
#     reference dir : <RESULT_BASE>/threshold_<thres>/early_ref_n_17
#     analyze dir   : <RESULT_BASE>/threshold_<thres>/analyze_n_283
#     cpg list      : <CPG_DIR>/220k_cpg_recall_<recall>.txt
#
# Outputs:
#     <RESULT_BASE>/output/threshold_<thres>_recall_<recall>/
#         _reference_zscore.tsv.gz
#         _analyze_zscore.tsv.gz
#
# The job name is also used to name the log file (logs/%x.log), so be sure to
# pass --job-name when submitting (submit.sh handles that automatically).

set -euo pipefail

threshold=${1:?Usage: sbatch run_beta_to_zscore.slurm <threshold> <recall>}
recall=${2:?Usage: sbatch run_beta_to_zscore.slurm <threshold> <recall>}

# Run from the directory containing beta_to_zscore.py (where sbatch was called).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260525-grid_search_240k_panel_240k_model
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_220k

REF_DIR="${RESULT_BASE}/threshold_${threshold}/220k_early_ref_n_16"
ANALYZE_DIR="${RESULT_BASE}/threshold_${threshold}/220k_analyze_n_26"
CPG_LIST="${CPG_DIR}/220k_cpg_recall_${recall}.txt"
OUT_DIR="${RESULT_BASE}/output_220k/threshold_${threshold}_recall_${recall}"

for path in "$REF_DIR" "$ANALYZE_DIR"; do
    if [ ! -d "$path" ]; then
        echo "ERROR: directory not found: $path" >&2
        exit 1
    fi
done
if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: CpG list not found: $CPG_LIST" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "[$(date)] Starting threshold=${threshold} recall=${recall} on host $(hostname)"
echo "  reference dir : $REF_DIR"
echo "  analyze dir   : $ANALYZE_DIR"
echo "  cpg list      : $CPG_LIST"
echo "  output dir    : $OUT_DIR"
echo "  cpus          : ${SLURM_CPUS_PER_TASK:-32}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 beta_to_zscore.py \
        --reference-beta-value-dir "$REF_DIR" \
        --analyze-beta-value-dir   "$ANALYZE_DIR" \
        --output-prefix            "${OUT_DIR}/" \
        --depth 30 \
        --cpg-list                 "$CPG_LIST" \
        --chr-list 1-22 \
        --ncpus "${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] Finished threshold=${threshold} recall=${recall}"
