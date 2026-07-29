#!/usr/bin/bash
#SBATCH -J val_b2z
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH -o logs/val_b2z_%a.log
#SBATCH -e logs/val_b2z_%a.log

set -euo pipefail

# Array over (threshold, recall) pairs listed in tasks file:
#   lines: threshold<TAB>recall

TASKS=${1:?tasks file threshold\\trecall}
ANALYZE_DIR=${2:?analyze beta dir}
OUT_ROOT=${3:?output root for threshold_X_recall_Y}
REF_DIR=${REF_DIR:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search/threshold_0.5/early_ref_n_17}
CPG_DIR=${CPG_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list}
B2Z=${B2Z:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/scripts/grid_search/beta_to_episcore.py}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

mapfile -t ROWS < <(awk -F'\t' 'NF>=2 {print $1"\t"$2}' "$TASKS")
TASK_ID=${SLURM_ARRAY_TASK_ID:?}
IFS=$'\t' read -r THR REC <<< "${ROWS[$TASK_ID]}"
THR_S=$(python3 -c "print(f'{float(\"$THR\"):g}')")
REC_S=$(python3 -c "print(f'{float(\"$REC\"):g}')")

# Use matching reference dir when available
CAND_REF="/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search/threshold_${THR_S}/early_ref_n_17"
if [[ -d "$CAND_REF" ]]; then
    REF_DIR="$CAND_REF"
fi

CPG_LIST="${CPG_DIR}/220k_cpg_recall_${REC_S}.txt"
OUT_DIR="${OUT_ROOT}/threshold_${THR_S}_recall_${REC_S}"
mkdir -p "$OUT_DIR"

if [[ -f "${OUT_DIR}/_analyze_zscore.tsv.gz" ]]; then
    echo "[$(date)] already done $OUT_DIR"
    exit 0
fi
[[ -f "$CPG_LIST" ]] || { echo "missing $CPG_LIST" >&2; exit 1; }
[[ -d "$ANALYZE_DIR" ]] || { echo "missing $ANALYZE_DIR" >&2; exit 1; }
[[ -d "$REF_DIR" ]] || { echo "missing $REF_DIR" >&2; exit 1; }

echo "[$(date)] val_b2z thr=$THR_S rec=$REC_S"
echo "  analyze=$ANALYZE_DIR"
echo "  ref=$REF_DIR"
echo "  out=$OUT_DIR"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 "$B2Z" \
        --reference-beta-value-dir "$REF_DIR" \
        --analyze-beta-value-dir "$ANALYZE_DIR" \
        --output-prefix "${OUT_DIR}/" \
        --depth 30 \
        --cpg-list "$CPG_LIST" \
        --chr-list 1-22 \
        --ncpus "${SLURM_CPUS_PER_TASK:-16}"

echo "[$(date)] OK $OUT_DIR"
