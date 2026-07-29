#!/usr/bin/bash
#SBATCH -J val_zscore
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH -o logs/val_zscore_%a.log
#SBATCH -e logs/val_zscore_%a.log

set -euo pipefail
export PYTHONUNBUFFERED=1

# Args:
#   1) zscore_tasks.tsv (score_type/sample/threshold/recall)
#   2) work dir with samplesheet.csv; writes zscore_results/ under here
#   3) optional: shared zscore root for beds/refs (default 20260513-grid_search)

TASKS_TSV=${1:?tasks tsv required}
WORK_DIR=${2:?work dir required}
ZSCORE_ROOT=${3:-/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search}
REF_NAME=${REF_NAME:-igtc_220k.early17}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
ANALYZER=${ANALYZER:-/lustre1/cqyi/AIPT_2.0/workflow/zscore/mq_zscore_analyzer.py}
SAMPLESHEET="${WORK_DIR}/samplesheet.csv"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

TASK_OFFSET=${TASK_OFFSET:-0}
TASK_ID=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))

mapfile -t ROWS < <(awk -F'\t' 'NR>1 {print $2"\t"$3"\t"$4}' "$TASKS_TSV")
N_TASKS=${#ROWS[@]}
if (( TASK_ID < 0 || TASK_ID >= N_TASKS )); then
    echo "[$(date)] array=${SLURM_ARRAY_TASK_ID} offset=${TASK_OFFSET} -> ${TASK_ID}: out of range (n=${N_TASKS})"
    exit 0
fi

IFS=$'\t' read -r SAMPLE THRESHOLD RECALL <<< "${ROWS[$TASK_ID]}"
THR_S=$(python3 -c "print(f'{float(\"$THRESHOLD\"):g}')")
REC_S=$(python3 -c "print(f'{float(\"$RECALL\"):g}')")

OUTDIR="${WORK_DIR}/zscore_results/recall.${REC_S}_cutoff.${THR_S}"
OUT_CSV="${OUTDIR}/${SAMPLE}.${THR_S}.1.0.220k_cpg_recall_${REC_S}.NoLen.zscore.csv"
RECALL_BED="${ZSCORE_ROOT}/recall_list_bed/220k_cpg_recall_${REC_S}.bed"
REF_CSV="${ZSCORE_ROOT}/reference/${REF_NAME}/${REF_NAME}.ref_list.${THR_S}.1.0.220k_cpg_recall_${REC_S}.NoLen.reference.csv"

if [[ -f "$OUT_CSV" ]] && [[ "$(wc -l < "$OUT_CSV")" -ge 23 ]]; then
    echo "[$(date)] Already complete: $OUT_CSV"
    exit 0
fi
[[ -f "$RECALL_BED" ]] || { echo "missing bed $RECALL_BED" >&2; exit 1; }
[[ -f "$REF_CSV" ]] || { echo "missing ref $REF_CSV" >&2; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "missing $SAMPLESHEET" >&2; exit 1; }

MQ_FILES=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $3}' "$SAMPLESHEET" | paste -sd',' -)
[[ -n "$MQ_FILES" ]] || { echo "no MQ for $SAMPLE" >&2; exit 1; }

mkdir -p "$OUTDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "[$(date)] val_zscore array=${TASK_ID} sample=$SAMPLE thr=$THR_S rec=$REC_S"

# Run from WORK_DIR so analyzer relative paths land under WORK_DIR/zscore_results
cd "$WORK_DIR"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 "$ANALYZER" cal_zscore \
        --reference "$REF_CSV" \
        --sample "$SAMPLE" \
        --cutoff "$THR_S" \
        --mTcount 1 \
        --recall "$RECALL_BED" \
        --mq_files "$MQ_FILES" \
        --gender unknown \
        --out_dir "$OUTDIR" \
        --threads 1

[[ -f "$OUT_CSV" ]] && [[ "$(wc -l < "$OUT_CSV")" -ge 23 ]] || {
    echo "ERROR incomplete $OUT_CSV" >&2
    exit 1
}
echo "[$(date)] OK $OUT_CSV"
