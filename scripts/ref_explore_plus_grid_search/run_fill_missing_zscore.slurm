#!/usr/bin/bash
#SBATCH -J fill_zscore
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH -o logs/fill_missing_zscore_%a.log
#SBATCH -e logs/fill_missing_zscore_%a.log

set -euo pipefail
export PYTHONUNBUFFERED=1

# Args:
#   1) missing_coverage.tsv (score_type/sample/threshold/recall)
#   2) zscore grid-search root (contains samplesheet.csv, recall_list_bed/, reference/, zscore_results/)

MISSING_TSV=${1:?missing_coverage.tsv required}
ZSCORE_ROOT=${2:-/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search}
REF_NAME=${REF_NAME:-igtc_220k.early17}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
ANALYZER=${ANALYZER:-/lustre1/cqyi/AIPT_2.0/workflow/zscore/mq_zscore_analyzer.py}
SAMPLESHEET=${SAMPLESHEET:-${ZSCORE_ROOT}/samplesheet.csv}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}

# Select zscore rows only; 0-based array index
mapfile -t ROWS < <(
    awk -F'\t' 'NR>1 && $1=="zscore" {print $2"\t"$3"\t"$4}' "$MISSING_TSV"
)
N_TASKS=${#ROWS[@]}
if (( TASK_ID < 0 || TASK_ID >= N_TASKS )); then
    echo "[$(date)] array=${TASK_ID}: out of range (n=${N_TASKS}), nothing to do"
    exit 0
fi

IFS=$'\t' read -r SAMPLE THRESHOLD RECALL <<< "${ROWS[$TASK_ID]}"
THR_S=$(python3 -c "print(f'{float(\"$THRESHOLD\"):g}')")
REC_S=$(python3 -c "print(f'{float(\"$RECALL\"):g}')")

OUTDIR="${ZSCORE_ROOT}/zscore_results/recall.${REC_S}_cutoff.${THR_S}"
OUT_CSV="${OUTDIR}/${SAMPLE}.${THR_S}.1.0.220k_cpg_recall_${REC_S}.NoLen.zscore.csv"
RECALL_BED="${ZSCORE_ROOT}/recall_list_bed/220k_cpg_recall_${REC_S}.bed"
REF_CSV="${ZSCORE_ROOT}/reference/${REF_NAME}/${REF_NAME}.ref_list.${THR_S}.1.0.220k_cpg_recall_${REC_S}.NoLen.reference.csv"

if [[ -f "$OUT_CSV" ]] && [[ "$(wc -l < "$OUT_CSV")" -ge 23 ]]; then
    echo "[$(date)] Already complete: $OUT_CSV"
    exit 0
fi

if [[ ! -f "$RECALL_BED" ]]; then
    echo "ERROR: missing recall bed: $RECALL_BED" >&2
    exit 1
fi
if [[ ! -f "$REF_CSV" ]]; then
    echo "ERROR: missing reference: $REF_CSV" >&2
    exit 1
fi
if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "ERROR: missing samplesheet: $SAMPLESHEET" >&2
    exit 1
fi

MQ_FILES=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $3}' "$SAMPLESHEET" | paste -sd',' -)
if [[ -z "$MQ_FILES" ]]; then
    echo "ERROR: no deconv_res rows for sample=$SAMPLE in $SAMPLESHEET" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] fill_missing_zscore array=${TASK_ID}/${N_TASKS}"
echo "  sample     : $SAMPLE"
echo "  threshold  : $THR_S"
echo "  recall     : $REC_S"
echo "  out csv    : $OUT_CSV"
echo "  mq files   : $MQ_FILES"

# Analyzer writes relative to --out_dir cwd; run from ZSCORE_ROOT so paths match historical layout
cd "$ZSCORE_ROOT"

# Use container PATH as-is (bedtools lives in /tools/bedtools_*/bin inside common_tools.sif).
# Do not wrap with bash -l — login shells can wipe the container PATH via host bashrc.
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

if [[ ! -f "$OUT_CSV" ]] || [[ "$(wc -l < "$OUT_CSV")" -lt 23 ]]; then
    echo "ERROR: expected output missing or incomplete: $OUT_CSV" >&2
    exit 1
fi

echo "[$(date)] Finished OK: $OUT_CSV"
