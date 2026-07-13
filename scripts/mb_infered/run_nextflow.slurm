#!/usr/bin/bash
#SBATCH -J mb_infered_nf
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run aipt_ref_40 scoring for methylbert-inferred deconv results (MB panel).
#
# Uses a fixed-combo grid-search asset dir built from the MB episcore reference
# matrix and MB zscore grid-search outputs (no ezscore assets).
#
# Usage:
#     sbatch run_nextflow.slurm
#
# Optional environment overrides:
#     OUTDIR, WORKDIR, SAMPLESHEET, GRID_SEARCH_RESULT, NF_RESUME (default 1)

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

nextflow=/appsnew/home/syfan/softwares/nextflow/nextflow
NF_MAIN=/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf
PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore
COMMON_TOOLS="${PROJECT_DIR}/containers/common_tools.sif"

SAMPLESHEET=${SAMPLESHEET:-/lustre1/cqyi/syfan/nf_methylbert_infer/data/220k_data_mb_vs_mq/20260712_inferred_summarized_samplesheet.csv}
GRID_SEARCH_RESULT=${GRID_SEARCH_RESULT:-/lustre1/cqyi/syfan/nf_methylbert_infer/data/220k_data_mb_vs_mq/20260712_grid_search_fixed_combo}
OUTDIR=${OUTDIR:-/lustre1/cqyi/syfan/nf_methylbert_infer/results/20260712-mb_infered}
WORKDIR=${WORKDIR:-/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir/mb_infered_20260712}
NF_RESUME=${NF_RESUME:-1}

EPISCORE_REF_MATRIX="/lustre1/cqyi/syfan/nf_methylbert_infer/results/old_model_downsampled_episcore/output/threshold_0.5_recall_0.66/_reference_zscore.tsv.gz"
ZSCORE_REF_DIR="/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260611-MB-grid_search/zscore_results/recall.0.92_cutoff.0.55"
EPISCORE_THRESHOLD=0.5
EPISCORE_RECALL=0.66
ZSCORE_THRESHOLD=0.55
ZSCORE_RECALL=0.92

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi
for f in "$EPISCORE_REF_MATRIX" "${PROJECT_DIR}/assets/hg38_only_chromsomes.fa" "${PROJECT_DIR}/assets/snps_in_220kprobes_v2.tsv"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: required asset not found: $f" >&2
        exit 1
    fi
done
if [ ! -d "$ZSCORE_REF_DIR" ]; then
    echo "ERROR: zscore reference dir not found: $ZSCORE_REF_DIR" >&2
    exit 1
fi

# Build grid-search assets when missing (episcore/zscore refs only; no ezscore).
if [ ! -f "${GRID_SEARCH_RESULT}/best_reference_matrix.tsv" ]; then
    echo "[$(date)] Building grid-search assets in ${GRID_SEARCH_RESULT}"
    mkdir -p "$GRID_SEARCH_RESULT"
    singularity exec -B /lustre1,/lustre2 "$COMMON_TOOLS" \
        python3 "${SCRIPT_DIR}/build_grid_search_assets.py" \
        --episcore-reference-matrix "$EPISCORE_REF_MATRIX" \
        --zscore-reference-dir "$ZSCORE_REF_DIR" \
        --output-dir "$GRID_SEARCH_RESULT" \
        --episcore-threshold "$EPISCORE_THRESHOLD" \
        --episcore-recall "$EPISCORE_RECALL" \
        --zscore-threshold "$ZSCORE_THRESHOLD" \
        --zscore-recall "$ZSCORE_RECALL"
fi

for f in \
    "${GRID_SEARCH_RESULT}/best_combo_episcore.csv" \
    "${GRID_SEARCH_RESULT}/best_combo_zscore.csv" \
    "${GRID_SEARCH_RESULT}/best_reference_matrix.tsv" \
    "${GRID_SEARCH_RESULT}/best_reference_samples.txt"
do
    if [ ! -f "$f" ]; then
        echo "ERROR: grid-search asset not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR" "$WORKDIR"

echo "[$(date)] Starting aipt_ref_40 (MB infer) on host $(hostname)"
echo "  samplesheet         : $SAMPLESHEET"
echo "  grid_search_result  : $GRID_SEARCH_RESULT"
echo "  outdir              : $OUTDIR"
echo "  work-dir            : $WORKDIR"
echo "  resume              : $NF_RESUME"

RESUME_FLAG=()
if [ "$NF_RESUME" = "1" ]; then
    RESUME_FLAG=(-resume)
fi

"$nextflow" run "$NF_MAIN" \
    -profile aipt_ref_40,alioth_slurm,singularity \
    -work-dir "$WORKDIR" \
    "${RESUME_FLAG[@]}" \
    --step aipt_ref_40 \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --grid_search_result "$GRID_SEARCH_RESULT" \
    2>&1 | tee "${OUTDIR}/nextflow_run.log"

echo "[$(date)] Finished aipt_ref_40 (MB infer)"
