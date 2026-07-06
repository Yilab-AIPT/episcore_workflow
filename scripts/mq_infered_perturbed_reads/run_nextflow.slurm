#!/usr/bin/bash
#SBATCH -J mq_perturbed_nf
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run episcore perturbed_res workflow for methylqueen-inferred perturbation reads.
#
# conf/perturbed_res.config targets methylbert-inferred results (220k recall_0.77
# panel, threshold 0.67). MQ inputs use the early panel instead (CpG_recall0.65,
# early reference matrix, threshold 0.5).
#
# Usage:
#     sbatch run_nextflow.slurm
#
# Optional environment overrides:
#     OUTDIR, WORKDIR, SAMPLESHEET, NF_RESUME (default 1: resume cached tasks)

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

nextflow=/appsnew/home/syfan/softwares/nextflow/nextflow
NF_MAIN=/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf
PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore

SAMPLESHEET=${SAMPLESHEET:-/lustre1/cqyi/syfan/nf_methylqueen_infer/results/methylation_conversion_experiment/20260703/mq_infered_result.csv}
OUTDIR=${OUTDIR:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260703-mq_infered_perturbed_reads}
WORKDIR=${WORKDIR:-/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir/mq_infered_perturbed_20260703}

# Early-panel params for MQ deconv results (override perturbed_res.config defaults).
CPG_LIST="${PROJECT_DIR}/assets/CpG_recall0.65.txt"
REFERENCE_EPISCORE_MATRIX="${PROJECT_DIR}/assets/early_reference_beta_zscore.tsv"
THRESHOLD=0.5
FF_PRECISION=0.0001
NF_RESUME=${NF_RESUME:-1}

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi
for f in "$CPG_LIST" "$REFERENCE_EPISCORE_MATRIX" "${PROJECT_DIR}/assets/hg38_only_chromsomes.fa" "${PROJECT_DIR}/assets/snps_in_220kprobes_v2.tsv"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: required asset not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR" "$WORKDIR"

echo "[$(date)] Starting perturbed_res (MQ infer) on host $(hostname)"
echo "  samplesheet : $SAMPLESHEET"
echo "  outdir      : $OUTDIR"
echo "  work-dir    : $WORKDIR"
echo "  cpg_list    : $CPG_LIST"
echo "  ref_matrix  : $REFERENCE_EPISCORE_MATRIX"
echo "  threshold   : $THRESHOLD"
echo "  ff_precision: $FF_PRECISION"
echo "  resume      : $NF_RESUME"

RESUME_FLAG=()
if [ "$NF_RESUME" = "1" ]; then
    RESUME_FLAG=(-resume)
fi

"$nextflow" run "$NF_MAIN" \
    -profile early,alioth_slurm,singularity \
    -work-dir "$WORKDIR" \
    "${RESUME_FLAG[@]}" \
    --step perturbed_res \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --cpg_list "$CPG_LIST" \
    --reference_episcore_matrix "$REFERENCE_EPISCORE_MATRIX" \
    --threshold "$THRESHOLD" \
    --ff_precision "$FF_PRECISION" \
    2>&1 | tee "${OUTDIR}/nextflow_run.log"

echo "[$(date)] Finished perturbed_res (MQ infer)"
