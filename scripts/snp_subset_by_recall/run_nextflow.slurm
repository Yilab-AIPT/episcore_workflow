#!/usr/bin/bash
#SBATCH -J snp_nf
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run episcore SNP FF workflow (est_ff_from_pileup) for one recall-level SNP panel.
#
# Usage:
#     sbatch --job-name=nf_r<recall> run_nextflow.slurm <recall>
#
# Example:
#     sbatch --job-name=nf_r0.05 run_nextflow.slurm 0.05

set -euo pipefail

recall=${1:?Usage: sbatch run_nextflow.slurm <recall>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

nextflow=/appsnew/home/syfan/softwares/nextflow/nextflow
NF_MAIN=/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf
SAMPLESHEET=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260604-snp_subset/pileup_samplesheet.csv
META_BASE=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260604-snp_subset
RESULTS_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260604-snp_subset
WORK_BASE=/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir

SNP_LIST="${META_BASE}/snp_subset/snp_for_recall_${recall}.tsv"
OUTDIR="${RESULTS_BASE}/snp_for_recall_${recall}"
WORKDIR="${WORK_BASE}/snp_recall_${recall}"

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi
if [ ! -f "$SNP_LIST" ]; then
    echo "ERROR: SNP panel not found: $SNP_LIST" >&2
    exit 1
fi

mkdir -p "$OUTDIR" "$WORKDIR"

echo "[$(date)] Starting nextflow recall=${recall} on host $(hostname)"
echo "  samplesheet : $SAMPLESHEET"
echo "  snp_list    : $SNP_LIST"
echo "  outdir      : $OUTDIR"
echo "  work-dir    : $WORKDIR"

"$nextflow" run "$NF_MAIN" \
    -profile early,alioth_slurm \
    -work-dir "$WORKDIR" \
    --step est_ff_from_pileup \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --snp_list "$SNP_LIST" \
    --ff_precision 0.001 \
    2>&1 | tee "${OUTDIR}/nextflow_run.log"

echo "[$(date)] Finished recall=${recall}"
