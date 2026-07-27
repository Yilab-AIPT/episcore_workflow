#!/usr/bin/bash
#SBATCH -J allo_nf
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log
#SBATCH --time=72:00:00

# Run Nextflow split_bam → episcore/FF for the 10 new 240k allosome samples.
# BAM-split threshold 0.5; full 240k CpG list so betas support the recall sweep.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples
OUTPUT_DIR=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples
PROJECT=/lustre1/cqyi/AIPT_2.0/workflow/episcore
SAMPLESHEET="${INPUT_DIR}/samplesheet_nf.csv"
CPG_LIST="${PROJECT}/assets/cpgs_in_240k_probes.txt"

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: missing $SAMPLESHEET — run prepare_inputs.py first" >&2
    exit 1
fi
if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: missing full 240k CpG list: $CPG_LIST" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

NEXTFLOW=${NEXTFLOW:-/appsnew/home/syfan/softwares/nextflow/nextflow}
if [ ! -x "$NEXTFLOW" ]; then
    echo "ERROR: nextflow not found: $NEXTFLOW" >&2
    exit 1
fi

echo "[$(date)] Starting Nextflow on $(hostname)"
echo "  samplesheet : $SAMPLESHEET"
echo "  outdir      : $OUTPUT_DIR"
echo "  threshold   : 0.5"
echo "  cpg_list    : $CPG_LIST"
echo "  nextflow    : $NEXTFLOW"

cd "$PROJECT"
"$NEXTFLOW" run "$PROJECT/main.nf" \
    -profile early_240k,alioth_slurm,singularity \
    --step split_bam \
    --input "$SAMPLESHEET" \
    --outdir "$OUTPUT_DIR" \
    --threshold 0.5 \
    --cpg_list "$CPG_LIST" \
    -resume \
    2>&1 | tee "${OUTPUT_DIR}/nextflow_run.log"

echo "[$(date)] Nextflow finished"
