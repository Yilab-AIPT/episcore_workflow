#!/usr/bin/bash
#SBATCH -J est_ff_bam
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o %x_%j.log
#SBATCH -e %x_%j.log

# Run episcore SNP FF workflow starting from clean BAM + deconv results.
#
# Usage:
#     sbatch run_est_ff_from_bam.slurm <samplesheet> <outdir>
#
# Samplesheet columns: sample, clean_bam, deconv_res
#
# Optional environment overrides:
#     NF_PROFILE   Nextflow profile (default: early,alioth_slurm,singularity)
#     FF_PRECISION FF search precision (default: 0.0001, from early.config)

set -euo pipefail

SAMPLESHEET=${1:?Usage: sbatch run_est_ff_from_bam.slurm <samplesheet> <outdir>}
OUTDIR=${2:?Usage: sbatch run_est_ff_from_bam.slurm <samplesheet> <outdir>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

NEXTFLOW="${NEXTFLOW:-/appsnew/usr/Nextflow/nextflow-24.10.3/nextflow}"
NF_MAIN="${NF_MAIN:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf}"
WORK_BASE="${WORK_BASE:-/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir}"
NF_PROFILE="${NF_PROFILE:-early,alioth_slurm,singularity}"
FF_PRECISION="${FF_PRECISION:-0.0001}"

OUTDIR="$(readlink -f "$OUTDIR")"
WORKDIR="${WORK_BASE}/est_ff_bam_$(basename "$OUTDIR")"

mkdir -p "$OUTDIR" "$WORKDIR"

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi
if [ ! -x "$NEXTFLOW" ] && ! command -v "$NEXTFLOW" >/dev/null 2>&1; then
    echo "ERROR: nextflow not found: $NEXTFLOW" >&2
    exit 1
fi

echo "[$(date)] Starting est_ff_from_bam on host $(hostname)"
echo "  samplesheet : $SAMPLESHEET"
echo "  outdir      : $OUTDIR"
echo "  work-dir    : $WORKDIR"
echo "  profile     : $NF_PROFILE"
echo "  ff_precision: $FF_PRECISION"

"$NEXTFLOW" run "$NF_MAIN" \
    -profile "$NF_PROFILE" \
    -work-dir "$WORKDIR" \
    --step est_ff_from_bam \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --ff_precision "$FF_PRECISION" \
    2>&1 | tee "${OUTDIR}/nextflow_run.log"

echo "[$(date)] Finished est_ff_from_bam"
