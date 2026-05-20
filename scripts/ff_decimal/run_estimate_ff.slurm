#!/usr/bin/bash
#SBATCH -J ff_dec
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run estimate_ff_with_higher_precision.py for one (sample, pileup) pair.
#
# Usage:
#     sbatch --job-name=ff_dec_<sample> \
#            run_estimate_ff.slurm <sample> <pileup_file> <output_dir>
#
# Inputs:
#     sample       : sample name (used as --output-prefix)
#     pileup_file  : path to the *_pileup.tsv.gz produced by snp_pileup
#     output_dir   : directory for results (<sample>_ff.tsv written here)
#
# Output:
#     <output_dir>/<sample>_ff.tsv

set -euo pipefail

sample=${1:?Usage: sbatch run_estimate_ff.slurm <sample> <pileup_file> <output_dir>}
pileup=${2:?Usage: sbatch run_estimate_ff.slurm <sample> <pileup_file> <output_dir>}
OUTPUT_DIR=${3:?Usage: sbatch run_estimate_ff.slurm <sample> <pileup_file> <output_dir>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

mkdir -p "$OUTPUT_DIR"

if [ ! -f "$pileup" ]; then
    echo "ERROR: pileup file not found: $pileup" >&2
    exit 1
fi

NCPUS="${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] sample=${sample} pileup=${pileup} output_dir=${OUTPUT_DIR} cpus=${NCPUS} host=$(hostname)"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 estimate_ff_with_higher_precision.py \
        --input-path     "$pileup" \
        --output-prefix  "${OUTPUT_DIR}/${sample}" \
        --ff-min         0 \
        --ff-max         0.3 \
        --ff-precision   0.0001 \
        --min-raw-depth  30 \
        --mode-list      "chr_only,chr_exclude,all" \
        --ncpus          "$NCPUS"

echo "[$(date)] Finished ${sample}"
