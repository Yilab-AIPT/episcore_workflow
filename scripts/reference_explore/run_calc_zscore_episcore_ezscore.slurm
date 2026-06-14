#!/usr/bin/bash
#SBATCH -J calc_ref40_scores
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Compute zscore / episcore / ezscore for the 20260607-ref_40 reference explore run.
#
# Usage:
#     sbatch run_calc_zscore_episcore_ezscore.slurm [<input_dir>] [<output_dir>]
#
# Defaults:
#     input_dir  : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40
#     output_dir : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40}
OUTPUT_DIR=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p logs

for f in meta.csv beta.csv percentage.csv; do
    if [ ! -f "${INPUT_DIR}/${f}" ]; then
        echo "ERROR: missing ${INPUT_DIR}/${f}" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] Starting calc_zscore_episcore_ezscore on host $(hostname)"
echo "  input dir  : $INPUT_DIR"
echo "  output dir : $OUTPUT_DIR"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 calc_zscore_episcore_ezscore.py \
        --input-dir  "$INPUT_DIR" \
        --output-dir "$OUTPUT_DIR"

echo "[$(date)] Finished calc_zscore_episcore_ezscore"
