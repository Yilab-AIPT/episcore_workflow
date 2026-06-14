#!/usr/bin/bash
# Submit SLURM array job for 1000 random ref-40 score repeats.
#
# Usage:
#     ./submit_random_ref40_scores.sh [-n|--dry-run] [--input-dir <path>] [--output-base <dir>]
#
# Defaults:
#     input_dir   : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40
#     output_base : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40

set -euo pipefail

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40
OUTPUT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260607-ref_40
DRY_RUN=${DRY_RUN:-0}

usage() {
    sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"

SLURM_SCRIPT=run_random_ref40_scores.slurm
JOB_NAME=random_ref40_scores

mkdir -p logs

for f in meta.csv beta.csv percentage.csv; do
    if [ ! -f "${INPUT_DIR}/${f}" ]; then
        echo "ERROR: missing ${INPUT_DIR}/${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_BASE}/randomly_select_ref_40"

echo "Submit script  : $SLURM_SCRIPT"
echo "Job name       : $JOB_NAME"
echo "Input dir      : $INPUT_DIR"
echo "Output base    : $OUTPUT_BASE"
echo "Array          : 0-19 (50 repeats each, total 1000)"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode (submit)  : DRY-RUN"
fi
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=${JOB_NAME} ${SLURM_SCRIPT} \\"
    echo "             '$INPUT_DIR' '$OUTPUT_BASE'"
    exit 0
fi

jobid=$(sbatch --parsable \
    --job-name="$JOB_NAME" \
    "$SLURM_SCRIPT" "$INPUT_DIR" "$OUTPUT_BASE")

echo "Submitted ${JOB_NAME} array  job_id=${jobid}  logs=logs/${JOB_NAME}_*.log"
