#!/usr/bin/bash
# Submit aipt_ref_40 vs meta comparison (FF skipped).
#
# Usage:
#   ./submit_aipt_ref40_check.sh [-n|--dry-run] [--input-dir <path>] [--output-dir <dir>]

set -euo pipefail

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260716-ezscore_check}
OUTPUT_DIR=${OUTPUT_DIR:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260716-ezscore_check}
DRY_RUN=${DRY_RUN:-0}

usage() { sed -n '2,6p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-dir) OUTPUT_DIR=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs "$OUTPUT_DIR"

for f in meta.csv samplesheet.csv; do
    if [ ! -e "${INPUT_DIR}/${f}" ]; then
        echo "ERROR: missing ${INPUT_DIR}/${f}" >&2
        exit 1
    fi
done

echo "Submit script : run_aipt_ref40_check.slurm"
echo "Input dir     : $INPUT_DIR"
echo "Output dir    : $OUTPUT_DIR"
echo "Log file      : logs/ezscore_check.log"
echo "Combo         : episcore 0.5/0.65 ; zscore 0.85/0.95 ; skip_ff=true"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode          : DRY-RUN"
fi
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=ezscore_check run_aipt_ref40_check.slurm '$INPUT_DIR' '$OUTPUT_DIR'"
    exit 0
fi

jobid=$(sbatch --parsable \
    --job-name=ezscore_check \
    run_aipt_ref40_check.slurm "$INPUT_DIR" "$OUTPUT_DIR")
echo "Submitted ezscore_check  job_id=${jobid}  log=logs/ezscore_check.log"
