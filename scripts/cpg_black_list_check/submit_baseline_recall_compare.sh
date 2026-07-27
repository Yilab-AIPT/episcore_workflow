#!/usr/bin/bash
# Submit baseline-0.6 / baseline-0.65 pipeline + ref_free compare.
#
# Usage:
#   ./submit_baseline_recall_compare.sh [-n|--dry-run]

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs
DRY_RUN=${DRY_RUN:-0}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown: $1" >&2; exit 2 ;;
    esac
done

echo "Will run: run_baseline_recall_compare.slurm"
echo "  baseline-0.65 ← copy from pipeline/baseline"
echo "  baseline-0.6  ← nextflow aipt_ref_40 (ep 0.5/0.6)"
echo "  then ref_free 10k + plots_baseline_recall/"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode: DRY-RUN"
    exit 0
fi

jobid=$(sbatch --parsable --job-name=bl_baseline_rec run_baseline_recall_compare.slurm)
echo "Submitted bl_baseline_rec job_id=${jobid}  log=logs/bl_baseline_rec.log"
