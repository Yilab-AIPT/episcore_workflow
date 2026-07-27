#!/usr/bin/bash
# Submit baseline-deeper-0.65 scoring + ref_free + 3-way baseline plot.
#
# Usage:
#   ./submit_baseline_deeper065.sh [-n|--dry-run]

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

echo "Will run: run_baseline_deeper065.slurm"
echo "  recall 0.65 ← replaced_deeper_recall_0.65 (same n as prod 0.65; meandiff from 0.6)"
echo "  shared -work-dir work_pipeline_shared (-resume)"
echo "  then ref_free + plots for baseline-0.6 / 0.65 / deeper-0.65"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode: DRY-RUN"
    exit 0
fi

jobid=$(sbatch --parsable --job-name=bl_deeper065 run_baseline_deeper065.slurm)
echo "Submitted bl_deeper065 job_id=${jobid}  log=logs/bl_deeper065.log"
