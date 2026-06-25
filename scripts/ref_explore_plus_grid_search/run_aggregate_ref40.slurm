#!/usr/bin/bash
#SBATCH -J aggregate_ref40
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Aggregate all repeats and select the best reference + combo.
#
# Usage:
#     sbatch run_aggregate_ref40.slurm [<output_base>]
#
# Default:
#     output_base : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng

set -euo pipefail

OUTPUT_BASE=${1:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng}
SELECT_SCORE=${SELECT_SCORE:-ezscore}
SELECT_METRIC=${SELECT_METRIC:-mean_dev_test}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "[$(date)] aggregate_ref40 host=$(hostname)"
echo "  output base  : $OUTPUT_BASE"
echo "  selection    : $SELECT_SCORE / $SELECT_METRIC"

singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 aggregate_results.py \
        --output-base "$OUTPUT_BASE" \
        --select-score "$SELECT_SCORE" \
        --select-metric "$SELECT_METRIC"

echo "[$(date)] Finished aggregate_ref40"
