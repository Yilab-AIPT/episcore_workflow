#!/usr/bin/bash
# Submit a single SLURM job that runs b2z_for_best_combo.py.
#
# That python script reads <best_combo_csv> (one row per chromosome with the
# best (threshold, recall)), pulls the per-chr aggregated beta values and
# CpG counts that beta_to_episcore.py already wrote under
# <output_base>/threshold_<t>_recall_<r>/, recomputes z_intra/z_inter/s_inter,
# and writes <output_base>/best_combo/_{reference,analyze}_zscore.tsv.gz.
#
# Usage:
#     ./submit_for_best_combo.sh [-n|--dry-run] \
#         [--best-combo-csv <path>] [--output-base <dir>]
#
# Defaults (matches the 20260508-grid_search run):
#     best_combo_csv : <RESULT_BASE>/best_combo_per_chr.csv
#     output_base    : <RESULT_BASE>/output

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
BEST_COMBO_CSV="${RESULT_BASE}/best_combo_per_chr.csv"
OUTPUT_BASE="${RESULT_BASE}/output"

DRY_RUN=${DRY_RUN:-0}

usage() {
    sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --best-combo-csv) BEST_COMBO_CSV=$2; shift 2 ;;
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

SLURM_SCRIPT=run_b2z_for_best_combo.slurm
JOB_NAME=b2z_best_combo

mkdir -p logs

if [ ! -f "$BEST_COMBO_CSV" ]; then
    echo "ERROR: best-combo CSV not found: $BEST_COMBO_CSV" >&2
    exit 1
fi
if [ ! -d "$OUTPUT_BASE" ]; then
    echo "ERROR: output base dir not found: $OUTPUT_BASE" >&2
    exit 1
fi

echo "Submit script : $SLURM_SCRIPT"
echo "Job name      : $JOB_NAME"
echo "Best-combo csv: $BEST_COMBO_CSV"
echo "Output base   : $OUTPUT_BASE"
echo "Output dir    : $OUTPUT_BASE/best_combo"
echo "Log file      : logs/${JOB_NAME}.log"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode          : DRY-RUN (no jobs will be submitted)"
fi
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=${JOB_NAME} ${SLURM_SCRIPT} \\"
    echo "             '$BEST_COMBO_CSV' '$OUTPUT_BASE'"
    exit 0
fi

jobid=$(sbatch --parsable \
    --job-name="$JOB_NAME" \
    "$SLURM_SCRIPT" "$BEST_COMBO_CSV" "$OUTPUT_BASE")

echo "Submitted ${JOB_NAME}  job_id=${jobid}  log=logs/${JOB_NAME}.log"
