#!/usr/bin/bash
# Submit a single SLURM job that runs b2z_for_enlarged_ref.py.
#
# Modes (--mode):
#   all (default) — original enlarged-reference sweep; writes report.csv.
#   isolated      — dev/test MCC with pool sampling; writes report.csv and
#                   top_100_reference_list.csv.
#
# Usage:
#     ./submit_for_enlarged_ref.sh [-n|--dry-run] \
#         [--combo-dir <path>] [--meta-csv <path>] \
#         [--output-dir <dir>] [--mode all|isolated]
#
# Defaults (matches the 20260508-grid_search run):
#     combo_dir   : <RESULT_BASE>/output/threshold_0.5_recall_0.65
#     meta_csv    : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv
#     output_dir  : <RESULT_BASE>/output/enlarged_reference  (all)
#                   <RESULT_BASE>/output/enlarged_reference_isolated  (isolated)
#     mode        : all

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
META_CSV_DEFAULT=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv

COMBO_DIR="${RESULT_BASE}/output/threshold_0.5_recall_0.65"
META_CSV="$META_CSV_DEFAULT"
MODE="all"
OUTPUT_DIR=""

DRY_RUN=${DRY_RUN:-0}

usage() {
    sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

default_output_dir() {
    case "$MODE" in
        all) echo "${RESULT_BASE}/output/enlarged_reference" ;;
        isolated) echo "${RESULT_BASE}/output/enlarged_reference_isolated" ;;
        *)
            echo "ERROR: mode must be 'all' or 'isolated', got: $MODE" >&2
            exit 1
            ;;
    esac
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --combo-dir) COMBO_DIR=$2; shift 2 ;;
        --meta-csv) META_CSV=$2; shift 2 ;;
        --output-dir) OUTPUT_DIR=$2; shift 2 ;;
        --mode) MODE=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$MODE" in
    all|isolated) ;;
    *)
        echo "ERROR: mode must be 'all' or 'isolated', got: $MODE" >&2
        exit 1
        ;;
esac

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(default_output_dir)
fi

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"

SLURM_SCRIPT=run_b2z_for_enlarged_ref.slurm
JOB_NAME=b2z_enlarged_ref

mkdir -p logs

if [ ! -d "$COMBO_DIR" ]; then
    echo "ERROR: combo dir not found: $COMBO_DIR" >&2
    exit 1
fi
if [ ! -f "$META_CSV" ]; then
    echo "ERROR: meta CSV not found: $META_CSV" >&2
    exit 1
fi
mkdir -p "$(dirname "$OUTPUT_DIR")"

echo "Submit script  : $SLURM_SCRIPT"
echo "Job name       : $JOB_NAME"
echo "Combo dir      : $COMBO_DIR"
echo "Meta csv       : $META_CSV"
echo "Output dir     : $OUTPUT_DIR"
echo "Mode           : $MODE"
echo "Log file       : logs/${JOB_NAME}.log"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode (submit)  : DRY-RUN (no jobs will be submitted)"
fi
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=${JOB_NAME} ${SLURM_SCRIPT} \\"
    echo "             '$COMBO_DIR' '$META_CSV' '$OUTPUT_DIR' '$MODE'"
    exit 0
fi

jobid=$(sbatch --parsable \
    --job-name="$JOB_NAME" \
    "$SLURM_SCRIPT" "$COMBO_DIR" "$META_CSV" "$OUTPUT_DIR" "$MODE")

echo "Submitted ${JOB_NAME}  job_id=${jobid}  log=logs/${JOB_NAME}.log"
