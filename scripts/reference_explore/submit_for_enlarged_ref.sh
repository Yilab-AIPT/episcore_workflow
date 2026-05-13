#!/usr/bin/bash
# Submit a single SLURM job that runs b2z_for_enlarged_ref.py.
#
# That python script reads the fixed (threshold, recall) combo's
# _analyze_zscore.tsv.gz and _reference_zscore.tsv.gz, joins with --meta-csv,
# and for each enlarged-reference size N (10..pool_size, with --runs random
# draws per N) recomputes per-sample s_inter and reports a per-run MCC under
# <output_base>/<output_subdir>/ (default subdir: enlarged_reference).
#
# Usage:
#     ./submit_for_enlarged_ref.sh [-n|--dry-run] \
#         [--combo-dir <path>] [--meta-csv <path>] \
#         [--output-base <dir>] [--output-subdir <name>]
#
# Defaults (matches the 20260508-grid_search run):
#     combo_dir      : <RESULT_BASE>/output/threshold_0.5_recall_0.65
#     meta_csv       : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv
#     output_base    : <RESULT_BASE>/output
#     output_subdir  : enlarged_reference

set -euo pipefail

RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
META_CSV_DEFAULT=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/meta.csv

COMBO_DIR="${RESULT_BASE}/output/threshold_0.5_recall_0.65"
META_CSV="$META_CSV_DEFAULT"
OUTPUT_BASE="${RESULT_BASE}/output"
OUTPUT_SUBDIR="enlarged_reference"

DRY_RUN=${DRY_RUN:-0}

usage() {
    sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --combo-dir) COMBO_DIR=$2; shift 2 ;;
        --meta-csv) META_CSV=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        --output-subdir) OUTPUT_SUBDIR=$2; shift 2 ;;
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
mkdir -p "$OUTPUT_BASE"

echo "Submit script  : $SLURM_SCRIPT"
echo "Job name       : $JOB_NAME"
echo "Combo dir      : $COMBO_DIR"
echo "Meta csv       : $META_CSV"
echo "Output base    : $OUTPUT_BASE"
echo "Output subdir  : $OUTPUT_SUBDIR"
echo "Output dir     : $OUTPUT_BASE/$OUTPUT_SUBDIR"
echo "Log file       : logs/${JOB_NAME}.log"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode           : DRY-RUN (no jobs will be submitted)"
fi
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=${JOB_NAME} ${SLURM_SCRIPT} \\"
    echo "             '$COMBO_DIR' '$META_CSV' '$OUTPUT_BASE' '$OUTPUT_SUBDIR'"
    exit 0
fi

jobid=$(sbatch --parsable \
    --job-name="$JOB_NAME" \
    "$SLURM_SCRIPT" "$COMBO_DIR" "$META_CSV" "$OUTPUT_BASE" "$OUTPUT_SUBDIR")

echo "Submitted ${JOB_NAME}  job_id=${jobid}  log=logs/${JOB_NAME}.log"
