#!/usr/bin/bash
# Submit the SLURM array for random ref-40 episcore/zscore/ezscore grid search,
# optionally chaining the aggregation job once all repeats finish.
#
# Usage:
#     ./submit_grid_search_ref40.sh [-n|--dry-run] [--no-aggregate] \
#         [--input-dir <path>] [--output-base <dir>] [--total-repeats <N>] \
#         [--mode dev_test_split|all] [--min-ff <float>]
#
# Defaults:
#     input_dir   : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
#     output_base : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng
#     total_repeats : 100
#     mode        : dev_test_split
#     min-ff      : 0

set -euo pipefail

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
OUTPUT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng
TOTAL_REPEATS=100
MODE=dev_test_split
MIN_FF=0
DRY_RUN=${DRY_RUN:-0}
AGGREGATE=1

usage() { sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --no-aggregate) AGGREGATE=0; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        --total-repeats) TOTAL_REPEATS=$2; shift 2 ;;
        --mode) MODE=$2; shift 2 ;;
        --min-ff) MIN_FF=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ "$MODE" != "dev_test_split" ] && [ "$MODE" != "all" ]; then
    echo "ERROR: --mode must be dev_test_split or all (got: $MODE)" >&2
    exit 2
fi

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs

for f in episcore_grid_search.parquet zscore_grid_search.parquet meta.csv \
         ref_pool_samples.txt ezscore_ref_samples.txt; do
    if [ ! -e "${INPUT_DIR}/${f}" ]; then
        echo "ERROR: missing ${INPUT_DIR}/${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_BASE}/randomly_select_ref_40"
LAST=$((TOTAL_REPEATS - 1))

echo "Input dir      : $INPUT_DIR"
echo "Output base    : $OUTPUT_BASE"
echo "Array          : 0-${LAST} (1 repeat each, total ${TOTAL_REPEATS})"
echo "Mode           : $MODE"
echo "Min FF         : $MIN_FF"
echo "Aggregate      : $([ "$AGGREGATE" = 1 ] && echo yes || echo no)"
[ "$DRY_RUN" = 1 ] && echo "Submit mode    : DRY-RUN"
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] TOTAL_REPEATS=${TOTAL_REPEATS} MODE=${MODE} MIN_FF=${MIN_FF} \\"
    echo "             sbatch --parsable --array=0-${LAST} run_grid_search_ref40.slurm \\"
    echo "             '$INPUT_DIR' '$OUTPUT_BASE'"
    if [ "$AGGREGATE" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --dependency=afterok:<array_jobid> \\"
        echo "             run_aggregate_ref40.slurm '$OUTPUT_BASE'"
    fi
    exit 0
fi

array_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" MODE="$MODE" MIN_FF="$MIN_FF" sbatch --parsable \
    --job-name=grid_search_ref40 \
    --array="0-${LAST}" \
    run_grid_search_ref40.slurm "$INPUT_DIR" "$OUTPUT_BASE")
echo "Submitted grid_search_ref40 array job_id=${array_jobid}  logs=logs/grid_search_ref40_*.log"

if [ "$AGGREGATE" = 1 ]; then
    agg_jobid=$(sbatch --parsable \
        --job-name=aggregate_ref40 \
        --dependency="afterok:${array_jobid}" \
        run_aggregate_ref40.slurm "$OUTPUT_BASE")
    echo "Submitted aggregate_ref40 job_id=${agg_jobid} (after array ${array_jobid})  logs=logs/aggregate_ref40_*.log"
fi
