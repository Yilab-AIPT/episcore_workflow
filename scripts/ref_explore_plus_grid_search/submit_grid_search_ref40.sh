#!/usr/bin/bash
# Submit the SLURM array for random ref-40 episcore/zscore/ezscore grid search,
# optionally chaining the aggregation job once all repeats finish.
#
# Usage:
#     ./submit_grid_search_ref40.sh [-n|--dry-run] [--no-aggregate] \
#         [--input-dir <path>] [--output-base <dir>] [--total-repeats <N>] \
#         [--split-mode dev_test_split|even_split|all] \
#         [--combo-mode flexible|fix] [--min-ff <float>] \
#         [--finalize] [--repeats-per-job <N>] [--max-array-jobs <N>]
#
# Defaults:
#     input_dir     : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
#     output_base   : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng
#     total_repeats : 100
#     split-mode    : dev_test_split
#     combo-mode    : flexible
#     min-ff        : 0

set -euo pipefail

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
OUTPUT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-ref_40_rebuild_consider_lib_ng
TOTAL_REPEATS=100
SPLIT_MODE=dev_test_split
COMBO_MODE=flexible
MIN_FF=0
DRY_RUN=${DRY_RUN:-0}
AGGREGATE=1
FINALIZE=0
N_EZSCORE_REF=20
EZSCORE_REPEATS=5000
SEED=42
SELECT_METRIC=mean_dev_test
REPEATS_PER_JOB=100
MAX_ARRAY_JOBS=100

usage() { sed -n '2,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

compute_array_size() {
    N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
    if [ "$N_JOBS" -gt "$MAX_ARRAY_JOBS" ]; then
        REPEATS_PER_JOB=$(( (TOTAL_REPEATS + MAX_ARRAY_JOBS - 1) / MAX_ARRAY_JOBS ))
        N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
    fi
    if [ "$N_JOBS" -lt 1 ]; then
        echo "ERROR: invalid array size for total_repeats=${TOTAL_REPEATS}" >&2
        exit 2
    fi
    ARRAY_LAST=$((N_JOBS - 1))
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --no-aggregate) AGGREGATE=0; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        --total-repeats) TOTAL_REPEATS=$2; shift 2 ;;
        --split-mode) SPLIT_MODE=$2; shift 2 ;;
        --combo-mode) COMBO_MODE=$2; shift 2 ;;
        --min-ff) MIN_FF=$2; shift 2 ;;
        --finalize) FINALIZE=1; shift ;;
        --n-ezscore-ref) N_EZSCORE_REF=$2; shift 2 ;;
        --ezscore-repeats) EZSCORE_REPEATS=$2; shift 2 ;;
        --select-metric) SELECT_METRIC=$2; shift 2 ;;
        --repeats-per-job) REPEATS_PER_JOB=$2; shift 2 ;;
        --max-array-jobs) MAX_ARRAY_JOBS=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ "$SPLIT_MODE" != "dev_test_split" ] && [ "$SPLIT_MODE" != "even_split" ] \
        && [ "$SPLIT_MODE" != "all" ]; then
    echo "ERROR: --split-mode must be dev_test_split, even_split, or all (got: $SPLIT_MODE)" >&2
    exit 2
fi
if [ "$COMBO_MODE" != "flexible" ] && [ "$COMBO_MODE" != "fix" ]; then
    echo "ERROR: --combo-mode must be flexible or fix (got: $COMBO_MODE)" >&2
    exit 2
fi

if [ "$SPLIT_MODE" = "all" ]; then
    SELECT_METRIC=all
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
compute_array_size

echo "Input dir      : $INPUT_DIR"
echo "Output base    : $OUTPUT_BASE"
echo "Total repeats  : ${TOTAL_REPEATS}"
echo "Array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB} repeats each, $((ARRAY_LAST + 1)) jobs)"
echo "Split mode     : $SPLIT_MODE"
echo "Combo mode     : $COMBO_MODE"
echo "Min FF         : $MIN_FF"
echo "Select metric  : $SELECT_METRIC"
echo "Aggregate      : $([ "$AGGREGATE" = 1 ] && echo yes || echo no)"
echo "Finalize       : $([ "$FINALIZE" = 1 ] && echo yes || echo no)"
[ "$DRY_RUN" = 1 ] && echo "Submit mode    : DRY-RUN"
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] TOTAL_REPEATS=${TOTAL_REPEATS} REPEATS_PER_JOB=${REPEATS_PER_JOB} \\"
    echo "             SPLIT_MODE=${SPLIT_MODE} COMBO_MODE=${COMBO_MODE} MIN_FF=${MIN_FF} \\"
    echo "             sbatch --parsable --array=0-${ARRAY_LAST} run_grid_search_ref40.slurm \\"
    echo "             '$INPUT_DIR' '$OUTPUT_BASE'"
    if [ "$AGGREGATE" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --dependency=afterok:<array_jobid> \\"
        echo "             run_aggregate_ref40.slurm '$OUTPUT_BASE'"
        if [ "$FINALIZE" = 1 ]; then
            echo "[DRY-RUN] sbatch --parsable --dependency=afterok:<agg_jobid> \\"
            echo "             run_finalize_ref40.slurm '$OUTPUT_BASE' '$INPUT_DIR'"
        fi
    fi
    exit 0
fi

array_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" REPEATS_PER_JOB="$REPEATS_PER_JOB" \
    SPLIT_MODE="$SPLIT_MODE" COMBO_MODE="$COMBO_MODE" MIN_FF="$MIN_FF" sbatch --parsable \
    --job-name=grid_search_ref40 \
    --array="0-${ARRAY_LAST}" \
    run_grid_search_ref40.slurm "$INPUT_DIR" "$OUTPUT_BASE")
echo "Submitted grid_search_ref40 array job_id=${array_jobid}  logs=logs/grid_search_ref40_*.log"

if [ "$AGGREGATE" = 1 ]; then
    agg_jobid=$(SELECT_SCORE=ezscore SELECT_METRIC="$SELECT_METRIC" sbatch --parsable \
        --job-name=aggregate_ref40 \
        --dependency="afterok:${array_jobid}" \
        run_aggregate_ref40.slurm "$OUTPUT_BASE")
    echo "Submitted aggregate_ref40 job_id=${agg_jobid} (after array ${array_jobid})  logs=logs/aggregate_ref40_*.log"

    if [ "$FINALIZE" = 1 ]; then
        fin_jobid=$(N_EZSCORE_REF="$N_EZSCORE_REF" N_REPEATS="$EZSCORE_REPEATS" SEED="$SEED" sbatch --parsable \
            --job-name=finalize_ref40 \
            --dependency="afterok:${agg_jobid}" \
            run_finalize_ref40.slurm "$OUTPUT_BASE" "$INPUT_DIR")
        echo "Submitted finalize_ref40 job_id=${fin_jobid} (after aggregate ${agg_jobid})  logs=logs/finalize_ref40_*.log"
    fi
fi
