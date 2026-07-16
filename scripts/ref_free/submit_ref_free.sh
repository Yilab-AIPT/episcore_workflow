#!/usr/bin/bash
# Submit the SLURM array for reference-free episcore/zscore abnormality sweep,
# optionally chaining aggregation once all repeats finish.
#
# Usage:
#     ./submit_ref_free.sh [-n|--dry-run] [--no-aggregate] \
#         [--input-dir <path>] [--output-base <dir>] [--total-repeats <N>] \
#         [--repeats-per-job <N>] [--max-array-jobs <N>]
#
# Defaults:
#     input_dir       : /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
#     output_base     : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260713-ref_free
#     total_repeats   : 10000
#     max_array_jobs  : 50

set -euo pipefail

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
OUTPUT_BASE=${OUTPUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260713-ref_free}
TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REF_N=${REF_N:-50}
CUTOFF=${CUTOFF:-3.0}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0}
COMBO_MODE=${COMBO_MODE:-all}
EP_THRESHOLD=${EP_THRESHOLD:-}
EP_RECALL=${EP_RECALL:-}
Z_THRESHOLD=${Z_THRESHOLD:-}
Z_RECALL=${Z_RECALL:-}
DRY_RUN=${DRY_RUN:-0}
AGGREGATE=1
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
MAX_ARRAY_JOBS=${MAX_ARRAY_JOBS:-50}

usage() { sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

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
        --ref-n) REF_N=$2; shift 2 ;;
        --combo-mode) COMBO_MODE=$2; shift 2 ;;
        --ep-threshold) EP_THRESHOLD=$2; shift 2 ;;
        --ep-recall) EP_RECALL=$2; shift 2 ;;
        --z-threshold) Z_THRESHOLD=$2; shift 2 ;;
        --z-recall) Z_RECALL=$2; shift 2 ;;
        --cutoff) CUTOFF=$2; shift 2 ;;
        --seed) SEED=$2; shift 2 ;;
        --min-ff) MIN_FF=$2; shift 2 ;;
        --repeats-per-job) REPEATS_PER_JOB=$2; shift 2 ;;
        --max-array-jobs) MAX_ARRAY_JOBS=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs

for f in episcore_grid_search.parquet zscore_grid_search.parquet meta.csv; do
    if [ ! -e "${INPUT_DIR}/${f}" ]; then
        echo "ERROR: missing ${INPUT_DIR}/${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_BASE}/ref_free"
compute_array_size

echo "Input dir      : $INPUT_DIR"
echo "Output base    : $OUTPUT_BASE"
echo "Total repeats  : ${TOTAL_REPEATS}"
echo "Array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB} repeats each, $((ARRAY_LAST + 1)) jobs)"
echo "Ref N          : $REF_N"
echo "Combo mode     : $COMBO_MODE"
if [ "$COMBO_MODE" = "fixed" ]; then
    echo "Episcore combo : threshold=$EP_THRESHOLD recall=$EP_RECALL"
    echo "Zscore combo   : threshold=$Z_THRESHOLD recall=$Z_RECALL"
else
    echo "Combos         : all available in episcore + zscore parquets"
fi
echo "Cutoff         : $CUTOFF"
echo "Aggregate      : $([ "$AGGREGATE" = 1 ] && echo yes || echo no)"
[ "$DRY_RUN" = 1 ] && echo "Submit mode    : DRY-RUN"
echo

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] TOTAL_REPEATS=${TOTAL_REPEATS} REPEATS_PER_JOB=${REPEATS_PER_JOB} \\"
    echo "             sbatch --parsable --array=0-${ARRAY_LAST} run_ref_free.slurm \\"
    echo "             '$INPUT_DIR' '$OUTPUT_BASE'"
    if [ "$AGGREGATE" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --dependency=afterok:<array_jobid> \\"
        echo "             run_aggregate_ref_free.slurm '$OUTPUT_BASE'"
    fi
    exit 0
fi

array_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" REPEATS_PER_JOB="$REPEATS_PER_JOB" \
    REF_N="$REF_N" COMBO_MODE="$COMBO_MODE" \
    EP_THRESHOLD="$EP_THRESHOLD" EP_RECALL="$EP_RECALL" \
    Z_THRESHOLD="$Z_THRESHOLD" Z_RECALL="$Z_RECALL" \
    CUTOFF="$CUTOFF" SEED="$SEED" MIN_FF="$MIN_FF" sbatch --parsable \
    --job-name=ref_free \
    --array="0-${ARRAY_LAST}" \
    run_ref_free.slurm "$INPUT_DIR" "$OUTPUT_BASE")
echo "Submitted ref_free array job_id=${array_jobid}  logs=logs/ref_free_*.log"

if [ "$AGGREGATE" = 1 ]; then
    agg_jobid=$(TOTAL_REPEATS="$TOTAL_REPEATS" sbatch --parsable \
        --job-name=aggregate_ref_free \
        --dependency="afterok:${array_jobid}" \
        run_aggregate_ref_free.slurm "$OUTPUT_BASE")
    echo "Submitted aggregate_ref_free job_id=${agg_jobid} (after array ${array_jobid})  logs=logs/aggregate_ref_free_*.log"
fi
