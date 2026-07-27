#!/usr/bin/bash
# Submit ref-free ezscore sweeps for baseline + 6 blacklist schemes.
#
# Method (see ref_free_scheme.py / scripts/ref_free/ref_free_ezscore.py):
#   - Pool: all set=dev / Normal (any ref_type; 96 when early_ref scored → 40+40)
#   - Each repeat: re-reference scheme episcore/zscore vs first half;
#     ezscore = z-normalize(ep+z) vs second half; flag eval if any chr > cutoff
#   - Eval: (dev trisomy | all test) & ff_before_mq > 0.01 (pool excluded)
#   - 10_000 repeats, fixed combo scores from pipeline/<scheme>/merge_scores
#
# Usage:
#   ./submit_ref_free_schemes.sh [-n|--dry-run]

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260718-cpg_black_list_check}
CHECK_OUT=${CHECK_OUT:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}
PIPELINE_ROOT=${PIPELINE_ROOT:-${CHECK_OUT}/pipeline}
REF_FREE_ROOT=${REF_FREE_ROOT:-${CHECK_OUT}/ref_free}
META_CSV=${META_CSV:-${INPUT_DIR}/meta.csv}

TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
MAX_ARRAY_JOBS=${MAX_ARRAY_JOBS:-50}
REF_N=${REF_N:-40}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0.01}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF_MIN=${EZ_CUTOFF_MIN:-3.0}
EZ_CUTOFF_MAX=${EZ_CUTOFF_MAX:-4.5}
EZ_CUTOFF_STEP=${EZ_CUTOFF_STEP:-0.1}
EZ_CUTOFF=${EZ_CUTOFF:-4.5}
DRY_RUN=${DRY_RUN:-0}

SCHEMES=(baseline 15-site1-J 15-site1-M 20-site1-M 15-site2-J 15-site2-M 20-site2-M)

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
if [ "$N_JOBS" -gt "$MAX_ARRAY_JOBS" ]; then
    REPEATS_PER_JOB=$(( (TOTAL_REPEATS + MAX_ARRAY_JOBS - 1) / MAX_ARRAY_JOBS ))
    N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
fi
ARRAY_LAST=$((N_JOBS - 1))

usage() { sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ ! -f "$META_CSV" ]; then
    echo "ERROR: missing $META_CSV" >&2
    exit 1
fi
for scheme in "${SCHEMES[@]}"; do
    if [ ! -d "${PIPELINE_ROOT}/${scheme}/merge_scores" ]; then
        echo "ERROR: missing ${PIPELINE_ROOT}/${scheme}/merge_scores" >&2
        exit 1
    fi
done

mkdir -p "$REF_FREE_ROOT"

echo "Meta           : $META_CSV"
echo "Pipeline root  : $PIPELINE_ROOT"
echo "Ref-free root  : $REF_FREE_ROOT"
echo "Schemes        : ${SCHEMES[*]}"
echo "Array          : 0-${ARRAY_LAST} (${REPEATS_PER_JOB} repeats/job, total ${TOTAL_REPEATS})"
echo "ref_n (request): $REF_N  (auto-clamped to pool//2 if needed)"
echo "Eval filter    : (dev trisomy | test) & ff_before_mq > ${MIN_FF}"
echo "Pool           : set=dev & Normal (any ref_type)"
echo "Ez scatter     : cutoff=${EZ_CUTOFF}"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode           : DRY-RUN"
    exit 0
fi

export TOTAL_REPEATS REPEATS_PER_JOB REF_N SEED MIN_FF CUTOFF
export EZ_CUTOFF_MIN EZ_CUTOFF_MAX EZ_CUTOFF_STEP EZ_CUTOFF

dep_ids=()
for scheme in "${SCHEMES[@]}"; do
    scheme_out="${REF_FREE_ROOT}/${scheme}"
    mkdir -p "$scheme_out"
    jobid=$(sbatch --parsable \
        --job-name="blrf_${scheme}" \
        --array="0-${ARRAY_LAST}" \
        --export=ALL,SCHEME="${scheme}" \
        run_ref_free_scheme.slurm \
        "$META_CSV" "$PIPELINE_ROOT" "$scheme_out")
    echo "Submitted ${scheme}  job_id=${jobid}"
    dep_ids+=("${jobid}")
done

# afterok for all scheme arrays
dep_arg=$(IFS=:; echo "${dep_ids[*]}")
agg_job=$(sbatch --parsable \
    --job-name=bl_rf_agg \
    --dependency="afterok:${dep_arg}" \
    --export=ALL,EZ_CUTOFF="${EZ_CUTOFF}" \
    run_ref_free_aggregate_plot.slurm \
    "$REF_FREE_ROOT")
echo "Submitted aggregate+plot  job_id=${agg_job}  dep=afterok:${dep_arg}"
echo "Plots will land in ${REF_FREE_ROOT}/plots/"
