#!/usr/bin/bash
# Submit episcore recall-grid jobs (threshold 0.5 betas → chrX s_inter).
#
#   recalls : 0.01 .. 0.99  (99 jobs)   |  --test → recall=0.01 only
#
# Pass --dry-run to print sbatch commands without submitting.

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
TEST=${TEST:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -t|--test) TEST=1 ;;
        -h|--help)
            sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "Unknown argument: $arg" >&2; exit 2 ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_240k
SAMPLES_META="${INPUT_DIR}/episcore_samples_meta.csv"
JOB_PREFIX=allo_epi_r
MAX_JOBS=40
SLEEP_BETWEEN=2
SLEEP_FULL=60
USER_NAME=$(whoami)
RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples/episcore_recall

if [ ! -f "$SAMPLES_META" ]; then
    echo "ERROR: missing $SAMPLES_META — run prepare_inputs.py after Nextflow betas exist" >&2
    exit 1
fi

mapfile -t RECALLS < <(awk 'BEGIN { for (i = 1; i <= 99; i += 1) printf "%g\n", i / 100 }')
if [ "$TEST" = 1 ]; then
    RECALLS=(0.01)
fi

echo "Episcore recalls : ${#RECALLS[@]} (${RECALLS[0]} .. ${RECALLS[-1]})"
echo "Samples meta     : $SAMPLES_META"
echo "Queue cap        : $MAX_JOBS (${JOB_PREFIX}*)"

count_my_jobs() {
    squeue -u "$USER_NAME" -h -o '%j' 2>/dev/null | grep -c "^${JOB_PREFIX}" || true
}

wait_for_slot() {
    while :; do
        n=$(count_my_jobs)
        n=${n:-0}
        if [ "$n" -lt "$MAX_JOBS" ]; then
            return
        fi
        echo "  [$(date +%H:%M:%S)] ${n} jobs queued; sleep ${SLEEP_FULL}s"
        sleep "$SLEEP_FULL"
    done
}

n_submitted=0
n_skipped=0
for recall in "${RECALLS[@]}"; do
    cpg_list="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
    if [ ! -f "$cpg_list" ]; then
        echo "WARN: missing $cpg_list" >&2
        n_skipped=$((n_skipped + 1))
        continue
    fi
    if [ -f "${RESULT_BASE}/recall_${recall}/_analyze_zscore.tsv.gz" ]; then
        echo "Skip recall=${recall} (output exists)"
        n_skipped=$((n_skipped + 1))
        continue
    fi
    job_name="${JOB_PREFIX}${recall}"
    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --job-name=${job_name} run_episcore_recall.slurm ${recall}"
        n_submitted=$((n_submitted + 1))
        continue
    fi
    wait_for_slot
    jobid=$(sbatch --parsable --job-name="$job_name" run_episcore_recall.slurm "$recall")
    echo "Submitted recall=${recall} job_id=${jobid}"
    n_submitted=$((n_submitted + 1))
    sleep "$SLEEP_BETWEEN"
done
echo "Submitted ${n_submitted}, skipped ${n_skipped}."
