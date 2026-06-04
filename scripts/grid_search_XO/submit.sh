#!/usr/bin/bash
# Submit run_xo_episcore.py SLURM jobs across the recall grid.
#
#   recalls : 0.01, 0.02, ..., 0.99   (99 jobs)
#
# Recall filenames strip trailing zeros (e.g. 0.10 -> 0.1, 0.20 -> 0.2),
# so we use awk's %g formatting which matches that convention.
#
# Pass --dry-run (or set DRY_RUN=1) to print sbatch commands without submitting.
# Pass --test (or set TEST=1) to submit only recall=0.01 for smoke testing.

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
TEST=${TEST:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -t|--test) TEST=1 ;;
        -h|--help)
            sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            echo
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run] [-t|--test]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run] [-t|--test]" >&2
            exit 2
            ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"

SLURM_SCRIPT=run_xo_episcore.slurm
META_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260526-grid_search_XO
RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260526-grid_search_XO
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260525-grid_search_240k_panel_240k_model/recall_list_240k

MAX_JOBS=100
SLEEP_BETWEEN=2
SLEEP_FULL=60
JOB_PREFIX=xo_b2z_

USER_NAME=$(whoami)

mkdir -p logs

if [ ! -f "${META_DIR}/samples_meta.csv" ]; then
    echo "ERROR: ${META_DIR}/samples_meta.csv not found." >&2
    echo "Run the notebook cell that writes samples_meta.csv first." >&2
    exit 1
fi

mapfile -t RECALLS < <(awk 'BEGIN { for (i = 1; i <= 99; i += 1) printf "%g\n", i / 100 }')
if [ "$TEST" = 1 ]; then
    RECALLS=(0.01)
fi

echo "Recalls      : ${#RECALLS[@]} value(s) (${RECALLS[0]} .. ${RECALLS[-1]})"
echo "Meta file    : ${META_DIR}/samples_meta.csv"
echo "Output base  : $RESULT_BASE"
echo "Queue cap    : $MAX_JOBS jobs (filtered by name prefix '$JOB_PREFIX')"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode         : DRY-RUN (no jobs will be submitted)"
elif [ "$TEST" = 1 ]; then
    echo "Mode         : TEST (single recall job only)"
fi
echo

count_my_jobs() {
    local n
    n=$(squeue -u "$USER_NAME" -h -o '%j' 2>/dev/null \
            | grep -c "^${JOB_PREFIX}" || true)
    echo "${n:-0}"
}

wait_for_slot() {
    local n
    while :; do
        n=$(count_my_jobs)
        if [ "$n" -lt "$MAX_JOBS" ]; then
            return
        fi
        echo "  [$(date +%H:%M:%S)] queue has $n '${JOB_PREFIX}*' jobs (>= $MAX_JOBS), sleeping ${SLEEP_FULL}s..."
        sleep "$SLEEP_FULL"
    done
}

n_submitted=0
n_skipped=0
for recall in "${RECALLS[@]}"; do
    cpg_list="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
    if [ ! -f "$cpg_list" ]; then
        echo "WARN: missing $cpg_list, skipping recall=${recall}" >&2
        n_skipped=$((n_skipped + 1))
        continue
    fi

    job_name="${JOB_PREFIX}r${recall}"

    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --job-name=${job_name} ${SLURM_SCRIPT} ${recall}"
        n_submitted=$((n_submitted + 1))
        continue
    fi

    wait_for_slot
    jobid=$(sbatch --parsable \
        --job-name="$job_name" \
        "$SLURM_SCRIPT" "$recall")
    echo "Submitted recall=${recall}  job_id=${jobid}  log=logs/${job_name}.log"
    n_submitted=$((n_submitted + 1))

    sleep "$SLEEP_BETWEEN"
done

echo
if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] Would submit ${n_submitted} jobs, skipped ${n_skipped}."
else
    echo "Submitted ${n_submitted} jobs, skipped ${n_skipped}."
fi
