#!/usr/bin/bash
# Submit beta_to_episcore.py SLURM jobs across the (threshold x recall) grid.
#
#   thresholds : 0.1, 0.33, 0.5, 0.67, 0.9       (5 values)
#   recalls    : 0.01, 0.02, ..., 0.99           (99 values)
#   total      : 5 * 99 = 495 jobs
#
# Recall filenames strip trailing zeros (e.g. 0.10 -> 0.1, 0.20 -> 0.2),
# so we use awk's %g formatting which matches that convention.
#
# To avoid swamping the scheduler we throttle submissions in two ways:
#   (a) wait until our queued/running jobs are < MAX_JOBS before submitting
#   (b) sleep SLEEP_BETWEEN seconds after every submission
#
# Pass --dry-run (or set DRY_RUN=1) to print the sbatch commands without
# submitting anything. In dry-run mode the queue check and the inter-job
# sleep are skipped so you get the full plan immediately.

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -h|--help)
            sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            echo
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run]" >&2
            exit 2
            ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"

SLURM_SCRIPT=run_beta_to_episcore.slurm
RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search
CPG_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list

THRESHOLDS=(0.1 0.33 0.5 0.67 0.9)

MAX_JOBS=100          # cap on simultaneously queued/running jobs (this user, this prefix)
SLEEP_BETWEEN=2       # seconds between consecutive sbatch calls
SLEEP_FULL=60         # seconds to wait when the queue is full
JOB_PREFIX=b2z_       # only count our own jobs against MAX_JOBS

USER_NAME=$(whoami)

mkdir -p logs

# Generate recall values 0.01, 0.02, ..., 0.99 with %g (no trailing zeros).
mapfile -t RECALLS < <(awk 'BEGIN { for (i = 1; i <= 99; i += 1) printf "%g\n", i / 100 }')

echo "Thresholds   : ${THRESHOLDS[*]}"
echo "Recalls      : ${#RECALLS[@]} values (${RECALLS[0]} .. ${RECALLS[-1]})"
echo "Total combos : $((${#THRESHOLDS[@]} * ${#RECALLS[@]}))"
echo "Queue cap    : $MAX_JOBS jobs (filtered by name prefix '$JOB_PREFIX')"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode         : DRY-RUN (no jobs will be submitted)"
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
for thres in "${THRESHOLDS[@]}"; do
    ref_dir="${RESULT_BASE}/threshold_${thres}/early_ref_n_17"
    analyze_dir="${RESULT_BASE}/threshold_${thres}/analyze_n_283"
    if [ ! -d "$ref_dir" ] || [ ! -d "$analyze_dir" ]; then
        echo "WARN: missing ref/analyze dirs under threshold_${thres}, skipping all recalls" >&2
        n_skipped=$((n_skipped + ${#RECALLS[@]}))
        continue
    fi

    for recall in "${RECALLS[@]}"; do
        cpg_list="${CPG_DIR}/220k_cpg_recall_${recall}.txt"
        if [ ! -f "$cpg_list" ]; then
            echo "WARN: missing $cpg_list, skipping thres=${thres} recall=${recall}" >&2
            n_skipped=$((n_skipped + 1))
            continue
        fi

        job_name="${JOB_PREFIX}t${thres}_r${recall}"

        if [ "$DRY_RUN" = 1 ]; then
            echo "[DRY-RUN] sbatch --parsable --job-name=${job_name} ${SLURM_SCRIPT} ${thres} ${recall}"
            n_submitted=$((n_submitted + 1))
            continue
        fi

        wait_for_slot
        jobid=$(sbatch --parsable \
            --job-name="$job_name" \
            "$SLURM_SCRIPT" "$thres" "$recall")
        echo "Submitted t=${thres} r=${recall}  job_id=${jobid}  log=logs/${job_name}.log"
        n_submitted=$((n_submitted + 1))

        sleep "$SLEEP_BETWEEN"
    done
done

echo
if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] Would submit ${n_submitted} jobs, skipped ${n_skipped}."
else
    echo "Submitted ${n_submitted} jobs, skipped ${n_skipped}."
fi
