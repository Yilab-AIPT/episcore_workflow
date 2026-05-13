#!/usr/bin/bash
# Submit estimate_ff_with_higher_precision.py SLURM jobs, one per sample.
#
# The script parses a CSV samplesheet of the form:
#     sample,pileup_file
#     HCPT0121,/lustre1/.../HCPT0121P_pileup.tsv.gz
#     ...
# and submits one job per row using run_estimate_ff.slurm.
#
# To avoid swamping the scheduler we throttle submissions in two ways:
#   (a) wait until our queued/running jobs are < MAX_JOBS before submitting
#   (b) sleep SLEEP_BETWEEN seconds after every submission
#
# Pass --dry-run (or set DRY_RUN=1) to print the sbatch commands without
# submitting anything.

set -euo pipefail

SAMPLESHEET=/lustre1/cqyi/syfan/nipt_article_plot/snp_pileup_result_samplesheet.csv
OUTPUT_DIR=/lustre1/cqyi/syfan/nipt_article_plot/dev_and_test_ff

DRY_RUN=${DRY_RUN:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -h|--help)
            sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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

SLURM_SCRIPT=run_estimate_ff.slurm

MAX_JOBS=100          # cap on simultaneously queued/running jobs (this user, this prefix)
SLEEP_BETWEEN=2       # seconds between consecutive sbatch calls
SLEEP_FULL=60         # seconds to wait when the queue is full
JOB_PREFIX=ff_dec_    # only count our own jobs against MAX_JOBS

USER_NAME=$(whoami)

mkdir -p logs
mkdir -p "$OUTPUT_DIR"

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi
if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "ERROR: slurm script not found: $SLURM_SCRIPT" >&2
    exit 1
fi

echo "Samplesheet  : $SAMPLESHEET"
echo "Output dir   : $OUTPUT_DIR"
echo "Slurm script : $SLURM_SCRIPT"
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

# Read header to find the right columns; fall back to "sample,pileup_file" order.
HEADER=$(head -n 1 "$SAMPLESHEET")
SAMPLE_IDX=1
PILEUP_IDX=2
IFS=',' read -r -a HEADER_COLS <<< "$HEADER"
for i in "${!HEADER_COLS[@]}"; do
    col=$(echo "${HEADER_COLS[$i]}" | tr -d '[:space:]"\r')
    case "$col" in
        sample)       SAMPLE_IDX=$((i + 1)) ;;
        pileup_file)  PILEUP_IDX=$((i + 1)) ;;
    esac
done

# Iterate over data rows, picking out the resolved columns.
while IFS=, read -r -a fields; do
    sample=$(echo "${fields[$((SAMPLE_IDX - 1))]:-}" | tr -d '"\r' | xargs || true)
    pileup=$(echo "${fields[$((PILEUP_IDX - 1))]:-}" | tr -d '"\r' | xargs || true)

    if [ -z "${sample}" ] || [ -z "${pileup}" ]; then
        continue
    fi
    if [ ! -f "$pileup" ]; then
        echo "WARN: pileup not found for ${sample}: ${pileup}" >&2
        n_skipped=$((n_skipped + 1))
        continue
    fi

    # Skip samples whose output already exists (idempotent re-runs).
    out_file="${OUTPUT_DIR}/${sample}_ff.tsv"
    if [ -f "$out_file" ]; then
        echo "SKIP: ${sample} (output already exists: ${out_file})"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    job_name="${JOB_PREFIX}${sample}"

    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --job-name=${job_name} ${SLURM_SCRIPT} ${sample} ${pileup}"
        n_submitted=$((n_submitted + 1))
        continue
    fi

    wait_for_slot
    jobid=$(sbatch --parsable \
        --job-name="$job_name" \
        "$SLURM_SCRIPT" "$sample" "$pileup")
    echo "Submitted ${sample}  job_id=${jobid}  log=logs/${job_name}.log"
    n_submitted=$((n_submitted + 1))

    sleep "$SLEEP_BETWEEN"
done < <(tail -n +2 "$SAMPLESHEET")

echo
if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] Would submit ${n_submitted} jobs, skipped ${n_skipped}."
else
    echo "Submitted ${n_submitted} jobs, skipped ${n_skipped}."
fi
