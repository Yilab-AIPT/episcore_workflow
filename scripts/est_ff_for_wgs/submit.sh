#!/usr/bin/bash
# Submit WGS FF estimation jobs (pileup + FF), one per samplesheet row.
#
# Samplesheet format:
#     sample,clean_bam
#
# Each row is keyed by the BAM basename (without .bam) so duplicate sample names
# with different BAMs get distinct job names and outputs.
#
# Pass --dry-run (or set DRY_RUN=1) to print sbatch commands without submitting.

set -euo pipefail

SAMPLESHEET=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-est_ff_for_wgs/samplesheet.csv
OUTPUT_DIR=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-est_ff_for_wgs

DRY_RUN=${DRY_RUN:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -h|--help)
            sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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

SLURM_SCRIPT=run_sample.slurm
MAX_JOBS=50
SLEEP_BETWEEN=2
SLEEP_FULL=60
JOB_PREFIX=wgs_ff_
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
echo "Queue cap    : $MAX_JOBS jobs (prefix '${JOB_PREFIX}')"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode         : DRY-RUN"
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

HEADER=$(head -n 1 "$SAMPLESHEET")
BAM_IDX=2
IFS=',' read -r -a HEADER_COLS <<< "$HEADER"
for i in "${!HEADER_COLS[@]}"; do
    col=$(echo "${HEADER_COLS[$i]}" | tr -d '[:space:]"\r')
    case "$col" in
        clean_bam|bam) BAM_IDX=$((i + 1)) ;;
    esac
done

while IFS=, read -r -a fields; do
    bam=$(echo "${fields[$((BAM_IDX - 1))]:-}" | tr -d '"\r' | xargs || true)
    if [ -z "${bam}" ]; then
        continue
    fi
    if [ ! -f "$bam" ]; then
        echo "WARN: BAM not found: ${bam}" >&2
        n_skipped=$((n_skipped + 1))
        continue
    fi

    bam_id=$(basename "$bam" .bam)
    out_file="${OUTPUT_DIR}/${bam_id}_ff.tsv"
    if [ -f "$out_file" ]; then
        echo "SKIP: ${bam_id} (output exists: ${out_file})"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    job_name="${JOB_PREFIX}${bam_id}"

    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --job-name=${job_name} ${SLURM_SCRIPT} ${bam_id} ${bam} ${OUTPUT_DIR}"
        n_submitted=$((n_submitted + 1))
        continue
    fi

    wait_for_slot
    jobid=$(sbatch --parsable \
        --job-name="$job_name" \
        "$SLURM_SCRIPT" "$bam_id" "$bam" "$OUTPUT_DIR")
    echo "Submitted ${bam_id}  job_id=${jobid}  log=logs/${job_name}.log"
    n_submitted=$((n_submitted + 1))
    sleep "$SLEEP_BETWEEN"
done < <(tail -n +2 "$SAMPLESHEET")

echo
if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] Would submit ${n_submitted} jobs, skipped ${n_skipped}."
else
    echo "Submitted ${n_submitted} jobs, skipped ${n_skipped}."
fi
