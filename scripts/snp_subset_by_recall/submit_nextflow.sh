#!/usr/bin/bash
# Submit est_ff_from_pileup Nextflow jobs for recalls 0.05 .. 0.65 (step 0.05).
#
# Pass --dry-run to print sbatch commands without submitting.
# Pass --recall <value> to submit one recall only (e.g. --recall 0.05).

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
SINGLE_RECALL=""

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --recall)
            SINGLE_RECALL=${2:?--recall requires a value}
            shift 2
            ;;
        --recall=*)
            SINGLE_RECALL="${1#--recall=}"
            shift
            ;;
        -h|--help)
            sed -n '2,5p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            echo
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run] [--recall <value>]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $(basename "${BASH_SOURCE[0]}") [-n|--dry-run] [--recall <value>]" >&2
            exit 2
            ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"

SLURM_SCRIPT=run_nextflow.slurm
JOB_PREFIX=nf_r

mkdir -p logs

if [ -n "$SINGLE_RECALL" ]; then
    mapfile -t RECALLS < <(printf '%s\n' "$SINGLE_RECALL")
else
    # 0.05, 0.10, ..., 0.65 with %g (0.10 -> 0.1) to match snp_for_recall_* filenames.
    mapfile -t RECALLS < <(awk 'BEGIN { for (r = 0.05; r <= 0.65 + 1e-9; r += 0.05) printf "%g\n", r }')
fi

echo "Recalls : ${#RECALLS[@]} value(s) (${RECALLS[0]} .. ${RECALLS[-1]})"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode    : DRY-RUN"
fi
echo

n_submitted=0
for recall in "${RECALLS[@]}"; do
    job_name="${JOB_PREFIX}${recall}"
    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --parsable --job-name=${job_name} ${SLURM_SCRIPT} ${recall}"
        n_submitted=$((n_submitted + 1))
        continue
    fi
    jobid=$(sbatch --parsable --job-name="$job_name" "$SLURM_SCRIPT" "$recall")
    echo "Submitted recall=${recall}  job_id=${jobid}  log=logs/${job_name}.log"
    n_submitted=$((n_submitted + 1))
    sleep 1
done

echo
echo "Submitted ${n_submitted} job(s)."
