#!/usr/bin/bash
# Submit build_snp_subset.py jobs for recalls 0.05 .. 0.65 (step 0.05).
#
# Pass --dry-run to print sbatch commands without submitting.

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -h|--help)
            sed -n '2,4p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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

SLURM_SCRIPT=run_build_snp_subset.slurm
JOB_PREFIX=snp_r

mkdir -p logs

# 0.05, 0.10, ..., 0.65 with %g (0.10 -> 0.1) to match recall filenames.
mapfile -t RECALLS < <(awk 'BEGIN { for (r = 0.05; r <= 0.65 + 1e-9; r += 0.05) printf "%g\n", r }')

echo "Recalls : ${#RECALLS[@]} values (${RECALLS[0]} .. ${RECALLS[-1]})"
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
