#!/usr/bin/bash
# Score early_ref Normals, then resubmit ref-free with updated pool/eval.
#
# Usage:
#   ./submit_ref_free_with_early_ref.sh [-n|--dry-run]

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

CHECK_OUT=${CHECK_OUT:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}
REF_FREE_ROOT=${REF_FREE_ROOT:-${CHECK_OUT}/ref_free}
DRY_RUN=${DRY_RUN:-0}

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown: $1" >&2; exit 2 ;;
    esac
done

echo "1) Score early_ref Normals for all schemes"
echo "2) Clear old ref_free outputs under $REF_FREE_ROOT (keep _smoke_*)"
echo "3) Resubmit ref_free (pool=all dev Normal, eval=dev trisomy|test, ff>0.01)"
if [ "$DRY_RUN" = 1 ]; then
    echo "Mode: DRY-RUN"
    exit 0
fi

early_job=$(sbatch --parsable --job-name=early_ref_score run_score_early_ref.slurm)
echo "Submitted early_ref_score job_id=${early_job}"

# Clear previous ref_free scheme dirs (not smoke)
if [ -d "$REF_FREE_ROOT" ]; then
    for d in "$REF_FREE_ROOT"/*; do
        base=$(basename "$d")
        case "$base" in
            _smoke*|plots) continue ;;
            *) rm -rf "$d" ;;
        esac
    done
    rm -rf "${REF_FREE_ROOT}/plots"
fi

# After early_ref scoring finishes, submit ref_free
ref_free_wrapper=$(sbatch --parsable \
    --job-name=bl_rf_launch \
    --dependency="afterok:${early_job}" \
    --wrap="cd '${SCRIPT_DIR}' && ./submit_ref_free_schemes.sh")
echo "Submitted ref_free launcher job_id=${ref_free_wrapper} dep=afterok:${early_job}"
