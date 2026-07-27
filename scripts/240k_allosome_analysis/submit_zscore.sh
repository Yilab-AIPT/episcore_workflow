#!/usr/bin/bash
# Submit zscore recall-grid jobs (prob_class_1 cutoff 0.85 → chrX zscore).
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
SAMPLES_META="${INPUT_DIR}/zscore_samples_meta.csv"
JOB_PREFIX=allo_z_r
MAX_JOBS=25
SLEEP_BETWEEN=2
SLEEP_FULL=60
USER_NAME=$(whoami)
RESULT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples/zscore_recall

if [ ! -f "$SAMPLES_META" ]; then
    echo "ERROR: missing $SAMPLES_META — run prepare_inputs.py first" >&2
    exit 1
fi

mapfile -t RECALLS < <(awk 'BEGIN { for (i = 1; i <= 99; i += 1) printf "%g\n", i / 100 }')
if [ "$TEST" = 1 ]; then
    RECALLS=(0.01)
fi

echo "Zscore recalls : ${#RECALLS[@]} (${RECALLS[0]} .. ${RECALLS[-1]})"
echo "Samples meta   : $SAMPLES_META"
echo "Cutoff         : 0.85"
echo "Queue cap      : $MAX_JOBS (${JOB_PREFIX}*)"

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

# Output is "good" only if chrX_percentage has any positive value.
# Bad/legacy runs mapped chr 23 → chr23 (no overlap with chrX CpGs) → all zeros.
zscore_output_ok() {
    local path=$1
    [ -f "$path" ] || return 1
    singularity exec -B /lustre1,/lustre2,/appsnew \
        /lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif \
        python3 - "$path" <<'PY'
import sys
import pandas as pd
path = sys.argv[1]
df = pd.read_csv(path, sep="\t")
if "chrX_percentage" not in df.columns:
    raise SystemExit(1)
raise SystemExit(0 if float(df["chrX_percentage"].max()) > 0 else 1)
PY
}

n_submitted=0
n_skipped=0
n_rerun=0
for recall in "${RECALLS[@]}"; do
    cpg_list="${CPG_DIR}/240k_cpg_recall_${recall}.txt"
    if [ ! -f "$cpg_list" ]; then
        echo "WARN: missing $cpg_list" >&2
        n_skipped=$((n_skipped + 1))
        continue
    fi
    out="${RESULT_BASE}/recall_${recall}/_analyze_zscore.tsv.gz"
    if zscore_output_ok "$out"; then
        echo "Skip recall=${recall} (valid chrX output)"
        n_skipped=$((n_skipped + 1))
        continue
    fi
    if [ -f "$out" ]; then
        echo "Rerun recall=${recall} (chrX_percentage all-zero / invalid)"
        rm -f "$out" "${RESULT_BASE}/recall_${recall}/_reference_percentage.tsv.gz"
        n_rerun=$((n_rerun + 1))
    fi
    job_name="${JOB_PREFIX}${recall}"
    if [ "$DRY_RUN" = 1 ]; then
        echo "[DRY-RUN] sbatch --job-name=${job_name} run_zscore_recall.slurm ${recall}"
        n_submitted=$((n_submitted + 1))
        continue
    fi
    wait_for_slot
    jobid=$(sbatch --parsable --job-name="$job_name" run_zscore_recall.slurm "$recall")
    echo "Submitted recall=${recall} job_id=${jobid}"
    n_submitted=$((n_submitted + 1))
    sleep "$SLEEP_BETWEEN"
done
echo "Submitted ${n_submitted}, skipped ${n_skipped}, invalidated-for-rerun ${n_rerun}."
