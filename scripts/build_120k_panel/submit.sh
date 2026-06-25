#!/usr/bin/bash
# Build a recall-0.65 CpG list with shallow sites swapped for deeper recall-0.6-only sites.
#
# Steps:
#   1. Randomly select 50 samples; write selected_samples.csv
#   2. SLURM array: compute per-site depth on recall 0.6 for each sample
#   3. SLURM job: summarize normalized depth, rank-match replace, write outputs
#
# Usage:
#     ./submit.sh [-n|--dry-run] [--n-samples 50] [--seed 42] [--skip-depth] [--replace-only] [--retry-missing]
#
# Outputs (under OUT_DIR):
#   selected_samples.csv
#   depth_per_sample/{sample}_cpg_depth.tsv.gz
#   cpg_mean_normalized_depth.tsv.gz
#   replaced_deeper_recall_0.65_sites.bed
#   replacement_pairs.tsv
#   replace_deeper_sites.log

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

OUT_DIR=${OUT_DIR:-/lustre1/cqyi/AIPT_2.0/results/small_panel}
RECALL_DIR=${RECALL_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list}
SAMPLESHEET=${SAMPLESHEET:-/lustre1/cqyi/syfan/nipt_article_plot/dev_and_test_mqres_samplesheet.csv}
CONTAINER=${CONTAINER:-/lustre1/cqyi/syfan/images/common_tools.sif}
N_SAMPLES=50
SEED=42
DRY_RUN=0
SKIP_DEPTH=0
REPLACE_ONLY=0
RETRY_MISSING=0

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --n-samples) N_SAMPLES=$2; shift 2 ;;
        --seed) SEED=$2; shift 2 ;;
        --out-dir) OUT_DIR=$2; shift 2 ;;
        --skip-depth) SKIP_DEPTH=1; shift ;;
        --replace-only) REPLACE_ONLY=1; shift ;;
        --retry-missing) RETRY_MISSING=1; REPLACE_ONLY=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

mkdir -p logs "${OUT_DIR}/depth_per_sample"

run_singularity() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry-run] singularity exec ... $*"
    else
        singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" "$@"
    fi
}

echo "=== build_120k_panel submit ==="
echo "  script_dir : $SCRIPT_DIR"
echo "  out_dir    : $OUT_DIR"
echo "  n_samples  : $N_SAMPLES"
echo "  seed       : $SEED"
echo "  container  : $CONTAINER"

if [ "$REPLACE_ONLY" -eq 0 ] && [ "$RETRY_MISSING" -eq 0 ]; then
    echo "[step 1] select ${N_SAMPLES} samples"
    run_singularity python3 "${SCRIPT_DIR}/replace_deeper_sites.py" \
        --samplesheet "$SAMPLESHEET" \
        --out-dir "$OUT_DIR" \
        --n-samples "$N_SAMPLES" \
        --seed "$SEED" \
        --selected-samples "${OUT_DIR}/selected_samples.csv" \
        --select-only
elif [ "$RETRY_MISSING" -eq 1 ]; then
    echo "[step 1] skipped (retry missing depth files)"
else
    echo "[step 1] skipped (replace-only)"
fi

n_selected=$(tail -n +2 "${OUT_DIR}/selected_samples.csv" | wc -l | tr -d ' ')
echo "  selected   : ${n_selected} samples"

missing_indices=()
while IFS= read -r sample_id; do
    depth_file="${OUT_DIR}/depth_per_sample/${sample_id}_cpg_depth.tsv.gz"
    if [ ! -f "$depth_file" ]; then
        idx=$(grep -n "^${sample_id}," "${OUT_DIR}/selected_samples.csv" | head -1 | cut -d: -f1)
        missing_indices+=($((idx - 1)))
    fi
done < <(tail -n +2 "${OUT_DIR}/selected_samples.csv" | cut -d, -f1)

if [ "${#missing_indices[@]}" -gt 0 ]; then
    echo "  missing depth: ${#missing_indices[@]} sample(s)"
    array_spec=$(IFS=,; echo "${missing_indices[*]}")
else
    echo "  missing depth: 0"
    array_spec=""
fi

if [ "$SKIP_DEPTH" -eq 0 ] && [ "$REPLACE_ONLY" -eq 0 ]; then
    if [ "$RETRY_MISSING" -eq 1 ]; then
        if [ -z "$array_spec" ]; then
            echo "[step 2] all depth files present; skipping depth array"
            DEPTH_JOB_ID=
        else
            echo "[step 2] submit depth array for missing tasks (${array_spec})"
            depth_cmd=(sbatch --array="${array_spec}" --export=ALL,OUT_DIR="${OUT_DIR}" \
                "${SCRIPT_DIR}/run_depth_array.slurm")
            if [ "$DRY_RUN" -eq 1 ]; then
                echo "[dry-run] ${depth_cmd[*]}"
                DEPTH_JOB_ID=0
            else
                depth_out=$("${depth_cmd[@]}")
                echo "$depth_out"
                DEPTH_JOB_ID=$(sed -n 's/.*Job \([0-9][0-9]*\).*/\1/p' <<<"$depth_out")
            fi
        fi
    else
        echo "[step 2] submit depth array (1-${n_selected})"
        depth_cmd=(sbatch --array="1-${n_selected}" --export=ALL,OUT_DIR="${OUT_DIR}" \
            "${SCRIPT_DIR}/run_depth_array.slurm")
        if [ "$DRY_RUN" -eq 1 ]; then
            echo "[dry-run] ${depth_cmd[*]}"
            DEPTH_JOB_ID=0
        else
            depth_out=$("${depth_cmd[@]}")
            echo "$depth_out"
            DEPTH_JOB_ID=$(sed -n 's/.*Job \([0-9][0-9]*\).*/\1/p' <<<"$depth_out")
        fi
    fi
else
    echo "[step 2] skipped depth array"
    DEPTH_JOB_ID=
fi

echo "[step 3] submit merge + replace"
replace_cmd=(sbatch --export=ALL,OUT_DIR="${OUT_DIR}",N_SAMPLES="${N_SAMPLES}",SEED="${SEED}" \
    "${SCRIPT_DIR}/run_replace_panel.slurm")
if [ -n "${DEPTH_JOB_ID:-}" ] && [ "$DEPTH_JOB_ID" != "0" ]; then
    replace_cmd=(sbatch --dependency=afterok:"${DEPTH_JOB_ID}" \
        --export=ALL,OUT_DIR="${OUT_DIR}",N_SAMPLES="${N_SAMPLES}",SEED="${SEED}" \
        "${SCRIPT_DIR}/run_replace_panel.slurm")
fi
if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry-run] ${replace_cmd[*]}"
else
    "${replace_cmd[@]}"
fi

echo "Done."
