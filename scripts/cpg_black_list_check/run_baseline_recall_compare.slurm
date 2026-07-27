#!/usr/bin/bash
#SBATCH -J bl_baseline_rec
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Build baseline-0.6 and baseline-0.65 (pipeline + ref_free + plots).
#
# - baseline-0.65: reuse existing pipeline/baseline merge_scores (already ep 0.65)
# - baseline-0.6: run aipt_ref_40 with ep combo 0.5/0.6 on all samples
# - then ref_free (10k) for both and write scatters under ref_free/plots_baseline_recall/
#
# Usage:
#   sbatch run_baseline_recall_compare.slurm

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

CHECK_OUT=${CHECK_OUT:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}
INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260718-cpg_black_list_check}
PIPELINE_ROOT=${PIPELINE_ROOT:-${CHECK_OUT}/pipeline}
RECALL_DIR=${RECALL_DIR:-${INPUT_DIR}/recall_list}
SAMPLESHEET=${SAMPLESHEET:-${CHECK_OUT}/prepared/samplesheet_all.csv}
WORKDIR=${WORKDIR:-${CHECK_OUT}/work_pipeline_shared}
SOURCE_ASSETS=${SOURCE_ASSETS:-${CHECK_OUT}/grid_search_assets}
REF_FREE_ROOT=${REF_FREE_ROOT:-${CHECK_OUT}/ref_free}
META_CSV=${META_CSV:-${INPUT_DIR}/meta.csv}

PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore
NEXTFLOW=${NEXTFLOW:-/appsnew/home/syfan/softwares/nextflow/nextflow}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}

TOTAL_REPEATS=${TOTAL_REPEATS:-10000}
REPEATS_PER_JOB=${REPEATS_PER_JOB:-200}
REF_N=${REF_N:-40}
SEED=${SEED:-42}
MIN_FF=${MIN_FF:-0.01}
CUTOFF=${CUTOFF:-3.0}
EZ_CUTOFF=${EZ_CUTOFF:-4.5}

run_py() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3 "$@"
}

echo "[$(date)] === baseline recall 0.6 / 0.65 compare ==="

# --- assets ---
ASSETS_06="${CHECK_OUT}/grid_search_assets_ep0.6"
ASSETS_065="${CHECK_OUT}/grid_search_assets_ep0.65"
run_py "${SCRIPT_DIR}/build_ep_recall_assets.py" \
    --source-assets "$SOURCE_ASSETS" \
    --output-dir "$ASSETS_06" \
    --ep-recall 0.6
run_py "${SCRIPT_DIR}/build_ep_recall_assets.py" \
    --source-assets "$SOURCE_ASSETS" \
    --output-dir "$ASSETS_065" \
    --ep-recall 0.65

# --- baseline-0.65: reuse existing pipeline/baseline ---
BL065_OUT="${PIPELINE_ROOT}/baseline-0.65"
mkdir -p "${BL065_OUT}/merge_scores"
if [ -d "${PIPELINE_ROOT}/baseline/merge_scores" ]; then
    echo "[$(date)] Linking baseline-0.65 merge_scores from pipeline/baseline"
    # copy to avoid broken relative symlinks across schemes
    cp -f "${PIPELINE_ROOT}/baseline/merge_scores/"*_scores.tsv "${BL065_OUT}/merge_scores/"
else
    echo "ERROR: missing ${PIPELINE_ROOT}/baseline/merge_scores" >&2
    exit 1
fi
n065=$(ls "${BL065_OUT}/merge_scores"/*_scores.tsv | wc -l)
echo "  baseline-0.65 scores: $n065"

# --- baseline-0.6: full scoring ---
BL06_OUT="${PIPELINE_ROOT}/baseline-0.6"
mkdir -p "$BL06_OUT" "$WORKDIR"
echo "[$(date)] Running aipt_ref_40 for baseline-0.6 (ep 0.5/0.6)"
"$NEXTFLOW" run "${PROJECT_DIR}/main.nf" \
    -profile aipt_ref_40,alioth_slurm,singularity \
    -work-dir "$WORKDIR" \
    -resume \
    --step aipt_ref_40 \
    --input "$SAMPLESHEET" \
    --outdir "$BL06_OUT" \
    --grid_search_result "$ASSETS_06" \
    --cpg_recall_dir "$RECALL_DIR" \
    --skip_ff true \
    2>&1 | tee "${BL06_OUT}/nextflow.log"
n06=$(ls "${BL06_OUT}/merge_scores"/*_scores.tsv | wc -l)
echo "  baseline-0.6 scores: $n06"

# --- ref_free for both (serial slices in this head job; 10k can be long) ---
# Prefer array submission for speed
export TOTAL_REPEATS REPEATS_PER_JOB REF_N SEED MIN_FF CUTOFF
export EZ_CUTOFF_MIN=3.0 EZ_CUTOFF_MAX=4.5 EZ_CUTOFF_STEP=0.1

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
ARRAY_LAST=$((N_JOBS - 1))
dep_ids=()
for scheme in baseline-0.6 baseline-0.65; do
    scheme_out="${REF_FREE_ROOT}/${scheme}"
    mkdir -p "$scheme_out"
    # clear previous slices for a clean aggregate
    rm -rf "${scheme_out}/ref_free_ezscore"
    jobid=$(sbatch --parsable \
        --job-name="blrf_${scheme}" \
        --array="0-${ARRAY_LAST}" \
        --export=ALL,SCHEME="${scheme}" \
        run_ref_free_scheme.slurm \
        "$META_CSV" "$PIPELINE_ROOT" "$scheme_out")
    echo "Submitted ref_free ${scheme} job_id=${jobid}"
    dep_ids+=("${jobid}")
done

dep_arg=$(IFS=:; echo "${dep_ids[*]}")
PLOT_DIR="${REF_FREE_ROOT}/plots_baseline_recall"
agg_job=$(sbatch --parsable \
    --dependency="afterok:${dep_arg}" \
    --export=ALL,EZ_CUTOFF="${EZ_CUTOFF}" \
    run_baseline_recall_aggregate.slurm \
    "$REF_FREE_ROOT")
echo "Submitted aggregate+plot job_id=${agg_job} dep=afterok:${dep_arg}"
echo "Plots → ${PLOT_DIR}"
echo "[$(date)] Head job done (pipeline finished; ref_free arrays running)"
