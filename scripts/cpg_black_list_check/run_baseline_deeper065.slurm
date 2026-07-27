#!/usr/bin/bash
#SBATCH -J bl_deeper065
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Score baseline-deeper-0.65 (ep combo 0.5/0.65,
# CpG list = replaced_deeper_recall_0.65 — same n as prod 0.65, ~15% sites swapped),
# then ref_free + 3-way baseline plot vs baseline-0.6 / baseline-0.65.
#
# Reuses work_pipeline_shared so SPLIT_BAM / MethylDackel / depth filter resume
# from prior baseline runs; only scoring sees the swapped 0.65 recall file.
#
# Usage:
#   sbatch run_baseline_deeper065.slurm

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

CHECK_OUT=${CHECK_OUT:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}
INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260718-cpg_black_list_check}
PIPELINE_ROOT=${PIPELINE_ROOT:-${CHECK_OUT}/pipeline}
RECALL_BUILD=${RECALL_BUILD:-${CHECK_OUT}/recall_build}
SAMPLESHEET=${SAMPLESHEET:-${CHECK_OUT}/prepared/samplesheet_all.csv}
WORKDIR=${WORKDIR:-${CHECK_OUT}/work_pipeline_shared}
SOURCE_ASSETS=${SOURCE_ASSETS:-${CHECK_OUT}/grid_search_assets}
ASSETS_065=${ASSETS_065:-${CHECK_OUT}/grid_search_assets_ep0.65}
REF_FREE_ROOT=${REF_FREE_ROOT:-${CHECK_OUT}/ref_free}
META_CSV=${META_CSV:-${INPUT_DIR}/meta.csv}
SCHEME=${SCHEME:-baseline-deeper-0.65}

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

echo "[$(date)] === ${SCHEME} (deeper 0.65 CpG list) ==="

# Ensure ep-0.65 assets exist (same combo as production baseline-0.65)
if [ ! -f "${ASSETS_065}/best_combo_episcore.csv" ]; then
    run_py "${SCRIPT_DIR}/build_ep_recall_assets.py" \
        --source-assets "$SOURCE_ASSETS" \
        --output-dir "$ASSETS_065" \
        --ep-recall 0.65
fi

echo "[$(date)] Build deeper-0.65 recall dir"
run_py "${SCRIPT_DIR}/build_deeper_recall_dir.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$RECALL_BUILD" \
    --scheme-name "$SCHEME"

RECALL_DIR="${RECALL_BUILD}/recall_dirs/${SCHEME}"
SCHEME_OUT="${PIPELINE_ROOT}/${SCHEME}"
mkdir -p "$SCHEME_OUT" "$WORKDIR"

echo "[$(date)] aipt_ref_40 scheme=${SCHEME}"
echo "  cpg_recall_dir=${RECALL_DIR}  (0.65 = deeper sites; rest symlinked)"
echo "  grid_search_result=${ASSETS_065}"
echo "  workDir=${WORKDIR} (-resume; upstream shared with baseline-0.6/0.65)"
"$NEXTFLOW" run "${PROJECT_DIR}/main.nf" \
    -profile aipt_ref_40,alioth_slurm,singularity \
    -work-dir "$WORKDIR" \
    -resume \
    --step aipt_ref_40 \
    --input "$SAMPLESHEET" \
    --outdir "$SCHEME_OUT" \
    --grid_search_result "$ASSETS_065" \
    --cpg_recall_dir "$RECALL_DIR" \
    --skip_ff true \
    2>&1 | tee "${SCHEME_OUT}/nextflow.log"

n=$(ls "${SCHEME_OUT}/merge_scores"/*_scores.tsv | wc -l)
echo "  ${SCHEME} scores: $n"

# --- ref_free only for the new scheme (0.6 / 0.65 already done) ---
export TOTAL_REPEATS REPEATS_PER_JOB REF_N SEED MIN_FF CUTOFF
export EZ_CUTOFF_MIN=3.0 EZ_CUTOFF_MAX=4.5 EZ_CUTOFF_STEP=0.1

N_JOBS=$(( (TOTAL_REPEATS + REPEATS_PER_JOB - 1) / REPEATS_PER_JOB ))
ARRAY_LAST=$((N_JOBS - 1))
scheme_rf="${REF_FREE_ROOT}/${SCHEME}"
mkdir -p "$scheme_rf"
rm -rf "${scheme_rf}/ref_free_ezscore"
jobid=$(sbatch --parsable \
    --job-name="blrf_${SCHEME}" \
    --array="0-${ARRAY_LAST}" \
    --export=ALL,SCHEME="${SCHEME}" \
    run_ref_free_scheme.slurm \
    "$META_CSV" "$PIPELINE_ROOT" "$scheme_rf")
echo "Submitted ref_free ${SCHEME} job_id=${jobid}"

agg_job=$(sbatch --parsable \
    --dependency="afterok:${jobid}" \
    --export=ALL,EZ_CUTOFF="${EZ_CUTOFF}" \
    run_baseline_recall_aggregate.slurm \
    "$REF_FREE_ROOT")
echo "Submitted aggregate+plot job_id=${agg_job}"
echo "Plots → ${REF_FREE_ROOT}/plots_baseline_recall"
echo "[$(date)] Head job done (pipeline finished; ref_free array running)"
