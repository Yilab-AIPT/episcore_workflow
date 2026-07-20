#!/usr/bin/bash
#SBATCH -J cpg_bl_check
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# CpG blacklist evaluation (aipt_ref_40, fixed combo, FF skipped).
#
# Shared upstream: one Nextflow -work-dir so MethylDackel / depth filter /
# PREPARE_INPUTS are computed once. Each blacklist scheme only changes
# --cpg_recall_dir (recall 0.65 / 0.95 with blacklist sites removed) and
# re-runs CALC_EPISCORE_FLEXIBLE + CALC_ZSCORE_FLEXIBLE + MERGE_SCORES.
#
# Usage:
#   sbatch run_blacklist_check.slurm [<input_dir>] [<output_dir>]

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260718-cpg_black_list_check}
OUTPUT_DIR=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore
NEXTFLOW=${NEXTFLOW:-/appsnew/home/syfan/softwares/nextflow/nextflow}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}
EZSCORE_CHECK_DIR="${PROJECT_DIR}/scripts/ezscore_workflow_check"

ASSETS_DIR="${OUTPUT_DIR}/grid_search_assets"
PREPARED_DIR="${OUTPUT_DIR}/prepared"
RECALL_BUILD_DIR="${OUTPUT_DIR}/recall_build"
PIPELINE_ROOT="${OUTPUT_DIR}/pipeline"
COMPARISON_DIR="${OUTPUT_DIR}/comparison"
INTERSECT_DIR="${OUTPUT_DIR}/intersect"
# Shared across all schemes — MethylDackel etc. resume here
WORKDIR="${OUTPUT_DIR}/work_pipeline_shared"

mkdir -p "$ASSETS_DIR" "$PREPARED_DIR" "$RECALL_BUILD_DIR" \
    "$PIPELINE_ROOT" "$COMPARISON_DIR" "$INTERSECT_DIR" "$WORKDIR"

run_py() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3 "$@"
}

echo "[$(date)] Step 0: intersect blacklists with submitted / production panels"
run_py "${SCRIPT_DIR}/intersect_blacklists.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$INTERSECT_DIR"

echo "[$(date)] Step 1: build original-combo assets (ep 0.5/0.65, z 0.85/0.95)"
run_py "${EZSCORE_CHECK_DIR}/build_original_combo_assets.py" \
    --output-dir "$ASSETS_DIR"

echo "[$(date)] Step 2: prepare samplesheet (meta analyze ∩ mqres)"
run_py "${SCRIPT_DIR}/prepare_samplesheet.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$PREPARED_DIR"

echo "[$(date)] Step 3: build baseline + per-blacklist recall dirs"
run_py "${SCRIPT_DIR}/build_filtered_recall_dirs.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$RECALL_BUILD_DIR"

SCHEMES_FILE="${RECALL_BUILD_DIR}/schemes.txt"
mapfile -t SCHEMES < "$SCHEMES_FILE"
echo "[$(date)] Schemes (${#SCHEMES[@]}): ${SCHEMES[*]}"
echo "[$(date)] Shared workDir: $WORKDIR"

for scheme in "${SCHEMES[@]}"; do
    recall_dir="${RECALL_BUILD_DIR}/recall_dirs/${scheme}"
    scheme_out="${PIPELINE_ROOT}/${scheme}"
    mkdir -p "$scheme_out"
    echo "[$(date)] Step 4: aipt_ref_40 scheme=${scheme}"
    echo "  cpg_recall_dir=$recall_dir"
    echo "  outdir=$scheme_out"
    "$NEXTFLOW" run "${PROJECT_DIR}/main.nf" \
        -profile aipt_ref_40,alioth_slurm,singularity \
        -work-dir "$WORKDIR" \
        -resume \
        --step aipt_ref_40 \
        --input "${PREPARED_DIR}/samplesheet.csv" \
        --outdir "$scheme_out" \
        --grid_search_result "$ASSETS_DIR" \
        --cpg_recall_dir "$recall_dir" \
        --skip_ff true \
        2>&1 | tee "${scheme_out}/nextflow.log"
done

echo "[$(date)] Step 5: compare all schemes vs meta (+ vs baseline)"
run_py "${SCRIPT_DIR}/compare_scores.py" \
    --input-dir "$INPUT_DIR" \
    --prepared-dir "$PREPARED_DIR" \
    --pipeline-root "$PIPELINE_ROOT" \
    --schemes-file "$SCHEMES_FILE" \
    --output-dir "$COMPARISON_DIR"

echo "[$(date)] Done."
echo "  Intersect : ${INTERSECT_DIR}/blacklist_source_report.tsv"
echo "  Rollup    : ${COMPARISON_DIR}/scheme_rollup.tsv"
echo "  Notebook  : ${PROJECT_DIR}/notebooks/aipt_2.0/cpg_black_list_check.ipynb"
