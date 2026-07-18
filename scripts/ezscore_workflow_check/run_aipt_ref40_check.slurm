#!/usr/bin/bash
#SBATCH -J ezscore_check
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# aipt_ref_40 vs meta ground-truth check (FF skipped).
#
# Steps:
#   1. Build original-combo assets (ep 0.5/0.65, z 0.85/0.95) from production refs
#   2. Prepare 10 check samples (incl. PTAY1183P)
#   3. Run aipt_ref_40 with --skip_ff
#   4. Compare merge_scores vs meta
#
# Usage:
#   sbatch run_aipt_ref40_check.slurm [<input_dir>] [<output_dir>]

set -euo pipefail

INPUT_DIR=${1:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260716-ezscore_check}
OUTPUT_DIR=${2:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260716-ezscore_check}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore
NEXTFLOW=${NEXTFLOW:-/appsnew/home/syfan/softwares/nextflow/nextflow}
SIF=${SIF:-/lustre1/cqyi/syfan/images/common_tools.sif}

ASSETS_DIR="${OUTPUT_DIR}/grid_search_assets"
PREPARED_DIR="${OUTPUT_DIR}/prepared"
PIPELINE_OUT="${OUTPUT_DIR}/pipeline"
COMPARISON_DIR="${OUTPUT_DIR}/comparison"
WORKDIR="${OUTPUT_DIR}/work_pipeline"

mkdir -p "$ASSETS_DIR" "$PREPARED_DIR" "$PIPELINE_OUT" "$COMPARISON_DIR" "$WORKDIR"

run_py() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3 "$@"
}

echo "[$(date)] Step 1: build original-combo assets"
run_py build_original_combo_assets.py --output-dir "$ASSETS_DIR"

echo "[$(date)] Step 2: prepare check samplesheet"
run_py prepare_check_samples.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$PREPARED_DIR"

echo "[$(date)] Step 3: aipt_ref_40 (skip_ff=true)"
"$NEXTFLOW" run "${PROJECT_DIR}/main.nf" \
    -profile aipt_ref_40,alioth_slurm,singularity \
    -work-dir "$WORKDIR" \
    -resume \
    --step aipt_ref_40 \
    --input "${PREPARED_DIR}/check_samplesheet.csv" \
    --outdir "$PIPELINE_OUT" \
    --grid_search_result "$ASSETS_DIR" \
    --skip_ff true \
    2>&1 | tee "${PIPELINE_OUT}/nextflow.log"

echo "[$(date)] Step 4: compare with meta"
run_py compare_with_meta.py \
    --input-dir "$INPUT_DIR" \
    --prepared-dir "$PREPARED_DIR" \
    --pipeline-outdir "$PIPELINE_OUT" \
    --output-dir "$COMPARISON_DIR"

echo "[$(date)] Done. Summary: ${COMPARISON_DIR}/check_comparison_summary.tsv"
