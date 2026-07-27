#!/usr/bin/bash
#SBATCH -J early_ref_score
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Score early_ref Normal samples under each blacklist/baseline recall dir
# so the ref-free pool can include all set=dev / Normal (96 → 40+40).
#
# Usage:
#   sbatch run_score_early_ref.slurm

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
mkdir -p logs

CHECK_OUT=${CHECK_OUT:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260718-cpg_black_list_check}
PIPELINE_ROOT=${PIPELINE_ROOT:-${CHECK_OUT}/pipeline}
ASSETS_DIR=${ASSETS_DIR:-${CHECK_OUT}/grid_search_assets}
RECALL_BUILD=${RECALL_BUILD:-${CHECK_OUT}/recall_build}
SAMPLESHEET=${SAMPLESHEET:-${CHECK_OUT}/prepared/early_ref_samplesheet.csv}
WORKDIR=${WORKDIR:-${CHECK_OUT}/work_pipeline_shared}
EARLY_OUT=${EARLY_OUT:-${CHECK_OUT}/pipeline_early_ref}

PROJECT_DIR=/lustre1/cqyi/AIPT_2.0/workflow/episcore
NEXTFLOW=${NEXTFLOW:-/appsnew/home/syfan/softwares/nextflow/nextflow}

SCHEMES=(baseline 15-site1-J 15-site1-M 20-site1-M 15-site2-J 15-site2-M 20-site2-M)

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: missing $SAMPLESHEET" >&2
    exit 1
fi

mkdir -p "$EARLY_OUT" "$WORKDIR"

echo "[$(date)] Scoring early_ref Normals for pool expansion"
echo "  samplesheet : $SAMPLESHEET"
echo "  shared work : $WORKDIR"

for scheme in "${SCHEMES[@]}"; do
    recall_dir="${RECALL_BUILD}/recall_dirs/${scheme}"
    scheme_out="${EARLY_OUT}/${scheme}"
    mkdir -p "$scheme_out"
    echo "[$(date)] scheme=${scheme}"
    "$NEXTFLOW" run "${PROJECT_DIR}/main.nf" \
        -profile aipt_ref_40,alioth_slurm,singularity \
        -work-dir "$WORKDIR" \
        -resume \
        --step aipt_ref_40 \
        --input "$SAMPLESHEET" \
        --outdir "$scheme_out" \
        --grid_search_result "$ASSETS_DIR" \
        --cpg_recall_dir "$recall_dir" \
        --skip_ff true \
        2>&1 | tee "${scheme_out}/nextflow.log"

    # Merge early_ref scores into the main scheme merge_scores dir
    dest="${PIPELINE_ROOT}/${scheme}/merge_scores"
    mkdir -p "$dest"
    src="${scheme_out}/merge_scores"
    if [ ! -d "$src" ]; then
        echo "ERROR: missing $src" >&2
        exit 1
    fi
    cp -f "${src}/"*_scores.tsv "$dest/"
    echo "[$(date)] copied $(ls "$src"/*_scores.tsv | wc -l) early_ref scores → $dest"
done

echo "[$(date)] Done scoring early_ref"
