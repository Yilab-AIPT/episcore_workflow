#!/usr/bin/bash
#SBATCH -J wgs_ff
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run WGS pileup + FF estimation for one BAM.
#
# Usage:
#     sbatch --job-name=wgs_ff_<bam_id> run_sample.slurm <bam_id> <input_bam> <output_dir>
#
# Outputs:
#     <output_dir>/<bam_id>_pileup.tsv.gz
#     <output_dir>/<bam_id>_ff.tsv

set -euo pipefail

bam_id=${1:?Usage: sbatch run_sample.slurm <bam_id> <input_bam> <output_dir>}
input_bam=${2:?Usage: sbatch run_sample.slurm <bam_id> <input_bam> <output_dir>}
OUTPUT_DIR=${3:?Usage: sbatch run_sample.slurm <bam_id> <input_bam> <output_dir>}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$SCRIPT_DIR"

KNOWN_SITES=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-est_ff_for_wgs/whole_genome_wide_af_filtered_ChinaMAP.tsv
CONTAINER="${REPO_ROOT}/containers/common_tools.sif"

mkdir -p logs
mkdir -p "$OUTPUT_DIR"

if [ ! -f "$input_bam" ]; then
    echo "ERROR: BAM not found: $input_bam" >&2
    exit 1
fi
if [ ! -f "$KNOWN_SITES" ]; then
    echo "ERROR: known sites not found: $KNOWN_SITES" >&2
    exit 1
fi
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: container not found: $CONTAINER" >&2
    exit 1
fi

NCPUS="${SLURM_CPUS_PER_TASK:-32}"
OUT_PREFIX="${OUTPUT_DIR}/${bam_id}"
PILEUP="${OUT_PREFIX}_pileup.tsv.gz"
FF_OUT="${OUT_PREFIX}_ff.tsv"

echo "[$(date)] bam_id=${bam_id} input_bam=${input_bam} output_dir=${OUTPUT_DIR} cpus=${NCPUS} host=$(hostname)"

run_py() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" python3 "$@"
}

if [ -f "$FF_OUT" ]; then
    echo "SKIP: output already exists: $FF_OUT"
    exit 0
fi

if [ ! -f "$PILEUP" ]; then
    echo "[$(date)] Running bam_to_pileup_wgs..."
    run_py bam_to_pileup_wgs.py \
        --input-bam    "$input_bam" \
        --known-sites  "$KNOWN_SITES" \
        --output       "$OUT_PREFIX" \
        --ncpus        "$NCPUS"
else
    echo "[$(date)] Reusing existing pileup: $PILEUP"
fi

echo "[$(date)] Running estimate_ff_wgs..."
run_py estimate_ff_wgs.py \
    --input-path    "$PILEUP" \
    --output-prefix "$OUT_PREFIX" \
    --ff-min        0 \
    --ff-max        0.3 \
    --ff-precision  0.0001 \
    --min-depth     3 \
    --mode-list     "all,chr_only,chr_exclude" \
    --ncpus         "$NCPUS"

echo "[$(date)] Finished ${bam_id}"
