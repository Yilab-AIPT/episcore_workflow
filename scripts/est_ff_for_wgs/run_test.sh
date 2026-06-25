#!/usr/bin/bash
# Quick end-to-end test: subsample 10,000 reads from one BAM, then run pileup + FF.
#
# Usage:
#     ./run_test.sh [input_bam]
#
# Default BAM: first row of the WGS samplesheet.
# Uses --max-sites 50000 on the known-sites panel so the test finishes in reasonable time.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$SCRIPT_DIR"

SAMPLESHEET=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-est_ff_for_wgs/samplesheet.csv
KNOWN_SITES=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-est_ff_for_wgs/whole_genome_wide_af_filtered_ChinaMAP.tsv
OUTPUT_DIR=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260621-est_ff_for_wgs/test_10k_reads
CONTAINER="${REPO_ROOT}/containers/common_tools.sif"

N_READS=10000
MAX_SITES=50000
NCPUS=8

if [ -n "${1:-}" ]; then
    INPUT_BAM="$1"
else
    INPUT_BAM=$(tail -n +2 "$SAMPLESHEET" | head -n 1 | cut -d',' -f2 | tr -d '"\r' | xargs)
fi

bam_id="test_${N_READS}reads_$(basename "$INPUT_BAM" .bam)"
TEST_DIR="${OUTPUT_DIR}/${bam_id}"
SUB_BAM="${TEST_DIR}/${bam_id}.bam"
OUT_PREFIX="${TEST_DIR}/${bam_id}"

mkdir -p "$TEST_DIR"

echo "Input BAM   : $INPUT_BAM"
echo "Subsample   : ${N_READS} reads -> $SUB_BAM"
echo "Known sites : $KNOWN_SITES (max ${MAX_SITES} for test)"
echo "Output      : $OUT_PREFIX"
echo "Container   : $CONTAINER"
echo

if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: BAM not found: $INPUT_BAM" >&2
    exit 1
fi
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: container not found: $CONTAINER" >&2
    exit 1
fi

run_py() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" python3 "$@"
}

run_samtools() {
    singularity exec -B /lustre1,/lustre2,/appsnew "$CONTAINER" samtools "$@"
}

if [ ! -f "$SUB_BAM" ]; then
    echo "[$(date)] Subsampling ${N_READS} reads (by read name, then sort)..."
    NAMES_FILE="${TEST_DIR}/read_names.txt"
    # Pick random QNAMEs; for paired-end BAMs each name yields two alignments.
    N_NAMES=$((N_READS / 2))
    run_samtools view "$INPUT_BAM" | cut -f1 | sort -u | shuf -n "$N_NAMES" > "$NAMES_FILE"
    run_samtools view -@ "$NCPUS" -b -N "$NAMES_FILE" "$INPUT_BAM" \
        | run_samtools sort -@ "$NCPUS" -o "$SUB_BAM" -
    run_samtools index "$SUB_BAM"
    echo "[$(date)] Subsampled BAM: $SUB_BAM ($(run_samtools view -c "$SUB_BAM") alignments)"
else
    echo "[$(date)] Reusing subsampled BAM: $SUB_BAM"
fi

echo "[$(date)] Running bam_to_pileup_wgs..."
run_py bam_to_pileup_wgs.py \
    --input-bam   "$SUB_BAM" \
    --known-sites "$KNOWN_SITES" \
    --output      "$OUT_PREFIX" \
    --max-sites   "$MAX_SITES" \
    --ncpus       "$NCPUS"

PILEUP="${OUT_PREFIX}_pileup.tsv.gz"
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

echo
echo "[$(date)] Test complete."
echo "  Pileup: ${PILEUP}"
echo "  FF:     ${OUT_PREFIX}_ff.tsv"
if [ -f "${OUT_PREFIX}_ff.tsv" ]; then
    cat "${OUT_PREFIX}_ff.tsv"
fi
