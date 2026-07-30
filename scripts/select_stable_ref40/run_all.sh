#!/usr/bin/env bash
# End-to-end: build source tables -> select ref_40 -> write updated meta.
set -euo pipefail

ROOT="/lustre1/cqyi/AIPT_2.0/workflow/episcore"
SIF="${ROOT}/containers/common_tools.sif"
SCRIPT_DIR="${ROOT}/scripts/select_stable_ref40"
OUT_DIR="/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260730-stable_ref40"

EPISCORE_SS="/lustre1/cqyi/syfan/nipt_article_plot/episcore_result_samplesheet.csv"
PCT_CSV="/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40/percentage.csv"
META_CSV="/lustre1/cqyi/syfan/nipt_article_plot/temporary_updated_samplesheet.csv"
EZ_REF="/lustre1/cqyi/myli/bert/analysis_nipt/multiomics/chr_stats_reference_samples.txt"

N_RANDOM="${N_RANDOM:-80000}"
N_SWAP="${N_SWAP:-15000}"
SEED="${SEED:-123}"
# Empty exclude keeps all Normal+dev eligible (incl. historical early_ref members)
EXCLUDE_SAMPLE="${EXCLUDE_SAMPLE:-}"

mkdir -p "${OUT_DIR}"

run_py() {
  singularity exec \
    -B /lustre1/cqyi:/lustre1/cqyi \
    -B /appsnew/home/myli:/appsnew/home/myli \
    -B /appsnew/home/syfan:/appsnew/home/syfan \
    "${SIF}" python3 "$@"
}

echo "=== 1/4 build_source_tables ==="
run_py "${SCRIPT_DIR}/build_source_tables.py" \
  --episcore-samplesheet "${EPISCORE_SS}" \
  --percentage-csv "${PCT_CSV}" \
  --meta-csv "${META_CSV}" \
  --output-dir "${OUT_DIR}"

# Provide meta.csv copy alongside beta/percentage for select_ref40 defaults
cp -f "${META_CSV}" "${OUT_DIR}/meta.csv"

echo "=== 2/4 select_ref40 ==="
run_py "${SCRIPT_DIR}/select_ref40.py" \
  --input-dir "${OUT_DIR}" \
  --meta-csv "${OUT_DIR}/meta.csv" \
  --output-dir "${OUT_DIR}" \
  --n-random "${N_RANDOM}" \
  --n-swap-rounds "${N_SWAP}" \
  --seed "${SEED}" \
  --exclude-sample "${EXCLUDE_SAMPLE}" \
  --ezscore-ref-samples "${EZ_REF}"

echo "=== 3/4 write_updated_meta ==="
run_py "${SCRIPT_DIR}/write_updated_meta.py" \
  --meta-csv "${META_CSV}" \
  --ref40-samples "${OUT_DIR}/ref40_samples.txt" \
  --score-tsv "${OUT_DIR}/ref40_score.tsv" \
  --output-csv "${OUT_DIR}/temporary_updated_samplesheet_ref40.csv"

cp -f "${OUT_DIR}/temporary_updated_samplesheet_ref40.csv" \
  /lustre1/cqyi/syfan/nipt_article_plot/temporary_updated_samplesheet_ref40.csv

echo "=== 4/4 render_plots ==="
run_py "${SCRIPT_DIR}/render_plots.py"

echo "Done. Outputs in ${OUT_DIR}"
cat "${OUT_DIR}/selection_summary.json"
echo
cat "${OUT_DIR}/zscore_fetch_summary.json"
