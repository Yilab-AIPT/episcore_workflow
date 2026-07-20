#!/usr/bin/bash
# Check grid parquet coverage, submit fill jobs for missing zscore cells, then patch parquet.
#
# Usage:
#   ./submit_fill_missing_zscore.sh [-n|--dry-run]
#
# Env overrides:
#   INPUT_DIR, ZSCORE_ROOT, PARQUET_PATH, SIF

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

DRY_RUN=false
if [[ "${1:-}" == "-n" || "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

INPUT_DIR=${INPUT_DIR:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
ZSCORE_ROOT=${ZSCORE_ROOT:-/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search}
PARQUET_PATH=${PARQUET_PATH:-${INPUT_DIR}/zscore_grid_search.parquet}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
COVER_DIR=${COVER_DIR:-${INPUT_DIR}/coverage_check}

echo "=== Coverage audit ==="
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 check_grid_coverage.py \
        --input-dir "$INPUT_DIR" \
        --output-dir "$COVER_DIR"

MISSING_TSV="${COVER_DIR}/missing_coverage.tsv"
N_EP=$(awk -F'\t' 'NR>1 && $1=="episcore" {c++} END{print c+0}' "$MISSING_TSV")
N_Z=$(awk -F'\t' 'NR>1 && $1=="zscore" {c++} END{print c+0}' "$MISSING_TSV")

echo "  episcore missing cells : $N_EP"
echo "  zscore   missing cells : $N_Z"

if (( N_EP > 0 )); then
    echo "WARN: episcore gaps detected. Fill is not automated here; inspect ${MISSING_TSV}" >&2
fi

# Drop zscore rows whose CSV already exists on disk (resume-safe).
REMAINING_TSV="${COVER_DIR}/missing_coverage_remaining.tsv"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 - <<PY
from pathlib import Path
import pandas as pd
from grid_coverage import fmt_float

miss = pd.read_csv("${MISSING_TSV}", sep="\\t")
zroot = Path("${ZSCORE_ROOT}") / "zscore_results"
rows = []
for _, r in miss.iterrows():
    if str(r["score_type"]) != "zscore":
        rows.append(r)
        continue
    sample = str(r["sample"])
    thr, rec = fmt_float(r["threshold"]), fmt_float(r["recall"])
    path = zroot / f"recall.{rec}_cutoff.{thr}" / f"{sample}.{thr}.1.0.220k_cpg_recall_{rec}.NoLen.zscore.csv"
    if path.is_file() and sum(1 for _ in open(path)) >= 23:
        continue
    rows.append(r)
out = pd.DataFrame(rows)
if not out.empty:
    out = out[["score_type", "sample", "threshold", "recall"]]
else:
    out = pd.DataFrame(columns=["score_type", "sample", "threshold", "recall"])
out.to_csv("${REMAINING_TSV}", sep="\\t", index=False)
print(f"  remaining after disk check: {len(out)}")
PY

N_Z=$(awk -F'\t' 'NR>1 && $1=="zscore" {c++} END{print c+0}' "$REMAINING_TSV")
if (( N_Z == 0 )); then
    echo "No remaining zscore gaps on disk; submitting patch only."
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] sbatch run_patch_zscore_parquet.slurm ${MISSING_TSV} ..."
        exit 0
    fi
    patch_jobid=$(sbatch --parsable \
        --job-name=patch_zscore_pq \
        run_patch_zscore_parquet.slurm \
        "$MISSING_TSV" "$PARQUET_PATH" "${ZSCORE_ROOT}/zscore_results")
    echo "Submitted patch job_id=${patch_jobid}"
    exit 0
fi

ARRAY_LAST=$((N_Z - 1))
echo "=== Submit fill array 0-${ARRAY_LAST} (remaining only) ==="

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] sbatch --array=0-${ARRAY_LAST} run_fill_missing_zscore.slurm \\"
    echo "           ${REMAINING_TSV} ${ZSCORE_ROOT}"
    echo "[DRY RUN] sbatch --dependency=afterok:<fill_jobid> run_patch_zscore_parquet.slurm \\"
    echo "           ${MISSING_TSV} ${PARQUET_PATH} ${ZSCORE_ROOT}/zscore_results"
    exit 0
fi

fill_jobid=$(sbatch --parsable \
    --array=0-"${ARRAY_LAST}" \
    --job-name=fill_zscore \
    run_fill_missing_zscore.slurm \
    "$REMAINING_TSV" "$ZSCORE_ROOT")
echo "Submitted fill array job_id=${fill_jobid}  logs=logs/fill_missing_zscore_*.log"

# Patch uses the original missing list (all cells that need to land in parquet).
patch_jobid=$(sbatch --parsable \
    --dependency=afterok:"${fill_jobid}" \
    --job-name=patch_zscore_pq \
    run_patch_zscore_parquet.slurm \
    "$MISSING_TSV" "$PARQUET_PATH" "${ZSCORE_ROOT}/zscore_results")
echo "Submitted patch job_id=${patch_jobid} (afterok:${fill_jobid})  logs=logs/patch_zscore_parquet_*.log"
echo "Cover dir: ${COVER_DIR}"
echo "Remaining tsv: ${REMAINING_TSV}"
