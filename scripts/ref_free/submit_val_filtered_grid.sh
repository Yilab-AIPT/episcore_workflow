#!/usr/bin/bash
# Build filtered-range episcore/zscore grids for val samples (blacklist applied),
# then (after current filtered 100k finishes) merge + re-run filtered with val
# and regenerate dual plots.
#
# Usage:
#   FILT_ARRAY_JOB=<id> ./submit_val_filtered_grid.sh
#   # or auto-detect running ref_free_filt array

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"
mkdir -p logs

TODAY_BASE=${TODAY_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260721-ref_free}
VAL_META=${VAL_META:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260702-ref_40_20260625_samples}
MAIN_INPUT=${MAIN_INPUT:-/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng}
BETA_ROOT=${BETA_ROOT:-/lustre1/cqyi/myli/bert/DNA_5mC_analysis_pipeline/output/20260625-XML/bwameth_results/zscore_downstream/beta_zscore}
SIF=${SIF:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif}
NF=${NF:-/appsnew/home/syfan/softwares/nextflow/nextflow}
PIPE=${PIPE:-/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf}

VAL_GRID="${TODAY_BASE}/val_filtered_grid"
Z_WORK="${VAL_GRID}/zscore"
EP_WORK="${VAL_GRID}/episcore"
ANALYZE_05="${EP_WORK}/analyze_thr0.5"
B2Z_OUT="${EP_WORK}/b2z_output"
MERGED_INPUT="${TODAY_BASE}/input_with_val_filtered"
FILTERED_VAL_BASE="${TODAY_BASE}/filtered_combos_with_val"

mkdir -p "$Z_WORK" "$EP_WORK" "$ANALYZE_05" "$B2Z_OUT" logs

# --- detect filt array dependency ---
FILT_ARRAY_JOB=${FILT_ARRAY_JOB:-}
if [[ -z "$FILT_ARRAY_JOB" ]]; then
    FILT_ARRAY_JOB=$(squeue -u "$USER" -n ref_free_filt -h -o '%i' | head -1 | cut -d_ -f1 || true)
fi
echo "Filt array dependency : ${FILT_ARRAY_JOB:-none}"
echo "Val grid root         : $VAL_GRID"

# --- 1) prepare zscore tasks + samplesheet ---
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
    python3 prepare_val_filtered_zscore_tasks.py \
        --val-meta-dir "$VAL_META" \
        --output-dir "$Z_WORK"

N_Z=$(awk 'NR>1' "$Z_WORK/zscore_tasks.tsv" | wc -l)
echo "Submitting val zscore (${N_Z} tasks; MaxArraySize=1001 → chunked)"
# Chunk 0: indices 0-999
z_job1=$(sbatch --parsable --array="0-999%50" \
    run_val_zscore.slurm "$Z_WORK/zscore_tasks.tsv" "$Z_WORK")
# Chunk 1: logical 1000..(N_Z-1) via TASK_OFFSET
Z_REM=$((N_Z - 1000))
if (( Z_REM > 0 )); then
    Z_LAST2=$((Z_REM - 1))
    z_job2=$(TASK_OFFSET=1000 sbatch --parsable --export=ALL,TASK_OFFSET=1000 \
        --array="0-${Z_LAST2}%50" \
        run_val_zscore.slurm "$Z_WORK/zscore_tasks.tsv" "$Z_WORK")
else
    z_job2=""
fi
z_job="${z_job1}${z_job2:+:${z_job2}}"
echo "  zscore jobs=$z_job"

# --- 2) symlink thr=0.5 betas (blacklist excluded via samplesheet samples) ---
mapfile -t VAL_SAMPLES < <(awk -F',' 'NR>1 {print $1}' "$Z_WORK/samplesheet.csv" | sort -u)
for s in "${VAL_SAMPLES[@]}"; do
    src="${BETA_ROOT}/${s}/extract_beta_value/${s}_beta_value.tsv.gz"
    dst="${ANALYZE_05}/${s}_beta_value.tsv.gz"
    if [[ -f "$src" ]]; then
        ln -sfn "$src" "$dst"
    else
        echo "WARN missing beta $src" >&2
    fi
done
echo "Linked $(ls "$ANALYZE_05"/*_beta_value.tsv.gz 2>/dev/null | wc -l) thr=0.5 betas"

# b2z tasks: thr 0.5 × recall 0.5..0.75
B2Z_TASKS="${EP_WORK}/b2z_tasks_thr0.5.tsv"
: > "$B2Z_TASKS"
for rec in $(awk 'BEGIN { for (i = 50; i <= 75; i++) printf "%g\n", i/100 }'); do
    printf '0.5\t%s\n' "$rec" >> "$B2Z_TASKS"
done
N_B2Z=$(wc -l < "$B2Z_TASKS")
B2Z_LAST=$((N_B2Z - 1))
echo "Submitting val b2z thr0.5 array 0-${B2Z_LAST}"
b2z05_job=$(sbatch --parsable --array="0-${B2Z_LAST}" \
    run_val_beta_to_episcore.slurm "$B2Z_TASKS" "$ANALYZE_05" "$B2Z_OUT")
echo "  b2z0.5 job=$b2z05_job"

# --- 3) nextflow extract_beta at thr 0.1 and 0.33 ---
# Build mqres without blacklist
MQRES_FILT="${VAL_GRID}/mqres_noblacklist.csv"
singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3 - <<PY
import pandas as pd
from pathlib import Path
from val_blacklist import VAL_BLACKLIST, drop_blacklisted
meta = pd.read_csv("${VAL_META}/meta.csv")
meta = meta[meta["label"].astype(str).eq("Normal") | meta["label"].astype(str).str.match(r"^T\\d")]
meta = drop_blacklisted(meta)
keep = set(meta["sample"].astype(str))
mq = pd.read_csv("${VAL_META}/mqres.csv")
mq = mq[mq["sample"].astype(str).isin(keep)]
mq.to_csv("${MQRES_FILT}", index=False)
print("mqres rows", len(mq), "samples", mq["sample"].nunique())
PY

submit_nf_thr() {
    local thr=$1
    local outdir="${EP_WORK}/nf_extract_thr${thr}"
    mkdir -p "$outdir"
    # launcher slurm so NF can submit children
    sbatch --parsable --job-name="val_ext_${thr}" \
        --partition=cn-long --cpus-per-task=4 --mem=8G --time=24:00:00 \
        --output="logs/val_extract_${thr}_%j.log" \
        --wrap="cd '${SCRIPT_DIR}' && \
            ${NF} run '${PIPE}' \
              -profile grid_search,alioth_slurm,singularity \
              --step grid_search \
              --input '${MQRES_FILT}' \
              --outdir '${outdir}' \
              --threshold ${thr} \
              -resume"
}

echo "Submitting nextflow extract for thr=0.1 and 0.33"
nf01=$(submit_nf_thr 0.1)
nf033=$(submit_nf_thr 0.33)
echo "  nf 0.1=$nf01  nf 0.33=$nf033"

# After extract: collect betas into analyze dirs and run b2z for thr 0.1 / 0.33
collect_and_b2z() {
    local thr=$1
    local nf_job=$2
    local analyze="${EP_WORK}/analyze_thr${thr}"
    local tasks="${EP_WORK}/b2z_tasks_thr${thr}.tsv"
    mkdir -p "$analyze"
    : > "$tasks"
    for rec in $(awk 'BEGIN { for (i = 50; i <= 75; i++) printf "%g\n", i/100 }'); do
        printf '%s\t%s\n' "$thr" "$rec" >> "$tasks"
    done
    local n
    n=$(wc -l < "$tasks")
    local last=$((n - 1))
    # collector job depends on nf
    local collect
    collect=$(sbatch --parsable --job-name="val_collect_${thr}" \
        --dependency="afterok:${nf_job}" \
        --partition=cn-long --cpus-per-task=2 --mem=4G --time=01:00:00 \
        --output="logs/val_collect_${thr}_%j.log" \
        --wrap="
set -euo pipefail
ANALYZE='${analyze}'
NFDIR='${EP_WORK}/nf_extract_thr${thr}'
mkdir -p \"\$ANALYZE\"
# betas published under extract_beta_value/
find \"\$NFDIR\" -name '*_beta_value.tsv.gz' | while read -r f; do
  bn=\$(basename \"\$f\")
  ln -sfn \"\$f\" \"\$ANALYZE/\$bn\"
done
echo linked \$(ls \"\$ANALYZE\"/*_beta_value.tsv.gz 2>/dev/null | wc -l) betas for thr=${thr}
")
    local b2z
    b2z=$(sbatch --parsable --job-name="val_b2z_${thr}" \
        --dependency="afterok:${collect}" \
        --array="0-${last}" \
        run_val_beta_to_episcore.slurm "$tasks" "$analyze" "$B2Z_OUT")
    echo "$b2z"
}

b2z01=$(collect_and_b2z 0.1 "$nf01")
b2z033=$(collect_and_b2z 0.33 "$nf033")
echo "  b2z 0.1=$b2z01  b2z 0.33=$b2z033"

# --- 4) after all grids + filt: merge, run filtered+val, plot dual ---
DEP_GRID="afterok:${z_job}:${b2z05_job}:${b2z01}:${b2z033}"
if [[ -n "${FILT_ARRAY_JOB}" ]]; then
    DEP_ALL="${DEP_GRID}:${FILT_ARRAY_JOB}"
else
    DEP_ALL="$DEP_GRID"
fi

merge_run=$(sbatch --parsable --job-name=val_merge_run \
    --dependency="$DEP_ALL" \
    --partition=cn-long --cpus-per-task=8 --mem=64G --time=48:00:00 \
    --output=logs/val_merge_run_%j.log \
    --wrap="
set -euo pipefail
cd '${SCRIPT_DIR}'
SIF='${SIF}'
singularity exec -B /lustre1,/lustre2,/appsnew \"\$SIF\" \
  python3 merge_val_filtered_parquets.py \
    --main-input '${MAIN_INPUT}' \
    --val-meta-dir '${VAL_META}' \
    --val-zscore-root '${Z_WORK}' \
    --val-episcore-root '${B2Z_OUT}' \
    --output-dir '${MERGED_INPUT}'

# Launch filtered+val 100k (reuse env from submit_ref_free_today)
export TOTAL_REPEATS=100000 REPEATS_PER_JOB=2000 REF_N=40 SEED=42
export COMBO_MODE=all STORE_PAIR_COUNTS=0 COMPRESS=1
export EZ_CUTOFF_MIN=3.0 EZ_CUTOFF_MAX=4.5 EZ_CUTOFF_STEP=0.1
export CUTOFF=3.0 MIN_FF=0 FF_MIN=0.01
export EP_THRESHOLD_MIN=0.1 EP_THRESHOLD_MAX=0.5
export EP_RECALL_MIN=0.5 EP_RECALL_MAX=0.75
export Z_THRESHOLD_MIN=0.8 Z_THRESHOLD_MAX=0.95
export Z_RECALL_MIN=0.9 Z_RECALL_MAX=0.99
unset EP_THRESHOLD EP_RECALL Z_THRESHOLD Z_RECALL

mkdir -p '${FILTERED_VAL_BASE}/plots'
rm -rf '${FILTERED_VAL_BASE}/ref_free_ezscore'
arr=\$(sbatch --parsable --job-name=ref_free_filt_val \
  --array=0-49 \
  run_ref_free_ezscore.slurm '${MERGED_INPUT}' '${FILTERED_VAL_BASE}')
echo filt_val_array=\$arr
plot=\$(sbatch --parsable --job-name=plot_filt_val \
  --dependency=afterok:\${arr} \
  run_aggregate_and_plot.slurm \
  '${FILTERED_VAL_BASE}' \
  '40+40 filtered combos + val' \
  0.01 0)
echo plot_filt_val=\$plot
")
echo "Submitted merge+filtered_val chain job_id=${merge_run} dep=${DEP_ALL}"

# --- 5) immediately replot fixed with blacklist + new sep metrics ---
if [[ -f "${TODAY_BASE}/fixed_combo/ref_free_ezscore/abnormality_signal_ratio.tsv" ]]; then
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
        python3 aggregate_ref_free_ezscore.py --output-base "${TODAY_BASE}/fixed_combo" --ff-min 0.01 || true
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
        python3 plot_ref_free_interactive.py \
            --result-dir "${TODAY_BASE}/fixed_combo" \
            --title "40+40 fixed combo (ep 0.5/0.65, z 0.85/0.95)" \
            --ff-min 0.01
    singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" \
        python3 plot_ref_free_dual_ezscore.py \
            --result-dir "${TODAY_BASE}/fixed_combo" \
            --title "40+40 fixed combo" \
            --ff-min 0.01
    echo "Replotted fixed dual + 3-panel"
fi

echo "Done submitting val filtered grid pipeline."
