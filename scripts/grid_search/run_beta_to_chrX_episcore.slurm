#!/usr/bin/bash
#SBATCH -J b2z_recall
#SBATCH --partition=cn-long
#SBATCH --cpus-per-task=32
#SBATCH -o logs/%x.log
#SBATCH -e logs/%x.log

# Run beta_to_episcore.py for one recall threshold.
#
# Usage:
#     sbatch --job-name=b2z_recall_<recall> run_beta_to_episcore.slurm <recall>
#
# The job name is also used to name the log file (logs/%x.log), so be sure to
# pass --job-name when submitting (submit.sh handles that automatically).

set -euo pipefail

recall=${1:?Usage: sbatch run_beta_to_episcore.slurm <recall>}

WORKDIR=/lustre1/cqyi/AIPT_2.0/tests/episcore/grid_search
cd "$WORKDIR"

mkdir -p logs output_chrX

REF_DIR=female_reference_beta_values
ANALYZE_DIR=analyze_beta_values
CPG_LIST=cpg_level_metamodel_240k/CpG_final_filtered_recall${recall}_allsite.txt

if [ ! -f "$CPG_LIST" ]; then
    echo "ERROR: CpG list not found: $CPG_LIST" >&2
    exit 1
fi

echo "[$(date)] Starting recall=${recall} on host $(hostname)"
echo "  reference dir : $REF_DIR"
echo "  analyze dir   : $ANALYZE_DIR"
echo "  cpg list      : $CPG_LIST"
echo "  cpus          : ${SLURM_CPUS_PER_TASK:-32}"

singularity exec -B /lustre1,/lustre2,/appsnew \
    /lustre1/cqyi/syfan/images/common_tools.sif \
    python3 beta_to_episcore.py \
        --reference-beta-value-dir "$REF_DIR" \
        --analyze-beta-value-dir   "$ANALYZE_DIR" \
        --output-prefix            "output_chrX/recall_${recall}" \
        --depth 30 \
        --cpg-list                 "$CPG_LIST" \
        --chr-list '1-22,X' \
        --ncpus "${SLURM_CPUS_PER_TASK:-32}"

echo "[$(date)] Finished recall=${recall}"
