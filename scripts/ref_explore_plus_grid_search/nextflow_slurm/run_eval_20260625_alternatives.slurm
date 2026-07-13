#!/usr/bin/bash
#SBATCH -J ref40_eval
#SBATCH --partition=cn-long
#SBATCH -o ref40_eval_20260625.log
#SBATCH -e ref40_eval_20260625.log
#SBATCH --cpus-per-task=8

set -euo pipefail

NF=/appsnew/home/syfan/softwares/nextflow/nextflow
PIPE=/lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf
INPUT=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260702-ref_40_20260625_samples/mqres.csv
BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260702-ref_40_20260625_samples

run_case() {
  local name="$1"
  local grid="$2"
  echo "=== Running ${name} ==="
  ${NF} run "${PIPE}" \
    -profile aipt_ref_40,alioth_slurm,singularity \
    --input "${INPUT}" \
    --outdir "${BASE}/${name}" \
    --grid_search_result "${grid}"
}

# run_case fixed_ez25 /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_all_fixed_ez25
# run_case even_split /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_even_split
run_case fixed_combo_all /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/fixed_combo_all

echo "All evaluation runs finished."
