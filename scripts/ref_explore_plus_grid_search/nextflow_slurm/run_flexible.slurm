#!/usr/bin/bash
#SBATCH -J flexible
#SBATCH --partition=cn-long
#SBATCH -o flexible.log
#SBATCH -e flexible.log
#SBATCH --cpus-per-task=16

# /appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
# -profile aipt_ref_40,alioth_slurm,singularity \
# --input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/mqres.csv \
# --outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/flexible_combo_all \
# --grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_all

# /appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
# -profile aipt_ref_40,alioth_slurm,singularity \
# --input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/mqres.csv \
# --outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/flexible_combo_even_split \
# --grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_even_split


/appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
-profile aipt_ref_40,alioth_slurm,singularity \
--input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/add_mqres.csv \
--outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/flexible_combo_all \
--grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_all

/appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
-profile aipt_ref_40,alioth_slurm,singularity \
--input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/add_mqres.csv \
--outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/flexible_combo_even_split \
--grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/flexible_combo_even_split