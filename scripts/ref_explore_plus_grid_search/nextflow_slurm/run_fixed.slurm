#!/usr/bin/bash
#SBATCH -J fixed
#SBATCH --partition=cn-long
#SBATCH -o fixed.log
#SBATCH -e fixed.log
#SBATCH --cpus-per-task=16

# /appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
# -profile aipt_ref_40,alioth_slurm,singularity \
# --input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/mqres.csv \
# --outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/fixed_combo_all \
# --grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/fixed_combo_all

# /appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
# -profile aipt_ref_40,alioth_slurm,singularity \
# --input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/mqres.csv \
# --outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/fixed_combo_even_split \
# --grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/fixed_combo_even_split


/appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
-profile aipt_ref_40,alioth_slurm,singularity \
--input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/add_mqres.csv \
--outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/fixed_combo_all \
--grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/fixed_combo_all

/appsnew/home/syfan/softwares/nextflow/nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
-profile aipt_ref_40,alioth_slurm,singularity \
--input /lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260626-ref_40_special_sample_check/add_mqres.csv \
--outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_special_sample_check/fixed_combo_even_split \
--grid_search_result /lustre1/cqyi/AIPT_2.0/workflow/episcore/assets/ref_40/fixed_combo_even_split