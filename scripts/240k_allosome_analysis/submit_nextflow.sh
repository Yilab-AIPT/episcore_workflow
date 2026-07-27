#!/usr/bin/bash
# Submit Nextflow for the 10 new 240k allosome samples.
#
#   ./submit_nextflow.sh [--dry-run]
#
# Prerequisites: prepare_inputs.py (writes samplesheet_nf.csv)

set -euo pipefail

DRY_RUN=${DRY_RUN:-0}
for arg in "$@"; do
    case "$arg" in
        -n|--dry-run) DRY_RUN=1 ;;
        -h|--help)
            sed -n '2,8p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "Unknown argument: $arg" >&2; exit 2 ;;
    esac
done

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs

SAMPLESHEET=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples/samplesheet_nf.csv
if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: $SAMPLESHEET not found. Run: singularity exec ... python3 prepare_inputs.py" >&2
    exit 1
fi

if [ "$DRY_RUN" = 1 ]; then
    echo "[DRY-RUN] sbatch --parsable --job-name=allo_nf run_nextflow.slurm"
    exit 0
fi

jobid=$(sbatch --parsable --job-name=allo_nf run_nextflow.slurm)
echo "Submitted allo_nf  job_id=${jobid}  log=logs/allo_nf.log"
