#!/usr/bin/bash
# Orchestrate 240k allosome analysis.
#
# Usage:
#   ./submit_all.sh prepare          # write samplesheets / meta
#   ./submit_all.sh nextflow         # extract beta + FF for new 10 (thres 0.5)
#   ./submit_all.sh prepare          # re-run after nextflow to fill new betas
#   ./submit_all.sh episcore [--test]
#   ./submit_all.sh zscore   [--test]
#   ./submit_all.sh collect          # chrY-FF + recall tables + HTML plots
#
# Or: ./submit_all.sh --help

set -euo pipefail

WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$WORKDIR"
mkdir -p logs

SIF=/lustre1/cqyi/AIPT_2.0/workflow/episcore/containers/common_tools.sif
PY=(singularity exec -B /lustre1,/lustre2,/appsnew "$SIF" python3)

cmd=${1:-}
shift || true

case "$cmd" in
    prepare)
        "${PY[@]}" prepare_inputs.py
        ;;
    nextflow)
        ./submit_nextflow.sh "$@"
        ;;
    episcore)
        ./submit_episcore.sh "$@"
        ;;
    zscore)
        ./submit_zscore.sh "$@"
        ;;
    collect)
        "${PY[@]}" collect_chry_ff.py
        "${PY[@]}" collect_recall_curves.py
        "${PY[@]}" plot_allosome_curves.py
        ;;
    -h|--help|help|"")
        sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
    *)
        echo "Unknown command: $cmd" >&2
        exit 2
        ;;
esac
