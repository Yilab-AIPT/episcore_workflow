#!/usr/bin/bash
# Full ref-40 pipeline: split-mode=all, combo-mode=fix -> aggregate -> finalize.
#
# Usage:
#     ./submit_ref40_best_of_all.sh [-n|--dry-run] [--no-finalize]
#
# Defaults:
#     output_base   : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260625-ref_40_fixed_combo_all
#     total_repeats : 10000
#     min-ff        : 0.01

set -euo pipefail

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
OUTPUT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260625-ref_40_fixed_combo_all
TOTAL_REPEATS=10000
SPLIT_MODE=all
COMBO_MODE=fix
MIN_FF=0.01
N_EZSCORE_REF=20
EZSCORE_REPEATS=5000
DRY_RUN=0
FINALIZE=1

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --no-finalize) FINALIZE=0; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        --total-repeats) TOTAL_REPEATS=$2; shift 2 ;;
        --min-ff) MIN_FF=$2; shift 2 ;;
        --n-ezscore-ref) N_EZSCORE_REF=$2; shift 2 ;;
        --ezscore-repeats) EZSCORE_REPEATS=$2; shift 2 ;;
        -h|--help)
            sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
args=(
    --input-dir "$INPUT_DIR"
    --output-base "$OUTPUT_BASE"
    --total-repeats "$TOTAL_REPEATS"
    --split-mode "$SPLIT_MODE"
    --combo-mode "$COMBO_MODE"
    --min-ff "$MIN_FF"
    --n-ezscore-ref "$N_EZSCORE_REF"
    --ezscore-repeats "$EZSCORE_REPEATS"
)
[ "$DRY_RUN" = 1 ] && args=(--dry-run "${args[@]}")
[ "$FINALIZE" = 1 ] && args+=(--finalize)

exec "$SCRIPT_DIR/submit_grid_search_ref40.sh" "${args[@]}"
