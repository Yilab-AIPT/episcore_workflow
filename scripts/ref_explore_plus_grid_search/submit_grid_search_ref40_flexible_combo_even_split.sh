#!/usr/bin/bash
# Submit random ref-40 grid search: split-mode=even_split, combo-mode=flexible.
#
# Usage:
#     ./submit_grid_search_ref40_flexible_combo_even_split.sh [-n|--dry-run]
#
# Defaults:
#     output_base   : /lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_flexible_combo_even_split
#     total_repeats : 10000
#     min-ff        : 0.01

set -euo pipefail

INPUT_DIR=/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng
OUTPUT_BASE=/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260626-ref_40_flexible_combo_even_split
TOTAL_REPEATS=10000
SPLIT_MODE=even_split
COMBO_MODE=flexible
MIN_FF=0.01
DRY_RUN=0

while [ $# -gt 0 ]; do
    case "$1" in
        -n|--dry-run) DRY_RUN=1; shift ;;
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --output-base) OUTPUT_BASE=$2; shift 2 ;;
        --total-repeats) TOTAL_REPEATS=$2; shift 2 ;;
        --min-ff) MIN_FF=$2; shift 2 ;;
        -h|--help)
            sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "$SCRIPT_DIR/submit_grid_search_ref40.sh" \
    $([ "$DRY_RUN" = 1 ] && echo -n "--dry-run ") \
    --input-dir "$INPUT_DIR" \
    --output-base "$OUTPUT_BASE" \
    --total-repeats "$TOTAL_REPEATS" \
    --split-mode "$SPLIT_MODE" \
    --combo-mode "$COMBO_MODE" \
    --min-ff "$MIN_FF" \
    --finalize
