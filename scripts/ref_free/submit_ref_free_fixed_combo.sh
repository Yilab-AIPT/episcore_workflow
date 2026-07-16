#!/usr/bin/bash
# Submit reference-free sweep with a fixed episcore/zscore combo pair.
#
# Fixed combo defaults:
#     episcore : threshold=0.5,  recall=0.65
#     zscore   : threshold=0.85, recall=0.95
#
# Usage:
#     ./submit_ref_free_fixed_combo.sh [-n|--dry-run] [--no-aggregate]

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export COMBO_MODE=fixed
export EP_THRESHOLD=0.5
export EP_RECALL=0.65
export Z_THRESHOLD=0.85
export Z_RECALL=0.95
export OUTPUT_BASE=${OUTPUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260713-ref_free_fixed_combo}
export MAX_ARRAY_JOBS=50
export REPEATS_PER_JOB=200

exec "${SCRIPT_DIR}/submit_ref_free.sh" "$@"
