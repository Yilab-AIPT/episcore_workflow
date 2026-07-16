#!/usr/bin/bash
# Submit 48+48 split ref_free_ezscore sweep with fixed combo pair.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export COMBO_MODE=fixed
export EP_THRESHOLD=0.5
export EP_RECALL=0.65
export Z_THRESHOLD=0.85
export Z_RECALL=0.95
export OUTPUT_BASE=${OUTPUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260715-ref_free_ezscore_fixed_combo}

exec "${SCRIPT_DIR}/submit_ref_free_ezscore_all.sh" "$@"
