#!/usr/bin/bash
# Submit 48+48 ref_free_ezscore sweep (fixed combo) with ezscore cutoff 4.5.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export CUTOFF=3.0
export EZ_CUTOFF=4.5
export OUTPUT_BASE=${OUTPUT_BASE:-/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260716-ref_free_ezscore_fixed_combo_ez45}

exec "${SCRIPT_DIR}/submit_ref_free_ezscore_fixed_combo.sh" "$@"
