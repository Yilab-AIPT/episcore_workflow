"""BAM path helpers shared by build_120k_panel scripts."""

from __future__ import annotations

import os
from pathlib import Path


def bam_index_path(bam_path: Path) -> Path | None:
    """Return path to an existing .bai for this BAM, or None."""
    bam_path = Path(bam_path)
    if not bam_path.is_file():
        return None
    real = Path(os.path.realpath(bam_path))
    for candidate in (Path(str(real) + ".bai"), Path(str(bam_path) + ".bai")):
        if candidate.is_file():
            return candidate
    return None


def bam_has_index(bam_path: Path) -> bool:
    """True when the BAM exists and a sidecar .bai is available for pileup."""
    return bam_index_path(bam_path) is not None
