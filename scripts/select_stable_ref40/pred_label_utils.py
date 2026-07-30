#!/usr/bin/env python3
"""Ezscore / score pred_label helpers (cutoff 4.5 strong, 3.0 gray)."""

from __future__ import annotations

from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd

CHR_NUMS = list(range(1, 23))
STRONG_CUTOFF = 4.5
GRAY_CUTOFF = 3.0


def assign_pred_label_from_scores(
    scores: Sequence[float],
    *,
    strong_cutoff: float = STRONG_CUTOFF,
    gray_cutoff: float = GRAY_CUTOFF,
) -> str:
    """Assign pred_label from per-chr scores (chr1..chr22 order).

    Matches ``update_samplesheet.py``:
      z > 4.5       -> T{n}
      3 <= z <= 4.5 -> Gray_T{n}
      otherwise     -> Normal (if no hits)
    """
    t_labels: List[str] = []
    gray_labels: List[str] = []
    for n, z in zip(CHR_NUMS, scores):
        if z is None or (isinstance(z, float) and np.isnan(z)):
            continue
        if z > strong_cutoff:
            t_labels.append(f"T{n}")
        elif z >= gray_cutoff:
            gray_labels.append(f"Gray_T{n}")
    parts = t_labels + gray_labels
    return ",".join(parts) if parts else "Normal"


def assign_pred_labels_matrix(
    score_matrix: np.ndarray,
    *,
    strong_cutoff: float = STRONG_CUTOFF,
    gray_cutoff: float = GRAY_CUTOFF,
) -> List[str]:
    """Assign pred_labels for an (n_samples, 22) score matrix."""
    mat = np.asarray(score_matrix, dtype=np.float64)
    n = mat.shape[0]
    out: List[str] = ["Normal"] * n
    # Precompute string fragments per chromosome
    t_names = [f"T{n_}" for n_ in CHR_NUMS]
    g_names = [f"Gray_T{n_}" for n_ in CHR_NUMS]
    strong = mat > strong_cutoff
    gray = (mat >= gray_cutoff) & (mat <= strong_cutoff)
    for i in range(n):
        parts = [t_names[j] for j in range(22) if strong[i, j]]
        parts.extend(g_names[j] for j in range(22) if gray[i, j])
        if parts:
            out[i] = ",".join(parts)
    return out


def parse_comma_scores(value: object, n_chr: int = 22) -> np.ndarray:
    """Parse meta ``*_zscores`` comma string into float array."""
    arr = np.full(n_chr, np.nan, dtype=np.float64)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return arr
    text = str(value).strip()
    if not text:
        return arr
    parts = text.split(",")
    for i, p in enumerate(parts[:n_chr]):
        p = p.strip()
        if p == "":
            continue
        try:
            arr[i] = float(p)
        except ValueError:
            arr[i] = np.nan
    return arr


def format_comma_scores(values: Iterable[float], precision: int = 6) -> str:
    """Format per-chr scores as comma-separated string for meta samplesheet."""
    fmt = f"{{:.{precision}f}}"
    out = []
    for v in values:
        if v is None or (isinstance(v, float) and (np.isnan(v) or not np.isfinite(v))):
            out.append("")
        else:
            out.append(fmt.format(float(v)))
    return ",".join(out)


def scores_dataframe_to_pred_labels(
    score_df: pd.DataFrame,
    prefix: str = "ezscore",
    *,
    strong_cutoff: float = STRONG_CUTOFF,
    gray_cutoff: float = GRAY_CUTOFF,
) -> pd.Series:
    """Build pred_label Series from wide score columns ``{prefix}_chr{n}``."""
    cols = [f"{prefix}_chr{n}" for n in CHR_NUMS]
    missing = [c for c in cols if c not in score_df.columns]
    if missing:
        raise KeyError(f"Missing score columns: {missing[:5]}")
    mat = score_df[cols].to_numpy(dtype=np.float64)
    return pd.Series(
        assign_pred_labels_matrix(
            mat, strong_cutoff=strong_cutoff, gray_cutoff=gray_cutoff
        ),
        index=score_df.index,
        name=f"pred_label_{prefix}",
    )
