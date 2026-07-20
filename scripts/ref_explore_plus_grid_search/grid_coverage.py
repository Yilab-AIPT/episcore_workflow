"""Shared episcore / zscore grid-search parquet coverage helpers."""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

N_AUTOSOMES = 22
Combo = Tuple[float, float]


def _as_float_pair(threshold, recall) -> Combo:
    return (float(threshold), float(recall))


def combo_key(threshold: float, recall: float) -> str:
    return f"threshold={threshold:g}, recall={recall:g}"


def majority_combos(
    df: pd.DataFrame,
    universe: Sequence[str],
    *,
    n_chr: int = N_AUTOSOMES,
    majority_frac: float = 0.95,
) -> List[Combo]:
    """Combos present with full chr coverage for >= majority_frac of universe."""
    if not universe:
        return []
    sub = df[df["sample"].astype(str).isin(set(universe))]
    g = sub.groupby(["sample", "threshold", "recall"], sort=False)["chr"].nunique()
    full = g[g >= n_chr].reset_index()
    if full.empty:
        return []
    counts = full.groupby(["threshold", "recall"])["sample"].nunique()
    keep = counts[counts >= majority_frac * len(universe)]
    return sorted((_as_float_pair(t, r) for t, r in keep.index), key=lambda x: (x[0], x[1]))


def full_coverage_index(
    df: pd.DataFrame,
    *,
    n_chr: int = N_AUTOSOMES,
) -> set:
    """Set of (sample, threshold, recall) with >= n_chr chromosomes."""
    g = df.groupby(["sample", "threshold", "recall"], sort=False)["chr"].nunique()
    full = g[g >= n_chr].reset_index()
    return {
        (str(s), float(t), float(r))
        for s, t, r in zip(full["sample"], full["threshold"], full["recall"])
    }


def find_missing_coverage(
    df: pd.DataFrame,
    universe: Sequence[str],
    combos: Optional[Sequence[Combo]] = None,
    *,
    n_chr: int = N_AUTOSOMES,
    majority_frac: float = 0.95,
) -> pd.DataFrame:
    """Return missing (sample, threshold, recall) rows vs expected combos."""
    if combos is None:
        combos = majority_combos(df, universe, n_chr=n_chr, majority_frac=majority_frac)
    present = full_coverage_index(df, n_chr=n_chr)
    rows = []
    for sample in universe:
        for thr, rec in combos:
            if (str(sample), float(thr), float(rec)) not in present:
                rows.append(
                    {
                        "sample": str(sample),
                        "threshold": float(thr),
                        "recall": float(rec),
                    }
                )
    return pd.DataFrame(rows, columns=["sample", "threshold", "recall"])


def format_missing_error(
    score_name: str,
    missing: pd.DataFrame,
    *,
    max_show: int = 20,
) -> str:
    if missing.empty:
        return ""
    lines = [
        f"{score_name} grid coverage incomplete: {len(missing)} missing "
        f"sample×combo (need {N_AUTOSOMES} chr each)."
    ]
    show = missing.head(max_show)
    for _, row in show.iterrows():
        lines.append(
            f"  - {row['sample']}: {combo_key(row['threshold'], row['recall'])}"
        )
    if len(missing) > max_show:
        lines.append(f"  ... and {len(missing) - max_show} more")
    return "\n".join(lines)


def assert_table_coverage(
    df: pd.DataFrame,
    universe: Sequence[str],
    score_name: str,
    combos: Optional[Sequence[Combo]] = None,
    *,
    n_chr: int = N_AUTOSOMES,
    majority_frac: float = 0.95,
) -> None:
    """Raise ValueError listing missing sample×combo coverage."""
    missing = find_missing_coverage(
        df, universe, combos, n_chr=n_chr, majority_frac=majority_frac
    )
    if not missing.empty:
        raise ValueError(format_missing_error(score_name, missing))


def missing_from_dense(
    values: np.ndarray,
    sample_names: Sequence[str],
    combos: Sequence[Combo],
) -> pd.DataFrame:
    """Detect sample×combo with all-NaN chromosomes in a dense [combo,chr,sample] array."""
    if values.ndim != 3:
        raise ValueError(f"Expected [n_combo, n_chr, n_sample], got shape {values.shape}")
    rows = []
    for ci, (thr, rec) in enumerate(combos):
        # all-chr NaN for a sample
        all_nan = np.all(~np.isfinite(values[ci]), axis=0)
        for si in np.flatnonzero(all_nan):
            rows.append(
                {
                    "sample": str(sample_names[si]),
                    "threshold": float(thr),
                    "recall": float(rec),
                }
            )
    return pd.DataFrame(rows, columns=["sample", "threshold", "recall"])


def assert_dense_coverage(
    values: np.ndarray,
    sample_names: Sequence[str],
    combos: Sequence[Combo],
    score_name: str,
) -> None:
    missing = missing_from_dense(values, sample_names, combos)
    if not missing.empty:
        raise ValueError(format_missing_error(score_name, missing))


def fmt_float(value: float) -> str:
    """Match zscore grid dir / file float formatting (``%g``)."""
    return f"{float(value):g}"
