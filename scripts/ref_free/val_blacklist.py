"""Shared val-set sample blacklist for ref_free downstream analysis."""

from __future__ import annotations

VAL_BLACKLIST = frozenset({"PTAY1559P9S1", "PTAY1236P9S1"})


def drop_blacklisted(df, sample_col: str = "sample"):
    """Return a copy of ``df`` without blacklisted sample IDs."""
    import pandas as pd

    if sample_col not in df.columns:
        return df
    mask = ~df[sample_col].astype(str).isin(VAL_BLACKLIST)
    return df.loc[mask].copy()
