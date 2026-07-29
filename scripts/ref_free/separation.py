#!/usr/bin/env python3
"""Separation indices between Normal and Trisomy on signal-ratio scores.

Primary (``sep`` / ``sep_auc``): ROC-AUC of ``signal_ratio`` among
``ff_before_mq >= ff_min``.

Secondary (``sep_youden``): max Youden's J = max_t (TPR − FPR) over thresholds
on the score — a single-cutoff operating-point separation.

Also reports Cohen's d (``sep_d``) and mean(T)−mean(N) (``mean_gap``).
"""

from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd


_TRISOMY_RE = re.compile(r"^T\d")


def is_trisomy_label(label: object) -> bool:
    return bool(_TRISOMY_RE.match(str(label)))


def roc_auc(scores: np.ndarray, y_pos: np.ndarray) -> float:
    """Mann–Whitney ROC-AUC; y_pos True = positive class (trisomy)."""
    scores = np.asarray(scores, dtype=float)
    y_pos = np.asarray(y_pos, dtype=bool)
    mask = np.isfinite(scores)
    scores, y_pos = scores[mask], y_pos[mask]
    n_pos = int(y_pos.sum())
    n_neg = int((~y_pos).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, dtype=float)
    ranks[order] = np.arange(1, scores.size + 1, dtype=float)
    sorted_scores = scores[order]
    i = 0
    while i < scores.size:
        j = i + 1
        while j < scores.size and sorted_scores[j] == sorted_scores[i]:
            j += 1
        if j - i > 1:
            avg = 0.5 * (i + 1 + j)
            ranks[order[i:j]] = avg
        i = j
    sum_pos = float(ranks[y_pos].sum())
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def youden_j(scores: np.ndarray, y_pos: np.ndarray) -> float:
    """Max Youden's J over unique score thresholds (positive = high score)."""
    scores = np.asarray(scores, dtype=float)
    y_pos = np.asarray(y_pos, dtype=bool)
    mask = np.isfinite(scores)
    scores, y_pos = scores[mask], y_pos[mask]
    n_pos = int(y_pos.sum())
    n_neg = int((~y_pos).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y_sorted = y_pos[order]
    tps = np.cumsum(y_sorted)
    fps = np.cumsum(~y_sorted)
    tpr = tps / n_pos
    fpr = fps / n_neg
    return float(np.max(tpr - fpr))


def cohens_d(scores: np.ndarray, y_pos: np.ndarray) -> float:
    """Cohen's d between trisomy and Normal score distributions."""
    scores = np.asarray(scores, dtype=float)
    y_pos = np.asarray(y_pos, dtype=bool)
    mask = np.isfinite(scores)
    scores, y_pos = scores[mask], y_pos[mask]
    pos = scores[y_pos]
    neg = scores[~y_pos]
    if pos.size < 2 or neg.size < 2:
        return float("nan")
    m1, m0 = float(np.mean(pos)), float(np.mean(neg))
    v1, v0 = float(np.var(pos, ddof=1)), float(np.var(neg, ddof=1))
    n1, n0 = pos.size, neg.size
    pooled = ((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2)
    if pooled <= 0:
        return float("nan")
    return (m1 - m0) / float(np.sqrt(pooled))


def separation_index(
    df: pd.DataFrame,
    score_col: str,
    *,
    ff_min: float = 0.01,
    ff_col: str = "ff_before_mq",
    label_col: str = "label",
) -> dict:
    """Compute AUC, Youden J, Cohen's d for ff-filtered Normal vs Trisomy."""
    work = df.copy()
    work[ff_col] = pd.to_numeric(work[ff_col], errors="coerce")
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work[work[ff_col] >= ff_min]
    work = work[
        work[label_col].astype(str).eq("Normal")
        | work[label_col].map(is_trisomy_label)
    ]
    empty = {
        "sep": float("nan"),
        "sep_auc": float("nan"),
        "sep_youden": float("nan"),
        "sep_d": float("nan"),
        "mean_gap": float("nan"),
        "n_normal": 0,
        "n_trisomy": 0,
        "score_col": score_col,
        "ff_min": ff_min,
    }
    if work.empty:
        return empty
    y = work[label_col].map(is_trisomy_label).to_numpy()
    scores = work[score_col].to_numpy(dtype=float)
    normals = scores[~y]
    tris = scores[y]
    mean_gap = (
        float(np.nanmean(tris) - np.nanmean(normals))
        if tris.size and normals.size
        else float("nan")
    )
    auc = roc_auc(scores, y)
    return {
        "sep": auc,  # alias for backward compatibility
        "sep_auc": auc,
        "sep_youden": youden_j(scores, y),
        "sep_d": cohens_d(scores, y),
        "mean_gap": mean_gap,
        "n_normal": int((~y).sum()),
        "n_trisomy": int(y.sum()),
        "score_col": score_col,
        "ff_min": ff_min,
    }


def separation_for_cutoffs(
    df: pd.DataFrame,
    cutoffs: Iterable[float],
    *,
    col_fmt: str = "ezscore_signal_ratio_{:g}",
    ff_min: float = 0.01,
) -> dict[float, dict]:
    out = {}
    for c in cutoffs:
        col = col_fmt.format(c)
        if col not in df.columns:
            continue
        out[float(c)] = separation_index(df, col, ff_min=ff_min)
    return out


def format_sep_pair(sep: dict, *, prefix: str = "") -> str:
    """Compact title fragment: AUC / Youden."""
    auc = sep.get("sep_auc", sep.get("sep", float("nan")))
    youden = sep.get("sep_youden", float("nan"))
    p = f"{prefix}" if prefix else ""
    return f"{p}AUC={auc:.3f} J={youden:.3f}"
