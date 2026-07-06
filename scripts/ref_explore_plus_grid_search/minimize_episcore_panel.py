#!/usr/bin/env python3
"""Greedy per-chr episcore combo shrink while keeping ezscore MCC=1.

Fixes the 40-sample reference draw, zscore combos, and ezscore reference list,
then walks chromosomes (lowest recall first) and raises episcore recall (smaller
CpG panel) whenever ezscore MCC stays at 1.0 on the analyze set.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

from grid_search_ref40 import (
    CHR_LIST,
    DEFAULT_CUTOFF,
    MAX_RECALL,
    _build_dense,
    _per_chr_metrics,
    _pick_chr_combo,
    _read_sample_list,
    compute_episcore,
    compute_zscore,
)
from ref40_score_eval import SCORE_CUTOFF, STRONG_CUTOFF

console = Console()


def _assign_pred_label_matrix(
    score_mat: np.ndarray,
    cutoff: float = SCORE_CUTOFF,
    strong: float = STRONG_CUTOFF,
) -> np.ndarray:
    """Return bool array (n_sample,) — True when sample is positive (T or Gray)."""
    above_strong = score_mat > strong
    above_cutoff = (score_mat > cutoff) & ~above_strong
    any_strong = above_strong.any(axis=0)
    any_gray = above_cutoff.any(axis=0)
    return any_strong | any_gray


def _mcc_from_binary(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    tp = int((y_true & y_pred).sum())
    tn = int((~y_true & ~y_pred).sum())
    fp = int((~y_true & y_pred).sum())
    fn = int((y_true & ~y_pred).sum())
    denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0


def build_score_mats(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    ez_idx: np.ndarray,
    n_sample: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ep_combo_index = {c: i for i, c in enumerate(ep_combos)}
    z_combo_index = {c: i for i, c in enumerate(z_combos)}
    ez_score_all = np.empty((len(CHR_LIST), n_sample), dtype=np.float64)
    episcore_best = np.empty((len(CHR_LIST), n_sample), dtype=np.float64)
    zscore_best = np.empty((len(CHR_LIST), n_sample), dtype=np.float64)
    for hi, chrom in enumerate(CHR_LIST):
        ep_vec = episcore_all[ep_combo_index[ep_best[chrom]], hi, :]
        z_vec = zscore_all[z_combo_index[z_best[chrom]], hi, :]
        episcore_best[hi] = ep_vec
        zscore_best[hi] = z_vec
        combined = z_vec + ep_vec
        with np.errstate(invalid="ignore"):
            ez_mean = np.nanmean(combined[ez_idx])
            ez_std = np.nanstd(combined[ez_idx], ddof=0)
        ez_mean = ez_mean if np.isfinite(ez_mean) else 0.0
        ez_std_safe = ez_std if ez_std > 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            ez_score_all[hi] = (combined - ez_mean) / ez_std_safe
    return ez_score_all, episcore_best, zscore_best


def eval_ezscore_mcc(
    ez_score_all: np.ndarray,
    episcore_best: np.ndarray,
    zscore_best: np.ndarray,
    eval_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    cutoff: float = SCORE_CUTOFF,
) -> Tuple[float, int, int, int, int]:
    labels = label_arr[eval_idx]
    ff = ff_arr[eval_idx]
    ez_sub = ez_score_all[:, eval_idx]
    ep_sub = episcore_best[:, eval_idx]
    z_sub = zscore_best[:, eval_idx]

    ep_pos = _assign_pred_label_matrix(ep_sub, cutoff)
    z_pos = _assign_pred_label_matrix(z_sub, cutoff)
    both_normal = ~(ep_pos | z_pos)

    any_ez_pos = (ez_sub > cutoff).any(axis=0)
    y_pred = any_ez_pos & ~both_normal

    y_true = np.array([str(l).startswith("T") for l in labels], dtype=bool)
    mask_t15 = (labels == "T15") & (ff < 0.01)
    y_pred[mask_t15] = False

    mcc = _mcc_from_binary(y_true, y_pred)
    tp = int((y_true & y_pred).sum())
    tn = int((~y_true & ~y_pred).sum())
    fp = int((~y_true & y_pred).sum())
    fn = int((y_true & ~y_pred).sum())
    return mcc, tp, tn, fp, fn


def load_combo_dict(path: Path) -> Tuple[Dict[str, Tuple[float, float]], Dict[str, bool], float]:
    df = pd.read_csv(path)
    combos = {
        str(row.chr): (float(row.threshold), float(row.recall))
        for row in df.itertuples(index=False)
    }
    has_target = {str(row.chr): bool(row.has_target) for row in df.itertuples(index=False)}
    min_recall = float(df["min_recall"].iloc[0])
    return combos, has_target, min_recall


def eval_ezscore_predictions(
    ez_score_all: np.ndarray,
    episcore_best: np.ndarray,
    zscore_best: np.ndarray,
    eval_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    sample_ids: Sequence[str],
    cutoff: float = SCORE_CUTOFF,
) -> pd.DataFrame:
    """Per-sample ezscore binary calls on the analyze set."""
    labels = label_arr[eval_idx]
    ff = ff_arr[eval_idx]
    ez_sub = ez_score_all[:, eval_idx]
    ep_sub = episcore_best[:, eval_idx]
    z_sub = zscore_best[:, eval_idx]

    ep_pos = _assign_pred_label_matrix(ep_sub, cutoff)
    z_pos = _assign_pred_label_matrix(z_sub, cutoff)
    both_normal = ~(ep_pos | z_pos)
    any_ez_pos = (ez_sub > cutoff).any(axis=0)
    y_pred = any_ez_pos & ~both_normal

    y_true = np.array([str(l).startswith("T") for l in labels], dtype=bool)
    mask_t15 = (labels == "T15") & (ff < 0.01)
    y_pred[mask_t15] = False

    status = np.full(len(eval_idx), "UNK", dtype=object)
    status[y_true & y_pred] = "TP"
    status[y_true & ~y_pred] = "FN"
    status[~y_true & y_pred] = "FP"
    status[~y_true & ~y_pred] = "TN"

    return pd.DataFrame(
        {
            "sample": [sample_ids[i] for i in eval_idx],
            "label": labels,
            "ff_before_mq": ff,
            "y_true": y_true,
            "y_pred": y_pred,
            "match_status_ezscore": status,
        }
    )


def build_episcore_per_chr_metrics(
    episcore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    cutoff: float = SCORE_CUTOFF,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, bool]]:
    """Per-chr episcore combo metrics on the analyze set (for _pick_chr_combo)."""
    per_chr_metrics: Dict[str, pd.DataFrame] = {}
    has_target: Dict[str, bool] = {}
    grid_labels = label_arr[analyze_idx]
    grid_is_normal = grid_labels == "Normal"

    for hi, chrom in enumerate(CHR_LIST):
        target_label = chrom.replace("chr", "T")
        is_target = grid_labels == target_label
        has_target[chrom] = bool(is_target.any())
        keep = (is_target | grid_is_normal) if has_target[chrom] else grid_is_normal
        scores_cs = episcore_all[:, hi, :][:, analyze_idx][:, keep]
        per_chr_metrics[chrom] = _per_chr_metrics(
            scores_cs,
            list(ep_combos),
            is_target[keep],
            grid_is_normal[keep],
            cutoff,
        )
    return per_chr_metrics, has_target


def sweep_chr9_recall_errors(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    per_chr_metrics: Dict[str, pd.DataFrame],
    has_target: Dict[str, bool],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    sample_ids: Sequence[str],
    recall_start: float = 0.11,
    recall_end: float = 0.65,
    recall_step: float = 0.01,
    cutoff: float = SCORE_CUTOFF,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Sweep chr9 episcore recall; record FP/FN samples vs baseline."""
    n_sample = episcore_all.shape[2]
    chrom = "chr9"
    targets = np.arange(recall_start, recall_end + recall_step / 2, recall_step)

    baseline = dict(ep_best)
    ez_b, ep_b, z_b = build_score_mats(
        baseline, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
    )
    base_pred = eval_ezscore_predictions(
        ez_b, ep_b, z_b, analyze_idx, label_arr, ff_arr, sample_ids, cutoff,
    )
    base_mcc, tp, tn, fp, fn = eval_ezscore_mcc(
        ez_b, ep_b, z_b, analyze_idx, label_arr, ff_arr, cutoff,
    )

    summary_rows: List[dict] = []
    error_rows: List[dict] = []

    for target_recall in targets:
        trial = dict(ep_best)
        picked = _pick_chr_combo(
            per_chr_metrics[chrom], has_target[chrom], float(target_recall),
        )
        if picked is None:
            summary_rows.append(
                {
                    "chr9_target_recall": round(float(target_recall), 2),
                    "chr9_threshold": np.nan,
                    "chr9_recall": np.nan,
                    "mcc": np.nan,
                    "tp": np.nan,
                    "tn": np.nan,
                    "fp": np.nan,
                    "fn": np.nan,
                    "n_fp_samples": np.nan,
                    "n_fn_samples": np.nan,
                    "status": "no_valid_combo",
                }
            )
            continue

        trial[chrom] = picked
        ez_t, ep_t, z_t = build_score_mats(
            trial, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
        )
        mcc, tp, tn, fp, fn = eval_ezscore_mcc(
            ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, cutoff,
        )
        pred = eval_ezscore_predictions(
            ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, sample_ids, cutoff,
        )
        merged = pred.merge(
            base_pred[["sample", "match_status_ezscore"]].rename(
                columns={"match_status_ezscore": "baseline_status"},
            ),
            on="sample",
            how="left",
        )
        merged["chr9_target_recall"] = round(float(target_recall), 2)
        merged["chr9_threshold"] = picked[0]
        merged["chr9_recall"] = picked[1]

        err = merged[
            merged["match_status_ezscore"].isin(["FP", "FN"])
        ].copy()
        err["error_type"] = err["match_status_ezscore"]
        err["baseline_was_correct"] = err["baseline_status"].isin(["TP", "TN"])
        error_rows.extend(err.to_dict("records"))

        summary_rows.append(
            {
                "chr9_target_recall": round(float(target_recall), 2),
                "chr9_threshold": picked[0],
                "chr9_recall": picked[1],
                "mcc": mcc,
                "tp": tp,
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "n_fp_samples": fp,
                "n_fn_samples": fn,
                "status": "ok",
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    error_df = pd.DataFrame(error_rows)
    if not error_df.empty:
        cols = [
            "chr9_target_recall", "chr9_threshold", "chr9_recall",
            "sample", "label", "ff_before_mq",
            "baseline_status", "match_status_ezscore", "error_type",
            "baseline_was_correct",
        ]
        error_df = error_df[[c for c in cols if c in error_df.columns]]

    summary_df.attrs["baseline_mcc"] = base_mcc
    return summary_df, error_df


def sweep_chr9_sample_chr_scores(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    per_chr_metrics: Dict[str, pd.DataFrame],
    has_target: Dict[str, bool],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    sample_ids: Sequence[str],
    recall_start: float = 0.11,
    recall_end: float = 0.65,
    recall_step: float = 0.01,
) -> pd.DataFrame:
    """Per-sample chr9 episcore/zscore/ezscore while sweeping chr9 recall."""
    n_sample = episcore_all.shape[2]
    chrom = "chr9"
    hi = CHR_LIST.index(chrom)
    targets = np.arange(recall_start, recall_end + recall_step / 2, recall_step)
    rows: List[dict] = []

    for target_recall in targets:
        trial = dict(ep_best)
        picked = _pick_chr_combo(
            per_chr_metrics[chrom], has_target[chrom], float(target_recall),
        )
        if picked is None:
            continue
        trial[chrom] = picked
        ez_t, ep_t, z_t = build_score_mats(
            trial, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
        )
        for i in analyze_idx:
            rows.append(
                {
                    "recall": round(float(target_recall), 2),
                    "chr9_threshold": picked[0],
                    "chr9_recall": picked[1],
                    "sample": sample_ids[i],
                    "label": label_arr[i],
                    "episcore_chr9": ep_t[hi, i],
                    "zscore_chr9": z_t[hi, i],
                    "ezscore_chr9": ez_t[hi, i],
                }
            )
    return pd.DataFrame(rows)


def sweep_min_recall_ezscore_mcc(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    per_chr_metrics: Dict[str, pd.DataFrame],
    has_target: Dict[str, bool],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    recall_start: float = 0.11,
    recall_end: float = 0.65,
    recall_step: float = 0.01,
    cutoff: float = SCORE_CUTOFF,
) -> pd.DataFrame:
    """Raise min_recall incrementally; update chrs at or below each new floor."""
    n_sample = episcore_all.shape[2]
    working = dict(ep_best)
    targets = np.arange(recall_start, recall_end + recall_step / 2, recall_step)
    rows: List[dict] = []

    for target_min_recall in targets:
        n_updated = 0
        blocked: List[str] = []
        for chrom in CHR_LIST:
            if working[chrom][1] > float(target_min_recall) + 1e-9:
                continue
            picked = _pick_chr_combo(
                per_chr_metrics[chrom], has_target[chrom], float(target_min_recall),
            )
            if picked is None:
                blocked.append(chrom)
                continue
            if picked != working[chrom]:
                n_updated += 1
            working[chrom] = picked

        actual_min_recall = min(v[1] for v in working.values())
        if blocked:
            rows.append(
                {
                    "target_min_recall": round(float(target_min_recall), 2),
                    "actual_min_recall": actual_min_recall,
                    "mcc": np.nan,
                    "tp": np.nan,
                    "tn": np.nan,
                    "fp": np.nan,
                    "fn": np.nan,
                    "n_chr_updated": n_updated,
                    "blocked_chrs": ",".join(blocked),
                    "status": "blocked",
                }
            )
            continue

        ez_t, ep_t, z_t = build_score_mats(
            working, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
        )
        mcc, tp, tn, fp, fn = eval_ezscore_mcc(
            ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, cutoff,
        )
        rows.append(
            {
                "target_min_recall": round(float(target_min_recall), 2),
                "actual_min_recall": actual_min_recall,
                "mcc": mcc,
                "tp": tp,
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "n_chr_updated": n_updated,
                "blocked_chrs": "",
                "status": "ok",
            }
        )

    return pd.DataFrame(rows)


def extract_episcore_at_min_recall(
    ep_best: Dict[str, Tuple[float, float]],
    per_chr_metrics: Dict[str, pd.DataFrame],
    has_target: Dict[str, bool],
    target_min_recall: float,
    recall_start: float = 0.11,
    recall_step: float = 0.01,
) -> Tuple[Dict[str, Tuple[float, float]], float]:
    """Apply incremental min_recall floor updates up to ``target_min_recall``.

    Mirrors the update loop in :func:`sweep_min_recall_ezscore_mcc` and returns
    the per-chr episcore combo state at the target floor.
    """
    working = dict(ep_best)
    targets = np.arange(recall_start, target_min_recall + recall_step / 2, recall_step)
    for target in targets:
        for chrom in CHR_LIST:
            if working[chrom][1] > float(target) + 1e-9:
                continue
            picked = _pick_chr_combo(
                per_chr_metrics[chrom], has_target[chrom], float(target),
            )
            if picked is None:
                raise ValueError(
                    f"No episcore combo for {chrom} at min_recall={target:.2f}"
                )
            working[chrom] = picked
    return working, min(v[1] for v in working.values())


def _shared_ep_combo(ep_best: Dict[str, Tuple[float, float]]) -> Tuple[float, float]:
    """Return the shared episcore combo for fixed-combo modes."""
    combo = ep_best[CHR_LIST[0]]
    for chrom in CHR_LIST[1:]:
        if ep_best[chrom] != combo:
            raise ValueError(f"Expected one shared episcore combo; {chrom} differs")
    return combo


def sweep_fixed_episcore_recall_mcc(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    recall_start: float,
    recall_end: float,
    recall_step: float = 0.01,
    cutoff: float = SCORE_CUTOFF,
) -> pd.DataFrame:
    """Sweep a shared episcore (threshold, recall) and evaluate ezscore MCC."""
    n_sample = episcore_all.shape[2]
    threshold, _ = _shared_ep_combo(ep_best)
    combo_set = set(ep_combos)
    targets = np.arange(recall_start, recall_end + recall_step / 2, recall_step)
    rows: List[dict] = []

    for target_recall in targets:
        recall = round(float(target_recall), 2)
        picked = (threshold, recall)
        if picked not in combo_set:
            rows.append(
                {
                    "target_recall": recall,
                    "threshold": threshold,
                    "recall": recall,
                    "mcc": np.nan,
                    "tp": np.nan,
                    "tn": np.nan,
                    "fp": np.nan,
                    "fn": np.nan,
                    "status": "missing_combo",
                }
            )
            continue

        trial = {chrom: picked for chrom in CHR_LIST}
        ez_t, ep_t, z_t = build_score_mats(
            trial, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
        )
        mcc, tp, tn, fp, fn = eval_ezscore_mcc(
            ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, cutoff,
        )
        rows.append(
            {
                "target_recall": recall,
                "threshold": picked[0],
                "recall": picked[1],
                "mcc": mcc,
                "tp": tp,
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "status": "ok",
            }
        )
    return pd.DataFrame(rows)


def shrink_fixed_episcore_recall(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    has_target: Dict[str, bool],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    recall_end: float = 0.75,
    recall_step: float = 0.01,
    cutoff: float = SCORE_CUTOFF,
) -> Tuple[Dict[str, Tuple[float, float]], pd.DataFrame]:
    """Raise shared episcore recall while keeping baseline ezscore MCC/FP/FN."""
    n_sample = episcore_all.shape[2]
    threshold, base_recall = _shared_ep_combo(ep_best)
    combo_set = set(ep_combos)
    working = dict(ep_best)

    ez_b, ep_b, z_b = build_score_mats(
        working, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
    )
    base_mcc, base_tp, base_tn, base_fp, base_fn = eval_ezscore_mcc(
        ez_b, ep_b, z_b, analyze_idx, label_arr, ff_arr, cutoff,
    )

    best_recall = base_recall
    log_rows: List[dict] = []
    for target_recall in np.arange(base_recall + recall_step, recall_end + recall_step / 2, recall_step):
        recall = round(float(target_recall), 2)
        picked = (threshold, recall)
        if picked not in combo_set:
            log_rows.append(
                {
                    "target_recall": recall,
                    "threshold": threshold,
                    "recall": recall,
                    "mcc": np.nan,
                    "status": "missing_combo",
                }
            )
            continue

        trial = {chrom: picked for chrom in CHR_LIST}
        ez_t, ep_t, z_t = build_score_mats(
            trial, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample,
        )
        mcc, tp, tn, fp, fn = eval_ezscore_mcc(
            ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, cutoff,
        )
        kept = (
            np.isclose(mcc, base_mcc)
            and tp == base_tp
            and tn == base_tn
            and fp == base_fp
            and fn == base_fn
        )
        log_rows.append(
            {
                "target_recall": recall,
                "threshold": picked[0],
                "recall": picked[1],
                "mcc": mcc,
                "tp": tp,
                "tn": tn,
                "fp": fp,
                "fn": fn,
                "status": "accepted" if kept else "rejected",
            }
        )
        if kept:
            best_recall = recall
            working = trial
        else:
            break

    optimized = {chrom: (threshold, best_recall) for chrom in CHR_LIST}
    log_df = pd.DataFrame(log_rows)
    log_df.attrs["baseline_mcc"] = base_mcc
    log_df.attrs["baseline_recall"] = base_recall
    log_df.attrs["optimized_recall"] = best_recall
    return optimized, log_df


def combo_df(
    ep_best: Dict[str, Tuple[float, float]],
    has_target: Dict[str, bool],
    min_recall: float,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chr": CHR_LIST,
            "threshold": [ep_best[c][0] for c in CHR_LIST],
            "recall": [ep_best[c][1] for c in CHR_LIST],
            "has_target": [has_target[c] for c in CHR_LIST],
            "min_recall": [min_recall] * len(CHR_LIST),
        }
    )


def greedy_shrink_episcore(
    ep_best: Dict[str, Tuple[float, float]],
    z_best: Dict[str, Tuple[float, float]],
    has_target: Dict[str, bool],
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    ep_combos: Sequence[Tuple[float, float]],
    z_combos: Sequence[Tuple[float, float]],
    ez_idx: np.ndarray,
    analyze_idx: np.ndarray,
    label_arr: np.ndarray,
    ff_arr: np.ndarray,
    cutoff: float = SCORE_CUTOFF,
    mcc_target: float = 1.0,
    max_rounds: int = 5,
) -> Tuple[Dict[str, Tuple[float, float]], pd.DataFrame]:
    n_sample = episcore_all.shape[2]
    working = dict(ep_best)
    log_rows: List[dict] = []

    ez_mat, ep_mat, z_mat = build_score_mats(
        working, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample
    )
    base_mcc, tp, tn, fp, fn = eval_ezscore_mcc(
        ez_mat, ep_mat, z_mat, analyze_idx, label_arr, ff_arr, cutoff
    )
    console.print(
        f"Baseline ezscore MCC={base_mcc:.6f}  TP={tp} TN={tn} FP={fp} FN={fn}  "
        f"min(recall)={min(v[1] for v in working.values()):.2f}"
    )
    if not np.isclose(base_mcc, mcc_target):
        raise click.ClickException(
            f"Baseline MCC {base_mcc:.6f} != target {mcc_target}; fix inputs before shrinking"
        )

    for round_idx in range(1, max_rounds + 1):
        chr_order = sorted(CHR_LIST, key=lambda c: working[c][1])
        round_updates = 0
        console.print(f"\n[bold]Round {round_idx}[/bold]  (chr order by recall)")
        for chrom in chr_order:
            cur_thr, cur_rec = working[chrom]
            candidates = [
                c for c in ep_combos
                if c[1] > cur_rec + 1e-9 and c[1] <= MAX_RECALL + 1e-9
            ]
            candidates.sort(key=lambda c: (-c[1], c[0]))
            if not candidates:
                continue

            best_candidate = None
            for thr, rec in candidates:
                trial = dict(working)
                trial[chrom] = (thr, rec)
                ez_t, ep_t, z_t = build_score_mats(
                    trial, z_best, episcore_all, zscore_all, ep_combos, z_combos, ez_idx, n_sample
                )
                mcc, tp, tn, fp, fn = eval_ezscore_mcc(
                    ez_t, ep_t, z_t, analyze_idx, label_arr, ff_arr, cutoff
                )
                if np.isclose(mcc, mcc_target) and fp == 0 and fn == 0:
                    best_candidate = (thr, rec)
                    break

            if best_candidate is None:
                continue

            working[chrom] = best_candidate
            round_updates += 1
            log_rows.append(
                {
                    "round": round_idx,
                    "chr": chrom,
                    "old_threshold": cur_thr,
                    "old_recall": cur_rec,
                    "new_threshold": best_candidate[0],
                    "new_recall": best_candidate[1],
                    "status": "updated",
                }
            )
            console.print(
                f"  {chrom}: recall {cur_rec:.2f} -> {best_candidate[1]:.2f}  "
                f"(thr {cur_thr:g} -> {best_candidate[0]:g})"
            )

        if round_updates == 0:
            console.print(f"  Round {round_idx}: no updates — stopping")
            break

    return working, pd.DataFrame(log_rows)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--result-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--min-ff", default=0.01, show_default=True, type=float)
@click.option("--cutoff", default=DEFAULT_CUTOFF, show_default=True, type=float)
@click.option("--output", default=None, type=click.Path(dir_okay=False),
              help="Write optimized best_combo_episcore.csv (default: stdout only)")
def main(
    input_dir: str,
    result_dir: str,
    min_ff: float,
    cutoff: float,
    output: str | None,
) -> None:
    """Shrink episcore CpG panels per chr while keeping ezscore MCC=1."""
    input_path = Path(input_dir)
    result_path = Path(result_dir)

    ep_best, has_target, _ = load_combo_dict(result_path / "best_combo_episcore.csv")
    z_best, _, _ = load_combo_dict(result_path / "best_combo_zscore.csv")
    ref_samples = _read_sample_list(result_path / "best_reference_samples.txt")
    ez_samples = _read_sample_list(result_path / "best_ezscore_ref_20_samples.txt")

    meta = pd.read_csv(input_path / "meta.csv").drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")

    console.print("[cyan]Loading parquets ...[/cyan]")
    ep_df = pd.read_parquet(input_path / "episcore_grid_search.parquet")
    z_df = pd.read_parquet(input_path / "zscore_grid_search.parquet")

    ep_samples = set(ep_df["sample"].astype(str).unique())
    z_samples = set(z_df["sample"].astype(str).unique())
    meta_samples = set(meta["sample"])
    ff_pass = set(meta.set_index("sample").index[meta["ff_before_mq"] > min_ff].astype(str))
    ref_set = set(ref_samples)
    ez_set = set(ez_samples)

    universe = sorted(
        (meta_samples & ep_samples & z_samples & ff_pass)
        | (ref_set & meta_samples & ep_samples & z_samples)
        | (ez_set & meta_samples & ep_samples & z_samples)
    )
    samples = universe
    sample_index = {s: i for i, s in enumerate(samples)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}

    ep_df = ep_df[ep_df["sample"].astype(str).isin(samples)]
    z_df = z_df[z_df["sample"].astype(str).isin(samples)]

    meta_idx = meta.set_index("sample").reindex(samples)
    label_arr = meta_idx["label"].astype(str).to_numpy()
    ff_arr = pd.to_numeric(meta_idx["ff_before_mq"], errors="coerce").to_numpy()
    analyze_ff_mask = ff_arr > min_ff

    ref_idx = np.array([sample_index[s] for s in ref_samples if s in sample_index], dtype=np.int64)
    ez_idx = np.array([sample_index[s] for s in ez_samples if s in sample_index], dtype=np.int64)
    if ref_idx.size != len(ref_samples):
        raise click.ClickException("Some reference samples missing from parquet universe")
    if ez_idx.size != len(ez_samples):
        raise click.ClickException("Some ezscore reference samples missing from parquet universe")

    n_sample = len(samples)
    ref_mask = np.zeros(n_sample, dtype=bool)
    ref_mask[ref_idx] = True
    analyze_idx = np.flatnonzero(~ref_mask & analyze_ff_mask)

    ep_combos, ep_arrays = _build_dense(
        ep_df,
        ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
        sample_index,
        chr_index,
    )
    z_combos, z_arrays = _build_dense(z_df, ["percentage"], sample_index, chr_index)

    console.print("[cyan]Computing score tensors ...[/cyan]")
    episcore_all = compute_episcore(ep_arrays[0], ep_arrays[1], ep_arrays[2], ep_arrays[3], ref_idx)
    zscore_all = compute_zscore(z_arrays[0], ref_idx)

    optimized, log_df = greedy_shrink_episcore(
        ep_best,
        z_best,
        has_target,
        episcore_all,
        zscore_all,
        ep_combos,
        z_combos,
        ez_idx,
        analyze_idx,
        label_arr,
        ff_arr,
        cutoff=cutoff,
    )

    new_min_recall = min(v[1] for v in optimized.values())
    out_df = combo_df(optimized, has_target, new_min_recall)

    console.print("\n[bold]Optimized combo summary[/bold]")
    console.print(out_df.to_string(index=False))
    console.print(f"\nmin(recall): {min(v[1] for v in ep_best.values()):.2f} -> {new_min_recall:.2f}")

    if output:
        out_path = Path(output)
        out_df.to_csv(out_path, index=False)
        log_path = out_path.with_name(out_path.stem + "_shrink_log.tsv")
        log_df.to_csv(log_path, sep="\t", index=False)
        console.print(f"[green]OK[/green] Wrote {out_path}")
        console.print(f"[green]OK[/green] Wrote {log_path}")


if __name__ == "__main__":
    try:
        main()
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
