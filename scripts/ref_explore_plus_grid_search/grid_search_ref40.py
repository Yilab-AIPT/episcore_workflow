#!/usr/bin/env python3
"""
Per-repeat random ref-40 grid search for episcore / zscore / ezscore.

For each repeat:
    1. Draw ``ref-n`` (default 40) reference samples from a candidate pool
       (``ref_pool_samples.txt``). Draws are deterministic and unique across the
       full sweep given a fixed ``--seed`` and ``--total-repeats``.
    2. Using that reference, compute reference-normalized scores for *every*
       (threshold, recall) combo and chromosome:
         - episcore (``s_inter``) from ``episcore_grid_search.parquet``
           (hypo/hyper ``z_intra`` + CpG counts), normalized by the 40-sample
           reference mean/std.
         - zscore (chromosome percentage) from ``zscore_grid_search.parquet``,
           normalized by the 40-sample reference mean/std.
    3. Run grid search independently for episcore and zscore on the
       grid-search analyze set. Combo selection strategy depends on ``--mode``:
         - ``dev_test_split`` / ``all`` ("method 3"): per-chromosome best
           (threshold, recall) combo (see
           ``notebooks/aipt_2.0/grid_search_for_smaller_panel_autosomes.ipynb``)
         - ``fix_combo_split`` / ``fix_combo_all``: one shared (threshold,
           recall) combo for all chromosomes
       Sample sets:
         - ``dev_test_split`` / ``fix_combo_split``: grid search on dev only,
           report dev + test
         - ``all`` / ``fix_combo_all``: all analyze samples (after ``--min-ff``)
    4. Build ezscore = z-normalize(zscore + episcore) where, per chromosome,
       zscore uses the best zscore combo and episcore the best episcore combo,
       and the mean/std are taken over the fixed ezscore reference sample list
       (``ezscore_ref_samples.txt``).
    5. Compute confusion (TP/TN/FP/FN) and MCC for episcore, zscore and ezscore.
       ``dev_test_split`` / ``fix_combo_split`` report dev and test sets;
       ``all`` / ``fix_combo_all`` report one ``all`` row. A sample is positive
       when any chromosome score exceeds ``--cutoff``
       (default 3.0); the negative class is ``Normal``.

Samples with ``ff_before_mq <= --min-ff`` are excluded from analyze/grid-search
and MCC evaluation. Reference-pool and ezscore-reference samples are always kept
in the data universe even when below ``--min-ff``, so reference normalization
remains well-defined.

Outputs (under ``<output-base>/randomly_select_ref_40/repeat_{i}/``):
    reference_samples.txt        the 40 drawn reference samples
    best_combo_episcore.csv      chr, threshold, recall (+ has_target)
    best_combo_zscore.csv        chr, threshold, recall (+ has_target)
    metrics.tsv                  score x set -> mcc/tp/tn/fp/fn (+min_recall)
    scores.tsv                   per analyze sample: episcore/zscore/ezscore per chr

A per-slice manifest (``metrics_manifest_{start}_{end}.tsv``) summarising every
repeat in the slice is written under ``<output-base>/randomly_select_ref_40/``.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

# Empty/all-NaN reference slices for sparse (combo, chr) cells are expected and
# handled explicitly; silence the resulting numpy warnings.
warnings.filterwarnings("ignore", category=RuntimeWarning)

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
DEFAULT_CUTOFF = 3.0
MAX_RECALL = 0.99


# ---------------------------------------------------------------------------
# Reference draws (mirrors scripts/reference_explore/calc_random_ref40_scores.py)
# ---------------------------------------------------------------------------

def generate_unique_ref_draws(
    pool_size: int,
    ref_n: int,
    n_repeats: int,
    rng: np.random.Generator,
) -> List[np.ndarray]:
    """Return ``n_repeats`` unique sorted index arrays into the candidate pool."""
    if ref_n <= 0 or ref_n > pool_size:
        raise ValueError(f"Invalid ref_n={ref_n} for pool_size={pool_size}")
    total = math.comb(int(pool_size), int(ref_n))
    if total <= n_repeats:
        raise ValueError(
            f"Pool too small for {n_repeats} unique draws: C({pool_size},{ref_n})={total}"
        )
    seen: set[bytes] = set()
    out: List[np.ndarray] = []
    max_attempts = n_repeats * 30 + 1000
    attempts = 0
    while len(out) < n_repeats and attempts < max_attempts:
        attempts += 1
        choice = rng.choice(pool_size, size=ref_n, replace=False)
        choice.sort()
        key = choice.tobytes()
        if key in seen:
            continue
        seen.add(key)
        out.append(choice.astype(np.int64, copy=False))
    if len(out) < n_repeats:
        raise RuntimeError(
            f"Could only generate {len(out)}/{n_repeats} unique reference draws "
            f"after {attempts} attempts"
        )
    return out


# ---------------------------------------------------------------------------
# Data loading -> dense [combo, chr, sample] arrays
# ---------------------------------------------------------------------------

def _read_sample_list(path: Path) -> List[str]:
    samples: List[str] = []
    with path.open() as handle:
        for line in handle:
            s = line.strip()
            if s and not s.startswith("#"):
                if "HCPT" in s:
                    s = s[:8]
                samples.append(s)
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def _build_dense(
    df: pd.DataFrame,
    value_cols: List[str],
    sample_index: Dict[str, int],
    chr_index: Dict[str, int],
) -> Tuple[List[Tuple[float, float]], List[np.ndarray]]:
    """Pivot a long grid-search table into dense [n_combo, n_chr, n_sample] arrays.

    Returns the ordered combo list (sorted by threshold then recall) and one
    dense float64 array per requested value column.
    """
    thr = df["threshold"].astype(float).to_numpy()
    rec = df["recall"].astype(float).to_numpy()
    combos = sorted({(float(t), float(r)) for t, r in zip(thr, rec)})
    combo_index = {c: i for i, c in enumerate(combos)}

    ci = np.fromiter((combo_index[(t, r)] for t, r in zip(thr, rec)),
                     dtype=np.int64, count=len(df))
    hi = df["chr"].astype(str).map(chr_index).to_numpy()
    si = df["sample"].astype(str).map(sample_index).to_numpy()

    n_combo = len(combos)
    n_chr = len(chr_index)
    n_samp = len(sample_index)

    keep = ~(np.isnan(hi.astype(np.float64)) | np.isnan(si.astype(np.float64)))
    if not keep.all():
        ci, hi, si = ci[keep], hi[keep], si[keep]
    flat = (ci.astype(np.int64) * n_chr + hi.astype(np.int64)) * n_samp + si.astype(np.int64)

    arrays: List[np.ndarray] = []
    for col in value_cols:
        vals = df[col].to_numpy(dtype=np.float64)
        if not keep.all():
            vals = vals[keep]
        arr = np.full(n_combo * n_chr * n_samp, np.nan, dtype=np.float64)
        arr[flat] = vals
        arrays.append(arr.reshape(n_combo, n_chr, n_samp))
    return combos, arrays


# ---------------------------------------------------------------------------
# Reference-normalized score computation
# ---------------------------------------------------------------------------

def _ref_mean_std(values: np.ndarray, ref_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Mean/std over the reference sample axis for [combo, chr, sample] arrays."""
    sub = values[:, :, ref_idx]
    with np.errstate(invalid="ignore"):
        means = np.nanmean(sub, axis=2)
        stds = np.nanstd(sub, axis=2, ddof=0)
    means = np.where(np.isfinite(means), means, 0.0)
    return means, stds


def compute_episcore(
    hypo_z: np.ndarray,
    hyper_z: np.ndarray,
    hypo_cnt: np.ndarray,
    hyper_cnt: np.ndarray,
    ref_idx: np.ndarray,
) -> np.ndarray:
    """episcore s_inter for every [combo, chr, sample]."""
    hypo_mean, hypo_std = _ref_mean_std(hypo_z, ref_idx)
    hyper_mean, hyper_std = _ref_mean_std(hyper_z, ref_idx)
    hypo_std_safe = np.where(hypo_std > 0, hypo_std, np.nan)[:, :, None]
    hyper_std_safe = np.where(hyper_std > 0, hyper_std, np.nan)[:, :, None]

    with np.errstate(divide="ignore", invalid="ignore"):
        hypo_z_inter = (hypo_z - hypo_mean[:, :, None]) / hypo_std_safe
        hyper_z_inter = (hyper_z - hyper_mean[:, :, None]) / hyper_std_safe

    w_hypo = np.sqrt(np.nan_to_num(hypo_cnt, nan=0.0))
    w_hyper = np.sqrt(np.nan_to_num(hyper_cnt, nan=0.0))
    total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)
    total_w = np.where(total_w > 0, total_w, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        s_inter = (hyper_z_inter * w_hyper - hypo_z_inter * w_hypo) / total_w
    return np.where(np.isnan(s_inter), 0.0, s_inter)


def compute_zscore(percentage: np.ndarray, ref_idx: np.ndarray) -> np.ndarray:
    """Reference z-score for chromosome percentage, per [combo, chr, sample]."""
    mean, std = _ref_mean_std(percentage, ref_idx)
    std_safe = np.where(std > 0, std, np.nan)[:, :, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        return (percentage - mean[:, :, None]) / std_safe


# ---------------------------------------------------------------------------
# Method 3 grid search (per chromosome best combo on the dev analyze set)
# ---------------------------------------------------------------------------

def _per_chr_metrics(
    scores_cs: np.ndarray,        # (n_combo, n_dev)
    combos: List[Tuple[float, float]],
    is_target: np.ndarray,        # (n_dev,) bool
    is_normal: np.ndarray,        # (n_dev,) bool
    cutoff: float,
) -> pd.DataFrame:
    pred = scores_cs > cutoff
    tp = (pred & is_target).sum(axis=1)
    fp = (pred & is_normal).sum(axis=1)
    fn = ((~pred) & is_target).sum(axis=1)
    tn = ((~pred) & is_normal).sum(axis=1)

    tp_scores = np.where(pred & is_target, scores_cs, np.nan)
    tn_scores = np.where((~pred) & is_normal, scores_cs, np.nan)
    tp_min = np.full(scores_cs.shape[0], np.nan)
    tn_max = np.full(scores_cs.shape[0], np.nan)
    valid_tp = ~np.all(np.isnan(tp_scores), axis=1)
    valid_tn = ~np.all(np.isnan(tn_scores), axis=1)
    with np.errstate(invalid="ignore"):
        if valid_tp.any():
            tp_min[valid_tp] = np.nanmin(tp_scores[valid_tp], axis=1)
        if valid_tn.any():
            tn_max[valid_tn] = np.nanmax(tn_scores[valid_tn], axis=1)
    sep_score = tp_min - tn_max

    fp_ext_mask = ((scores_cs > cutoff) | (scores_cs < -cutoff)) & is_normal
    fp_ext = fp_ext_mask.sum(axis=1)
    normal_scores = np.where(is_normal, scores_cs, np.nan)
    with np.errstate(invalid="ignore"):
        normal_var = np.nanvar(normal_scores, axis=1, ddof=0)

    tpf, fpf, fnf, tnf = tp.astype(float), fp.astype(float), fn.astype(float), tn.astype(float)
    denom = np.sqrt((tpf + fpf) * (tpf + fnf) * (tnf + fpf) * (tnf + fnf))
    with np.errstate(divide="ignore", invalid="ignore"):
        mcc = np.where(denom > 0, (tpf * tnf - fpf * fnf) / denom, 0.0)

    return pd.DataFrame(
        {
            "threshold": [c[0] for c in combos],
            "recall": [c[1] for c in combos],
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "mcc": mcc, "sep_score": sep_score,
            "fp_ext": fp_ext, "normal_var": normal_var,
        }
    )


def _pick_chr_combo(metrics: pd.DataFrame, has_target: bool, min_recall: float) -> Optional[Tuple[float, float]]:
    sub = metrics[(metrics["recall"] >= min_recall) & (metrics["recall"] <= MAX_RECALL)]
    if sub.empty:
        return None
    use_case_b = has_target and bool((sub["tp"] > 0).any())
    if use_case_b:
        by, asc = ["mcc", "sep_score", "recall"], [False, False, False]
    else:
        by, asc = ["fp_ext", "normal_var", "recall"], [True, True, False]
    best = sub.sort_values(by=by, ascending=asc, kind="mergesort").iloc[0]
    return float(best["threshold"]), float(best["recall"])


def _overall_confusion(
    score_all: np.ndarray,                       # (n_combo, n_chr, n_sample)
    combo_index: Dict[Tuple[float, float], int],
    chr_combos: Dict[str, Tuple[float, float]],
    eval_idx: np.ndarray,                        # sample indices to evaluate
    is_normal_eval: np.ndarray,                  # (len(eval_idx),) bool
    cutoff: float,
) -> Dict[str, float]:
    if eval_idx.size == 0:
        return {"tp": 0, "tn": 0, "fp": 0, "fn": 0, "mcc": 0.0}
    cols = []
    for hi, chrom in enumerate(CHR_LIST):
        ci = combo_index[chr_combos[chrom]]
        cols.append(score_all[ci, hi, eval_idx])
    z_wide = np.vstack(cols).T  # (n_eval, n_chr)
    any_pos = (z_wide > cutoff).any(axis=1)
    is_normal = is_normal_eval
    fp = int((is_normal & any_pos).sum())
    tn = int((is_normal & ~any_pos).sum())
    tp = int((~is_normal & any_pos).sum())
    fn = int((~is_normal & ~any_pos).sum())
    denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0
    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn, "mcc": mcc}


def method3_grid_search(
    score_all: np.ndarray,
    combos: List[Tuple[float, float]],
    grid_idx: np.ndarray,
    grid_labels: np.ndarray,
    cutoff: float,
) -> Tuple[Dict[str, Tuple[float, float]], Dict[str, bool], float]:
    """Pick per-chr best combo by sweeping min_recall on the grid-search set.

    Returns (best_combo_per_chr, has_target_per_chr, best_min_recall).
    """
    combo_index = {c: i for i, c in enumerate(combos)}
    per_chr_metrics: Dict[str, pd.DataFrame] = {}
    has_target: Dict[str, bool] = {}
    grid_is_normal = grid_labels == "Normal"

    for hi, chrom in enumerate(CHR_LIST):
        target_label = chrom.replace("chr", "T")
        is_target = grid_labels == target_label
        has_target[chrom] = bool(is_target.any())
        if has_target[chrom]:
            keep = is_target | grid_is_normal
        else:
            keep = grid_is_normal
        scores_cs = score_all[:, hi, :][:, grid_idx][:, keep]
        per_chr_metrics[chrom] = _per_chr_metrics(
            scores_cs, combos, is_target[keep], grid_is_normal[keep], cutoff
        )

    min_recalls = [i / 100 for i in range(1, 100)]
    best_min_recall = None
    best_mcc = -np.inf
    for min_recall in min_recalls:
        chr_combos: Dict[str, Tuple[float, float]] = {}
        valid = True
        for chrom in CHR_LIST:
            combo = _pick_chr_combo(per_chr_metrics[chrom], has_target[chrom], min_recall)
            if combo is None:
                valid = False
                break
            chr_combos[chrom] = combo
        if not valid:
            continue
        grid_conf = _overall_confusion(
            score_all, combo_index, chr_combos, grid_idx, grid_is_normal, cutoff
        )
        # method 3 tie-break: higher grid MCC, then larger min_recall.
        if grid_conf["mcc"] > best_mcc or (
            np.isclose(grid_conf["mcc"], best_mcc) and (best_min_recall is None or min_recall > best_min_recall)
        ):
            best_mcc = grid_conf["mcc"]
            best_min_recall = min_recall

    if best_min_recall is None:
        raise RuntimeError("method3_grid_search found no valid min_recall")

    best_combos = {
        chrom: _pick_chr_combo(per_chr_metrics[chrom], has_target[chrom], best_min_recall)
        for chrom in CHR_LIST
    }
    return best_combos, has_target, best_min_recall


def fixed_combo_grid_search(
    score_all: np.ndarray,
    combos: List[Tuple[float, float]],
    grid_idx: np.ndarray,
    grid_labels: np.ndarray,
    cutoff: float,
) -> Tuple[Dict[str, Tuple[float, float]], Dict[str, bool], float]:
    """Pick one (threshold, recall) combo shared by all chromosomes.

    Each combo is scored by overall sample-level MCC on the grid-search set
    (positive when any chromosome exceeds ``cutoff``). Returns the same combo
    for every chromosome plus ``has_target`` flags per chromosome.
    """
    combo_index = {c: i for i, c in enumerate(combos)}
    grid_is_normal = grid_labels == "Normal"
    has_target = {
        chrom: bool((grid_labels == chrom.replace("chr", "T")).any())
        for chrom in CHR_LIST
    }

    best_combo: Optional[Tuple[float, float]] = None
    best_mcc = -np.inf
    for combo in combos:
        chr_combos = {chrom: combo for chrom in CHR_LIST}
        conf = _overall_confusion(
            score_all, combo_index, chr_combos, grid_idx, grid_is_normal, cutoff
        )
        mcc = conf["mcc"]
        better = False
        if mcc > best_mcc:
            better = True
        elif np.isclose(mcc, best_mcc):
            if best_combo is None:
                better = True
            elif combo[1] > best_combo[1]:
                better = True
            elif np.isclose(combo[1], best_combo[1]) and combo[0] < best_combo[0]:
                better = True
        if better:
            best_mcc = mcc
            best_combo = combo

    if best_combo is None:
        raise RuntimeError("fixed_combo_grid_search found no valid combo")

    best_combos = {chrom: best_combo for chrom in CHR_LIST}
    return best_combos, has_target, best_combo[1]


def _resolve_mode(mode: str) -> Tuple[str, bool]:
    """Return (sample_split, use_fixed_combo)."""
    if mode == "all":
        return "all", False
    if mode == "dev_test_split":
        return "dev_test_split", False
    if mode == "fix_combo_all":
        return "all", True
    if mode == "fix_combo_split":
        return "dev_test_split", True
    raise ValueError(f"Unknown mode: {mode}")


# ---------------------------------------------------------------------------
# Per-repeat driver
# ---------------------------------------------------------------------------

def run_repeat(
    repeat_index: int,
    ref_idx: np.ndarray,
    samples: List[str],
    sample_index: Dict[str, int],
    set_arr: np.ndarray,
    label_arr: np.ndarray,
    ez_idx: np.ndarray,
    ep_combos: List[Tuple[float, float]],
    ep_arrays: List[np.ndarray],
    z_combos: List[Tuple[float, float]],
    z_array: np.ndarray,
    cutoff: float,
    mode: str,
    analyze_ff_mask: np.ndarray,
    out_dir: Path,
) -> Dict[str, object]:
    n_sample = len(samples)
    ref_mask = np.zeros(n_sample, dtype=bool)
    ref_mask[ref_idx] = True
    analyze_mask = ~ref_mask & analyze_ff_mask

    analyze_idx = np.flatnonzero(analyze_mask)
    sample_split, use_fixed_combo = _resolve_mode(mode)
    if sample_split == "all":
        grid_idx = analyze_idx
        eval_sets: List[Tuple[str, np.ndarray]] = [("all", analyze_idx)]
    else:
        dev_idx = np.flatnonzero(analyze_mask & (set_arr == "dev"))
        test_idx = np.flatnonzero(analyze_mask & (set_arr == "test"))
        grid_idx = dev_idx
        eval_sets = [("dev", dev_idx), ("test", test_idx)]

    if grid_idx.size == 0:
        raise ValueError(f"mode={mode}: no analyze samples available for grid search")

    # --- episcore / zscore scores for all combos ---------------------------
    episcore_all = compute_episcore(ep_arrays[0], ep_arrays[1], ep_arrays[2], ep_arrays[3], ref_idx)
    zscore_all = compute_zscore(z_array, ref_idx)

    grid_labels = label_arr[grid_idx]
    grid_search = fixed_combo_grid_search if use_fixed_combo else method3_grid_search
    ep_best, ep_has_target, ep_min_recall = grid_search(
        episcore_all, ep_combos, grid_idx, grid_labels, cutoff
    )
    z_best, z_has_target, z_min_recall = grid_search(
        zscore_all, z_combos, grid_idx, grid_labels, cutoff
    )

    ep_combo_index = {c: i for i, c in enumerate(ep_combos)}
    z_combo_index = {c: i for i, c in enumerate(z_combos)}

    # --- ezscore: z-normalize(zscore + episcore) over ezscore reference ----
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

    # --- confusion / MCC for each score x set ------------------------------
    def _conf_from_best(best_mat: np.ndarray, eval_idx: np.ndarray, is_normal: np.ndarray) -> Dict[str, float]:
        if eval_idx.size == 0:
            return {"tp": 0, "tn": 0, "fp": 0, "fn": 0, "mcc": 0.0}
        sub = best_mat[:, eval_idx]  # (n_chr, n_eval)
        any_pos = (sub > cutoff).any(axis=0)
        fp = int((is_normal & any_pos).sum())
        tn = int((is_normal & ~any_pos).sum())
        tp = int((~is_normal & any_pos).sum())
        fn = int((~is_normal & ~any_pos).sum())
        denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
        mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0
        return {"tp": tp, "tn": tn, "fp": fp, "fn": fn, "mcc": mcc}

    metrics_rows = []
    metric_lookup: Dict[str, Dict[str, float]] = {}
    for score_name, best_mat in (("episcore", episcore_best), ("zscore", zscore_best), ("ezscore", ez_score_all)):
        for set_name, eidx in eval_sets:
            isn = label_arr[eidx] == "Normal"
            conf = _conf_from_best(best_mat, eidx, isn)
            row = {"score": score_name, "set": set_name, **conf}
            metrics_rows.append(row)
            metric_lookup[f"{score_name}_{set_name}"] = conf

    # --- write per-repeat outputs ------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_samples = [samples[i] for i in ref_idx]
    (out_dir / "reference_samples.txt").write_text("\n".join(ref_samples) + "\n")

    pd.DataFrame(
        {
            "chr": CHR_LIST,
            "threshold": [ep_best[c][0] for c in CHR_LIST],
            "recall": [ep_best[c][1] for c in CHR_LIST],
            "has_target": [ep_has_target[c] for c in CHR_LIST],
            "min_recall": [ep_min_recall] * len(CHR_LIST),
        }
    ).to_csv(out_dir / "best_combo_episcore.csv", index=False)
    pd.DataFrame(
        {
            "chr": CHR_LIST,
            "threshold": [z_best[c][0] for c in CHR_LIST],
            "recall": [z_best[c][1] for c in CHR_LIST],
            "has_target": [z_has_target[c] for c in CHR_LIST],
            "min_recall": [z_min_recall] * len(CHR_LIST),
        }
    ).to_csv(out_dir / "best_combo_zscore.csv", index=False)

    metrics_df = pd.DataFrame(metrics_rows)[["score", "set", "mcc", "tp", "tn", "fp", "fn"]]
    metrics_df.insert(0, "repeat_index", repeat_index)
    metrics_df.to_csv(out_dir / "metrics.tsv", sep="\t", index=False, float_format="%.6f")

    # Per analyze-sample scores under the repeat's best combos.
    score_data = {
        "sample": [samples[i] for i in analyze_idx],
        "set": set_arr[analyze_idx],
        "label": label_arr[analyze_idx],
    }
    for hi, chrom in enumerate(CHR_LIST):
        num = chrom.removeprefix("chr")
        score_data[f"episcore_chr{num}"] = episcore_best[hi, analyze_idx]
        score_data[f"zscore_chr{num}"] = zscore_best[hi, analyze_idx]
        score_data[f"ezscore_chr{num}"] = ez_score_all[hi, analyze_idx]
    pd.DataFrame(score_data).to_csv(out_dir / "scores.tsv", sep="\t", index=False, float_format="%.6f")

    summary = {
        "repeat_index": repeat_index,
        "reference_list": ",".join(ref_samples),
        "episcore_min_recall": ep_min_recall,
        "zscore_min_recall": z_min_recall,
    }
    for key, conf in metric_lookup.items():
        for m in ("mcc", "tp", "tn", "fp", "fn"):
            summary[f"{key}_{m}"] = conf[m]
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Dir with episcore_grid_search.parquet, zscore_grid_search.parquet, meta.csv, ref_pool_samples.txt, ezscore_ref_samples.txt")
@click.option("--output-base", required=True, type=click.Path(file_okay=False),
              help="Output base; repeat dirs land under randomly_select_ref_40/")
@click.option("--total-repeats", default=100, show_default=True, type=int)
@click.option("--repeat-start", default=0, show_default=True, type=int)
@click.option("--repeat-end", default=None, type=int, help="Exclusive. Default: total-repeats")
@click.option("--ref-n", default=40, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--cutoff", default=DEFAULT_CUTOFF, show_default=True, type=float,
              help="Positive-call z-score threshold")
@click.option("--ref-pool-file", default="ref_pool_samples.txt", show_default=True,
              help="File name (under input-dir) or path of the candidate reference pool")
@click.option("--ezscore-ref-file", default="ezscore_ref_samples.txt", show_default=True,
              help="File name (under input-dir) or path of the ezscore reference list")
@click.option("--mode", default="dev_test_split", show_default=True,
              type=click.Choice(["dev_test_split", "all", "fix_combo_all", "fix_combo_split"]),
              help="dev_test_split/all: per-chr combo (method 3); "
                   "fix_combo_all/fix_combo_split: one combo for all chromosomes")
@click.option("--min-ff", default=0.0, show_default=True, type=float,
              help="Keep samples with ff_before_mq > min_ff")
def main(
    input_dir: str,
    output_base: str,
    total_repeats: int,
    repeat_start: int,
    repeat_end: Optional[int],
    ref_n: int,
    seed: int,
    cutoff: float,
    ref_pool_file: str,
    ezscore_ref_file: str,
    mode: str,
    min_ff: float,
) -> None:
    """Run random ref-40 episcore/zscore/ezscore grid search repeats."""
    input_path = Path(input_dir)
    out_root = Path(output_base) / "randomly_select_ref_40"
    out_root.mkdir(parents=True, exist_ok=True)

    if repeat_end is None:
        repeat_end = total_repeats
    if repeat_start < 0 or repeat_end > total_repeats or repeat_end <= repeat_start:
        raise click.ClickException(
            f"Repeat slice [{repeat_start}, {repeat_end}) must lie within [0, {total_repeats})"
        )

    console.rule("[bold blue]Random ref-40 grid search")
    console.print(f"  Input dir    : {input_path}")
    console.print(f"  Output root  : {out_root}")
    console.print(f"  Repeat range : [{repeat_start}, {repeat_end}) of {total_repeats}")
    console.print(f"  ref-n / seed : {ref_n} / {seed}")
    console.print(f"  cutoff       : {cutoff}")
    console.print(f"  mode         : {mode}")
    console.print(f"  min-ff       : {min_ff}")

    # --- meta ---------------------------------------------------------------
    meta = pd.read_csv(input_path / "meta.csv")
    for col in ("sample", "set", "label", "ff_before_mq"):
        if col not in meta.columns:
            raise click.ClickException(f"meta.csv missing column: {col}")
    meta = meta.drop_duplicates("sample", keep="first").copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")

    pool_path = Path(ref_pool_file)
    if not pool_path.is_absolute() and not pool_path.exists():
        pool_path = input_path / ref_pool_file
    pool_samples = set(_read_sample_list(pool_path))

    ez_path = Path(ezscore_ref_file)
    if not ez_path.is_absolute() and not ez_path.exists():
        ez_path = input_path / ezscore_ref_file
    ez_samples = set(_read_sample_list(ez_path))

    # --- grid-search parquets ----------------------------------------------
    console.print("[cyan]Loading parquets ...[/cyan]")
    ep_df = pd.read_parquet(input_path / "episcore_grid_search.parquet")
    z_df = pd.read_parquet(input_path / "zscore_grid_search.parquet")

    ep_samples = set(ep_df["sample"].astype(str).unique())
    z_samples = set(z_df["sample"].astype(str).unique())
    meta_samples = set(meta["sample"])

    meta_ff = meta.set_index("sample")
    ff_pass = set(meta_ff.index[meta_ff["ff_before_mq"] > min_ff].astype(str))

    # Universe: ff-passing samples + ref pool + ezscore ref (for normalization).
    universe = sorted(
        (meta_samples & ep_samples & z_samples & ff_pass)
        | (pool_samples & meta_samples & ep_samples & z_samples)
        | (ez_samples & meta_samples & ep_samples & z_samples)
    )
    if not universe:
        raise click.ClickException("No samples remain for scoring after filters")
    samples = universe
    keep_samples = set(samples)
    ep_df = ep_df[ep_df["sample"].astype(str).isin(keep_samples)]
    z_df = z_df[z_df["sample"].astype(str).isin(keep_samples)]
    sample_index = {s: i for i, s in enumerate(samples)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}

    meta_idx = meta.set_index("sample").reindex(samples)
    set_arr = meta_idx["set"].astype(str).to_numpy()
    label_arr = meta_idx["label"].astype(str).to_numpy()
    ff_arr = pd.to_numeric(meta_idx["ff_before_mq"], errors="coerce").to_numpy()
    analyze_ff_mask = ff_arr > min_ff

    n_analyze_ff = int(analyze_ff_mask.sum())
    console.print(f"  samples (universe) : {len(samples)}")
    console.print(f"  analyze (ff>{min_ff:g}) : {n_analyze_ff}")
    console.print("[cyan]Building dense arrays ...[/cyan]")
    ep_combos, ep_arrays = _build_dense(
        ep_df, ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
        sample_index, chr_index,
    )
    z_combos, z_arrays = _build_dense(z_df, ["percentage"], sample_index, chr_index)
    z_array = z_arrays[0]
    console.print(f"  episcore combos  : {len(ep_combos)}")
    console.print(f"  zscore combos    : {len(z_combos)}")

    # --- reference pool + ezscore reference --------------------------------
    pool_idx = np.array([sample_index[s] for s in pool_samples if s in sample_index], dtype=np.int64)
    if pool_idx.size < ref_n:
        raise click.ClickException(f"Reference pool ({pool_idx.size}) smaller than ref-n ({ref_n})")

    ez_idx = np.array([sample_index[s] for s in ez_samples if s in sample_index], dtype=np.int64)
    if ez_idx.size == 0:
        raise click.ClickException("No ezscore reference samples found in data")
    console.print(f"  reference pool   : {pool_idx.size}")
    console.print(f"  ezscore ref      : {ez_idx.size}")

    rng = np.random.default_rng(seed)
    all_draws = generate_unique_ref_draws(pool_idx.size, ref_n, total_repeats, rng)
    draws = all_draws[repeat_start:repeat_end]

    summaries: List[Dict[str, object]] = []
    failures: List[str] = []
    for offset, local_draw in enumerate(draws):
        repeat_index = repeat_start + offset
        ref_idx = pool_idx[np.asarray(local_draw, dtype=np.int64)]
        repeat_dir = out_root / f"repeat_{repeat_index}"
        try:
            summary = run_repeat(
                repeat_index, ref_idx, samples, sample_index, set_arr, label_arr, ez_idx,
                ep_combos, ep_arrays, z_combos, z_array, cutoff, mode, analyze_ff_mask, repeat_dir,
            )
            summaries.append(summary)
            sample_split, _ = _resolve_mode(mode)
            if sample_split == "all":
                console.print(
                    f"  repeat {repeat_index}: ez all MCC={summary['ezscore_all_mcc']:.3f}"
                )
            else:
                console.print(
                    f"  repeat {repeat_index}: "
                    f"ez dev MCC={summary['ezscore_dev_mcc']:.3f} "
                    f"test MCC={summary['ezscore_test_mcc']:.3f}"
                )
        except Exception as exc:  # noqa: BLE001
            failures.append(f"repeat_{repeat_index}: {exc}")
            console.print_exception()

    if summaries:
        manifest = pd.DataFrame(summaries)
        manifest_path = out_root / f"metrics_manifest_{repeat_start}_{repeat_end}.tsv"
        manifest.to_csv(manifest_path, sep="\t", index=False, float_format="%.6f")
        console.print(f"[green]OK[/green] Wrote {manifest_path}")

    if failures:
        for msg in failures[:10]:
            console.print(f"[red]{msg}[/red]")
        raise click.ClickException(f"{len(failures)} repeats failed")

    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
