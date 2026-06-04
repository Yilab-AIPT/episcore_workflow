#!/usr/bin/env python3
"""
Enlarged-reference per-chromosome Z-score recomputation.

Goal:
    Test how the per-chromosome ``s_inter`` values respond to growing the
    reference set with additional dev-set Normal samples. For each enlarged
    reference size ``N`` (and many random draws per ``N``) we recompute
    ``s_inter`` for a *fixed* evaluation set and report a per-run MCC.

Inputs:
    --combo-dir : directory containing ``_analyze_zscore.tsv.gz`` and
                  ``_reference_zscore.tsv.gz`` produced by ``beta_to_zscore.py``
                  for the fixed (threshold, recall) combo.
    --meta-csv  : sample-level metadata. Must provide
                  ``sample, label, ref_type, set, ff_before_mq``.

Reference candidate pool
    ``set == "dev"`` AND ``label == "Normal"``. For each ``N`` from
    ``--n-min`` (default 10) to the pool size, ``--runs`` random draws
    (default 10000) are generated. Combinations are guaranteed to be
    unique within a given ``N``; when ``C(pool_size, N) <= runs`` the
    enumeration is exhaustive (so ``N == pool_size`` uses exactly one run).

Evaluation set (fixed across all runs)
    All samples that are *not* in the candidate pool, i.e.
    ``set == "test"`` OR (``set == "dev"`` AND ``label != "Normal"``).
    Holding the eval set constant makes MCC values directly comparable
    across runs and across ``N``. ``_classify`` further restricts the
    contribution to rows whose label is either ``Normal`` or starts with
    ``T`` (e.g. T21).

Modes (``--mode``):

    all (default)
        Same as the original workflow: one dev+Normal pool, fixed combined
        eval set (test + dev non-Normal), single MCC per run.

    isolated
        Isolate dev and test evaluation. For each of 100 random 50-sample
        reference pools drawn from dev Normal, draw 50 unique reference sets
        per ``N`` in [5, 45], compute ``MCC_dev`` on remaining dev samples
        (outside the pool) and ``MCC_test`` on the full test set. Maintain a
        running top-100 list (by ``MCC_dev + MCC_test``) for ``N`` in [35, 45].

Outputs (under ``--output-dir`` or ``--output-base/<output-subdir>/``):

    report.csv

    all mode columns:
        ``ref_n, run_index, MCC, TP, TN, FP, FN``

    isolated mode columns:
        ``pool_run_index, ref_run_index, ref_n, MCC_dev, MCC_test,
         TP_dev, TN_dev, FP_dev, FN_dev, TP_test, TN_test, FP_test, FN_test``

    top_100_reference_list.csv (isolated mode only):
        ``ref_n, MCC_dev, MCC_test, reference_list``

The classification rule for the MCC counts is:
    Positive : ``label`` starts with ``T`` (i.e. trisomy)
    Negative : ``label == "Normal"``
    Predicted positive : ``any(chr_s_inter > --s-inter-cutoff)``
    Predicted negative : otherwise
"""

from __future__ import annotations

import itertools
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

console = Console()


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

# Only the columns required for recomputing s_inter and classifying.
_NEEDED_SUFFIXES: Tuple[str, ...] = (
    "hypo_z_intra",
    "hyper_z_intra",
    "hypo_cpgs_count",
    "hyper_cpgs_count",
)


def _load_combo_table(path: Path, keep_cols: List[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing combo output: {path}")
    return pd.read_csv(path, sep="\t", compression="gzip", usecols=keep_cols)


def _stack_per_chr(
    df: pd.DataFrame,
    chr_list: List[str],
    suffix: str,
    *,
    dtype: type = np.float64,
) -> np.ndarray:
    """Build an ``(n_samples, n_chr)`` matrix from ``chr{n}_{suffix}`` columns."""
    cols = [f"{c}_{suffix}" for c in chr_list]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Columns missing: {missing}")
    arr = df.loc[:, cols].to_numpy(copy=False)
    if np.issubdtype(dtype, np.integer):
        arr = np.nan_to_num(arr, nan=0).astype(dtype, copy=False)
    else:
        arr = arr.astype(dtype, copy=False)
    return arr


def _load_combined(
    combo_dir: Path,
    chr_list: List[str],
) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Concat reference + analyze TSVs and return only the metric arrays we need."""
    keep = ["sample"] + [
        f"{c}_{suf}" for c in chr_list for suf in _NEEDED_SUFFIXES
    ]
    ref_df = _load_combo_table(combo_dir / "_reference_zscore.tsv.gz", keep)
    an_df = _load_combo_table(combo_dir / "_analyze_zscore.tsv.gz", keep)

    combined = pd.concat([ref_df, an_df], axis=0, ignore_index=True)
    del ref_df, an_df
    if combined["sample"].duplicated().any():
        dups = combined["sample"][combined["sample"].duplicated()].tolist()
        raise ValueError(
            f"Duplicate sample IDs across reference + analyze: {dups[:5]} "
            f"(total {len(dups)})"
        )

    samples = combined["sample"].astype(str).tolist()
    arrays = {
        "hypo_z_intra": _stack_per_chr(combined, chr_list, "hypo_z_intra"),
        "hyper_z_intra": _stack_per_chr(combined, chr_list, "hyper_z_intra"),
        "hypo_counts": _stack_per_chr(
            combined, chr_list, "hypo_cpgs_count", dtype=np.int64,
        ),
        "hyper_counts": _stack_per_chr(
            combined, chr_list, "hyper_cpgs_count", dtype=np.int64,
        ),
    }
    return samples, arrays


def _load_meta(meta_csv: str, samples: List[str]) -> pd.DataFrame:
    """Load meta and align to ``samples`` (one row per sample, in order)."""
    needed = ["sample", "label", "ref_type", "set", "ff_before_mq"]
    meta = pd.read_csv(meta_csv, usecols=needed)
    meta = meta.drop_duplicates(subset="sample", keep="first")
    aligned = meta.set_index("sample").reindex(samples).reset_index()
    missing = aligned["label"].isna() & aligned["ref_type"].isna() & aligned["set"].isna()
    if missing.any():
        miss_samples = aligned.loc[missing, "sample"].tolist()
        console.print(
            f"[yellow]Warning[/yellow] {len(miss_samples)} samples have no row "
            f"in meta (e.g. {miss_samples[:3]}); they will be excluded from "
            "MCC and from the reference candidate pool."
        )
    return aligned


# ---------------------------------------------------------------------------
# Combo generation (one unique sample-set per "run")
# ---------------------------------------------------------------------------

def _generate_unique_combos(
    pool_size: int,
    N: int,
    runs: int,
    rng: np.random.Generator,
) -> List[np.ndarray]:
    """Return up to ``runs`` *unique* sorted combinations of ``N`` indices.

    Each returned array contains indices in ``[0, pool_size)``.

    If ``C(pool_size, N) <= runs`` the full set of combinations is enumerated
    (so e.g. ``N == pool_size`` yields exactly one combo). Otherwise unique
    random draws are produced via rejection sampling.
    """
    if N <= 0 or N > pool_size:
        raise ValueError(f"Invalid N={N} for pool_size={pool_size}.")
    total = math.comb(int(pool_size), int(N))
    if total <= runs:
        return [
            np.fromiter(c, dtype=np.int64, count=N)
            for c in itertools.combinations(range(pool_size), int(N))
        ]
    seen: set = set()
    out: List[np.ndarray] = []
    max_attempts = runs * 20 + 1000
    attempts = 0
    while len(out) < runs and attempts < max_attempts:
        attempts += 1
        choice = rng.choice(pool_size, size=N, replace=False)
        choice.sort()
        key = choice.tobytes()
        if key in seen:
            continue
        seen.add(key)
        out.append(choice.astype(np.int64, copy=False))
    if len(out) < runs:
        raise RuntimeError(
            f"Could only generate {len(out)}/{runs} unique combos for N={N} "
            f"after {attempts} attempts (pool_size={pool_size})."
        )
    return out


# ---------------------------------------------------------------------------
# Vectorized batch evaluation
# ---------------------------------------------------------------------------

def _evaluate_batch(
    chosen_batch: np.ndarray,            # (B, N) int indices into pool arrays
    pool_hypo_z: np.ndarray,             # (P, C)
    pool_hyper_z: np.ndarray,            # (P, C)
    eval_hypo_z: np.ndarray,             # (E, C)
    eval_hyper_z: np.ndarray,            # (E, C)
    eval_w_hypo: np.ndarray,             # (E, C)
    eval_w_hyper: np.ndarray,            # (E, C)
    eval_total_w: np.ndarray,            # (E, C)
    eval_is_normal: np.ndarray,          # (E,)  bool
    eval_is_trisomy: np.ndarray,         # (E,)  bool
    s_inter_cutoff: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized evaluation of B runs sharing the same N.

    Returns (TP, TN, FP, FN, MCC) each shaped (B,).
    """
    ref_hypo = pool_hypo_z[chosen_batch]    # (B, N, C)
    ref_hyper = pool_hyper_z[chosen_batch]  # (B, N, C)

    with np.errstate(invalid="ignore"):
        hypo_means = np.nanmean(ref_hypo, axis=1)   # (B, C)
        hypo_stds = np.nanstd(ref_hypo, axis=1, ddof=0)
        hyper_means = np.nanmean(ref_hyper, axis=1)
        hyper_stds = np.nanstd(ref_hyper, axis=1, ddof=0)
    del ref_hypo, ref_hyper

    hypo_means = np.where(np.isfinite(hypo_means), hypo_means, 0.0)
    hyper_means = np.where(np.isfinite(hyper_means), hyper_means, 0.0)
    hypo_std_safe = np.where(hypo_stds > 0, hypo_stds, np.nan)
    hyper_std_safe = np.where(hyper_stds > 0, hyper_stds, np.nan)

    eval_hypo_z_b = eval_hypo_z[None, :, :]    # (1, E, C)
    eval_hyper_z_b = eval_hyper_z[None, :, :]
    hypo_means_b = hypo_means[:, None, :]      # (B, 1, C)
    hyper_means_b = hyper_means[:, None, :]
    hypo_std_b = hypo_std_safe[:, None, :]
    hyper_std_b = hyper_std_safe[:, None, :]

    with np.errstate(divide="ignore", invalid="ignore"):
        hypo_z_inter = (eval_hypo_z_b - hypo_means_b) / hypo_std_b   # (B, E, C)
        hyper_z_inter = (eval_hyper_z_b - hyper_means_b) / hyper_std_b
        s_inter = (hyper_z_inter * eval_w_hyper - hypo_z_inter * eval_w_hypo) / eval_total_w
    del hypo_z_inter, hyper_z_inter
    s_inter = np.where(np.isnan(s_inter), 0.0, s_inter)             # (B, E, C)

    with np.errstate(invalid="ignore"):
        s_inter_max = np.max(s_inter, axis=2)                       # (B, E)
    del s_inter
    any_pos = s_inter_max > s_inter_cutoff                          # (B, E) bool

    is_normal_b = eval_is_normal[None, :]
    is_trisomy_b = eval_is_trisomy[None, :]
    tp = np.sum(any_pos & is_trisomy_b, axis=1).astype(np.int64)
    fn = np.sum((~any_pos) & is_trisomy_b, axis=1).astype(np.int64)
    fp = np.sum(any_pos & is_normal_b, axis=1).astype(np.int64)
    tn = np.sum((~any_pos) & is_normal_b, axis=1).astype(np.int64)

    num = tp.astype(np.float64) * tn - fp.astype(np.float64) * fn
    denom_sq = (
        (tp + fp).astype(np.float64)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        mcc = np.where(denom_sq > 0, num / np.sqrt(denom_sq), np.nan)
    return tp, tn, fp, fn, mcc


def _mcc_from_counts(
    tp: np.ndarray,
    tn: np.ndarray,
    fp: np.ndarray,
    fn: np.ndarray,
) -> np.ndarray:
    """Return MCC shaped like ``tp`` (supports batch dimension)."""
    num = tp.astype(np.float64) * tn - fp.astype(np.float64) * fn
    denom_sq = (
        (tp + fp).astype(np.float64)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(denom_sq > 0, num / np.sqrt(denom_sq), np.nan)


def _confusion_for_eval(
    s_inter_max: np.ndarray,          # (B, E)
    eval_is_normal: np.ndarray,       # (E,)
    eval_is_trisomy: np.ndarray,      # (E,)
    s_inter_cutoff: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (TP, TN, FP, FN) each shaped (B,)."""
    any_pos = s_inter_max > s_inter_cutoff
    is_normal_b = eval_is_normal[None, :]
    is_trisomy_b = eval_is_trisomy[None, :]
    tp = np.sum(any_pos & is_trisomy_b, axis=1).astype(np.int64)
    fn = np.sum((~any_pos) & is_trisomy_b, axis=1).astype(np.int64)
    fp = np.sum(any_pos & is_normal_b, axis=1).astype(np.int64)
    tn = np.sum((~any_pos) & is_normal_b, axis=1).astype(np.int64)
    return tp, tn, fp, fn


def _evaluate_batch_dual(
    chosen_batch: np.ndarray,
    pool_hypo_z: np.ndarray,
    pool_hyper_z: np.ndarray,
    dev_hypo_z: np.ndarray,
    dev_hyper_z: np.ndarray,
    dev_w_hypo: np.ndarray,
    dev_w_hyper: np.ndarray,
    dev_total_w: np.ndarray,
    dev_is_normal: np.ndarray,
    dev_is_trisomy: np.ndarray,
    test_hypo_z: np.ndarray,
    test_hyper_z: np.ndarray,
    test_w_hypo: np.ndarray,
    test_w_hyper: np.ndarray,
    test_total_w: np.ndarray,
    test_is_normal: np.ndarray,
    test_is_trisomy: np.ndarray,
    s_inter_cutoff: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized dev + test MCC for B reference draws.

    Returns
        mcc_dev, mcc_test,
        tp_dev, tn_dev, fp_dev, fn_dev,
        tp_test, tn_test, fp_test, fn_test
    each shaped (B,).
    """
    ref_hypo = pool_hypo_z[chosen_batch]
    ref_hyper = pool_hyper_z[chosen_batch]

    with np.errstate(invalid="ignore"):
        hypo_means = np.nanmean(ref_hypo, axis=1)
        hypo_stds = np.nanstd(ref_hypo, axis=1, ddof=0)
        hyper_means = np.nanmean(ref_hyper, axis=1)
        hyper_stds = np.nanstd(ref_hyper, axis=1, ddof=0)
    del ref_hypo, ref_hyper

    hypo_means = np.where(np.isfinite(hypo_means), hypo_means, 0.0)
    hyper_means = np.where(np.isfinite(hyper_means), hyper_means, 0.0)
    hypo_std_safe = np.where(hypo_stds > 0, hypo_stds, np.nan)
    hyper_std_safe = np.where(hyper_stds > 0, hyper_stds, np.nan)

    hypo_means_b = hypo_means[:, None, :]
    hyper_means_b = hyper_means[:, None, :]
    hypo_std_b = hypo_std_safe[:, None, :]
    hyper_std_b = hyper_std_safe[:, None, :]

    def _s_inter_max(
        eval_hypo_z: np.ndarray,
        eval_hyper_z: np.ndarray,
        eval_w_hypo: np.ndarray,
        eval_w_hyper: np.ndarray,
        eval_total_w: np.ndarray,
    ) -> np.ndarray:
        eval_hypo_z_b = eval_hypo_z[None, :, :]
        eval_hyper_z_b = eval_hyper_z[None, :, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            hypo_z_inter = (eval_hypo_z_b - hypo_means_b) / hypo_std_b
            hyper_z_inter = (eval_hyper_z_b - hyper_means_b) / hyper_std_b
            s_inter = (
                hyper_z_inter * eval_w_hyper[None, :, :]
                - hypo_z_inter * eval_w_hypo[None, :, :]
            ) / eval_total_w[None, :, :]
        s_inter = np.where(np.isnan(s_inter), 0.0, s_inter)
        with np.errstate(invalid="ignore"):
            return np.max(s_inter, axis=2)

    dev_s_max = _s_inter_max(dev_hypo_z, dev_hyper_z, dev_w_hypo, dev_w_hyper, dev_total_w)
    test_s_max = _s_inter_max(
        test_hypo_z, test_hyper_z, test_w_hypo, test_w_hyper, test_total_w,
    )

    tp_dev, tn_dev, fp_dev, fn_dev = _confusion_for_eval(
        dev_s_max, dev_is_normal, dev_is_trisomy, s_inter_cutoff,
    )
    tp_test, tn_test, fp_test, fn_test = _confusion_for_eval(
        test_s_max, test_is_normal, test_is_trisomy, s_inter_cutoff,
    )
    mcc_dev = _mcc_from_counts(tp_dev, tn_dev, fp_dev, fn_dev)
    mcc_test = _mcc_from_counts(tp_test, tn_test, fp_test, fn_test)
    return (
        mcc_dev, mcc_test,
        tp_dev, tn_dev, fp_dev, fn_dev,
        tp_test, tn_test, fp_test, fn_test,
    )


class _Top100ReferenceTracker:
    """Keep top-100 reference lists per ref_n by sum(MCC_dev + MCC_test)."""

    def __init__(self, ref_n_values: List[int], limit: int = 100) -> None:
        self._limit = limit
        self._entries: Dict[int, List[Dict[str, object]]] = {n: [] for n in ref_n_values}

    def consider(
        self,
        ref_n: int,
        mcc_dev: float,
        mcc_test: float,
        reference_list: str,
    ) -> None:
        if ref_n not in self._entries:
            return
        score = mcc_dev + mcc_test
        if not np.isfinite(score):
            return
        entry = {
            "ref_n": int(ref_n),
            "MCC_dev": float(mcc_dev),
            "MCC_test": float(mcc_test),
            "reference_list": reference_list,
            "_score": float(score),
        }
        entries = self._entries[ref_n]
        if len(entries) < self._limit:
            entries.append(entry)
            if len(entries) == self._limit:
                entries.sort(key=lambda e: e["_score"], reverse=True)
            return
        if score <= entries[-1]["_score"]:
            return
        entries.append(entry)
        entries.sort(key=lambda e: e["_score"], reverse=True)
        del entries[self._limit:]

    def to_dataframe(self) -> pd.DataFrame:
        rows: List[Dict[str, object]] = []
        for ref_n in sorted(self._entries):
            for entry in self._entries[ref_n]:
                rows.append(
                    {
                        "ref_n": entry["ref_n"],
                        "MCC_dev": entry["MCC_dev"],
                        "MCC_test": entry["MCC_test"],
                        "reference_list": entry["reference_list"],
                    }
                )
        return pd.DataFrame(
            rows,
            columns=["ref_n", "MCC_dev", "MCC_test", "reference_list"],
        )


def _label_masks(labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    is_normal = (labels == "Normal")
    is_trisomy = np.array(
        [
            isinstance(l, str) and l.startswith("T") and l != "Normal"
            for l in labels
        ],
        dtype=bool,
    )
    return is_normal, is_trisomy


def _build_eval_weights(
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    w_hypo = np.sqrt(hypo_counts.astype(np.float64))
    w_hyper = np.sqrt(hyper_counts.astype(np.float64))
    total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)
    total_w = np.where(total_w > 0, total_w, np.nan)
    return w_hypo, w_hyper, total_w


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--combo-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Directory containing the source _analyze_zscore.tsv.gz and "
        "_reference_zscore.tsv.gz produced by beta_to_zscore.py for the "
        "fixed (threshold, recall) combo (default 0.5 / 0.65)."
    ),
)
@click.option(
    "--meta-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="Sample meta CSV providing sample, label, ref_type, set, ff_before_mq.",
)
@click.option(
    "--output-base",
    default=None,
    type=click.Path(file_okay=False, dir_okay=True),
    help=(
        "Directory under which the output subdir is created. Required unless "
        "--output-dir is set."
    ),
)
@click.option(
    "--output-subdir",
    default="enlarged_reference",
    type=str,
    help="Subdirectory name under --output-base (default: enlarged_reference).",
)
@click.option(
    "--output-dir",
    default=None,
    type=click.Path(file_okay=False, dir_okay=True),
    help=(
        "Final output directory (overrides --output-base/--output-subdir). "
        "Created if missing."
    ),
)
@click.option(
    "--mode",
    default="all",
    type=click.Choice(["all", "isolated"], case_sensitive=False),
    help=(
        "Workflow mode: 'all' (default, original combined eval MCC) or "
        "'isolated' (separate dev/test MCC with pool sampling)."
    ),
)
@click.option(
    "--chr-spec",
    default="1-22",
    type=str,
    help='Chromosome list to use (default "1-22").',
)
@click.option(
    "--n-min",
    default=10,
    type=int,
    help="Minimum reference size N (inclusive). Default 10.",
)
@click.option(
    "--n-max",
    default=None,
    type=int,
    help=(
        "Maximum reference size N (inclusive). Default: full size of the "
        "dev+Normal pool."
    ),
)
@click.option(
    "--n-step",
    default=1,
    type=int,
    help="Step between consecutive N values. Default 1.",
)
@click.option(
    "--runs",
    default=10000,
    type=int,
    help=(
        "Number of unique random combinations per N (default 10000). "
        "If C(pool_size, N) <= runs the full set of combinations is "
        "enumerated, so N == pool_size always yields exactly one run."
    ),
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Master RNG seed used to draw the unique combinations.",
)
@click.option(
    "--s-inter-cutoff",
    default=3.0,
    type=float,
    help="s_inter cutoff used for the TP/TN/FP/FN classification (default 3.0).",
)
@click.option(
    "--threads",
    default=None,
    type=int,
    help=(
        "Worker threads for batched run evaluation. Defaults to "
        "$SLURM_CPUS_PER_TASK if set, otherwise os.cpu_count()."
    ),
)
@click.option(
    "--batch-size",
    default=100,
    type=int,
    help=(
        "Number of runs per vectorized batch (default 100). Larger values "
        "amortize threading overhead at the cost of peak memory."
    ),
)
def main(
    combo_dir: str,
    meta_csv: str,
    output_base: Optional[str],
    output_subdir: str,
    output_dir: Optional[str],
    mode: str,
    chr_spec: str,
    n_min: int,
    n_max: Optional[int],
    n_step: int,
    runs: int,
    seed: int,
    s_inter_cutoff: float,
    threads: Optional[int],
    batch_size: int,
) -> None:
    """Enlarged-reference sweep (--mode all) or isolated dev/test MCC (--mode isolated)."""
    console.rule("[bold blue]Enlarged-reference Z-score recomputation")

    if output_dir:
        out_dir = Path(output_dir)
    elif output_base:
        out_dir = Path(output_base) / output_subdir
    else:
        raise click.ClickException(
            "Either --output-dir or --output-base must be provided."
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    combo_dir_path = Path(combo_dir)
    chr_list = _parse_chr_spec(chr_spec)
    mode = mode.lower()

    if threads is None or threads <= 0:
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "0") or 0) or (os.cpu_count() or 1)
    batch_size = max(1, int(batch_size))

    console.print("\n[bold]Input parameters[/bold]")
    console.print(f"  Mode            : {mode}")
    console.print(f"  Combo dir       : {combo_dir_path}")
    console.print(f"  Meta CSV        : {meta_csv}")
    console.print(f"  Output dir      : {out_dir}")
    console.print(f"  Chromosomes     : {len(chr_list)}  ({chr_list[0]}..{chr_list[-1]})")
    if mode == "all":
        console.print(f"  N range         : [{n_min}, {n_max if n_max else 'pool_size'}] step {n_step}")
        console.print(f"  Runs per N      : up to {runs} (capped by C(pool, N))")
    console.print(f"  Seed            : {seed}")
    console.print(f"  s_inter cutoff  : {s_inter_cutoff}")
    console.print(f"  Threads         : {threads}")
    console.print(f"  Batch size      : {batch_size}")

    try:
        if mode == "isolated":
            _run_isolated_mode(
                combo_dir_path=combo_dir_path,
                meta_csv=meta_csv,
                out_dir=out_dir,
                chr_list=chr_list,
                seed=seed,
                s_inter_cutoff=s_inter_cutoff,
                threads=threads,
                batch_size=batch_size,
            )
            console.rule("[bold green]Enlarged-reference recomputation complete")
            return

        _run_all_mode(
            combo_dir_path=combo_dir_path,
            meta_csv=meta_csv,
            out_dir=out_dir,
            chr_list=chr_list,
            n_min=n_min,
            n_max=n_max,
            n_step=n_step,
            runs=runs,
            seed=seed,
            s_inter_cutoff=s_inter_cutoff,
            threads=threads,
            batch_size=batch_size,
        )
        console.rule("[bold green]Enlarged-reference recomputation complete")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


def _run_all_mode(
    *,
    combo_dir_path: Path,
    meta_csv: str,
    out_dir: Path,
    chr_list: List[str],
    n_min: int,
    n_max: Optional[int],
    n_step: int,
    runs: int,
    seed: int,
    s_inter_cutoff: float,
    threads: int,
    batch_size: int,
) -> None:
    """Original workflow: combined eval set, single MCC per run."""
    # --------------------------------------------------------------
    # Step 1: load combo TSVs (only needed columns)
    # --------------------------------------------------------------
    console.print("\n[bold cyan]Step 1: Loading combo TSVs[/bold cyan]")
    samples, arrays = _load_combined(combo_dir_path, chr_list)
    n_samples = len(samples)
    console.print(
        f"[green]OK[/green] Loaded {n_samples} samples x {len(chr_list)} chrs "
        f"(reference + analyze concatenated, kept only z_intra & CpG counts)"
    )

    # --------------------------------------------------------------
    # Step 2: meta join + build pool / fixed eval set
    # --------------------------------------------------------------
    console.print("\n[bold cyan]Step 2: Joining meta and partitioning samples[/bold cyan]")
    meta = _load_meta(meta_csv, samples)
    labels = meta["label"].fillna("").astype(str).to_numpy()
    sets = meta["set"].fillna("").astype(str).to_numpy()

    pool_mask = (sets == "dev") & (labels == "Normal")
    pool_indices = np.flatnonzero(pool_mask)
    pool_size = int(pool_indices.size)
    console.print(
        f"[green]OK[/green] Reference candidate pool (set=dev, label=Normal): "
        f"{pool_size} samples"
    )
    if pool_size < n_min:
        raise ValueError(
            f"Reference pool ({pool_size}) is smaller than --n-min ({n_min})."
        )

    eval_mask = ~pool_mask
    eval_indices = np.flatnonzero(eval_mask)
    eval_labels = labels[eval_indices]
    eval_is_normal = (eval_labels == "Normal")
    eval_is_trisomy = np.array(
        [
            isinstance(l, str) and l.startswith("T") and l != "Normal"
            for l in eval_labels
        ],
        dtype=bool,
    )
    n_eval = int(eval_indices.size)
    n_eval_pos = int(eval_is_trisomy.sum())
    n_eval_neg = int(eval_is_normal.sum())
    console.print(
        f"[green]OK[/green] Fixed eval set (test set + dev non-Normal): "
        f"{n_eval} samples ({n_eval_pos} positive / {n_eval_neg} negative, "
        f"{n_eval - n_eval_pos - n_eval_neg} ignored)"
    )

    if n_max is None:
        n_max = pool_size
    n_max = min(int(n_max), pool_size)
    if n_max < n_min:
        raise ValueError(f"--n-max ({n_max}) is smaller than --n-min ({n_min}).")

    n_values = list(range(n_min, n_max + 1, max(1, n_step)))
    if n_values and n_values[-1] != n_max:
        n_values.append(n_max)
    console.print(f"  N values        : {n_values[0]}..{n_values[-1]} ({len(n_values)} values)")

    # --------------------------------------------------------------
    # Pre-slice arrays into pool and eval views (much smaller)
    # --------------------------------------------------------------
    pool_hypo_z = np.ascontiguousarray(arrays["hypo_z_intra"][pool_indices])
    pool_hyper_z = np.ascontiguousarray(arrays["hyper_z_intra"][pool_indices])
    eval_hypo_z = np.ascontiguousarray(arrays["hypo_z_intra"][eval_indices])
    eval_hyper_z = np.ascontiguousarray(arrays["hyper_z_intra"][eval_indices])
    eval_hypo_counts = arrays["hypo_counts"][eval_indices].astype(np.float64)
    eval_hyper_counts = arrays["hyper_counts"][eval_indices].astype(np.float64)
    # Drop the full-sample arrays; everything from here on uses the sliced views.
    del arrays

    eval_w_hypo = np.sqrt(eval_hypo_counts)
    eval_w_hyper = np.sqrt(eval_hyper_counts)
    eval_total_w = np.sqrt(eval_w_hypo ** 2 + eval_w_hyper ** 2)
    eval_total_w = np.where(eval_total_w > 0, eval_total_w, np.nan)
    del eval_hypo_counts, eval_hyper_counts

    # --------------------------------------------------------------
    # Step 3: plan runs and dispatch in vectorized batches across threads
    # --------------------------------------------------------------
    rng = np.random.default_rng(seed)

    # Plan the per-N run counts so we can show one global progress bar.
    per_n_run_count: Dict[int, int] = {}
    for N in n_values:
        total_combos = math.comb(int(pool_size), int(N))
        per_n_run_count[N] = int(min(runs, total_combos))
    total_runs = sum(per_n_run_count.values())

    console.print(
        f"\n[bold cyan]Step 3: Recomputing s_inter across "
        f"{len(n_values)} N values, {total_runs} total runs "
        f"(threads={threads}, batch={batch_size})[/bold cyan]"
    )

    report_rows: List[Dict[str, float]] = []

    progress_columns = (
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    )

    with Progress(*progress_columns, console=console) as progress:
        task = progress.add_task("Runs", total=total_runs)

        for N in n_values:
            combos = _generate_unique_combos(pool_size, N, runs, rng)
            if not combos:
                continue

            # Stack combos for this N into a contiguous (R, N) array.
            chosen_full = np.stack(combos, axis=0)
            n_runs_this_N = chosen_full.shape[0]

            # Build batch slices.
            batches: List[Tuple[int, int]] = []
            for start in range(0, n_runs_this_N, batch_size):
                end = min(start + batch_size, n_runs_this_N)
                batches.append((start, end))

            with ThreadPoolExecutor(max_workers=threads) as ex:
                fut_to_range = {}
                for (start, end) in batches:
                    fut = ex.submit(
                        _evaluate_batch,
                        chosen_full[start:end],
                        pool_hypo_z,
                        pool_hyper_z,
                        eval_hypo_z,
                        eval_hyper_z,
                        eval_w_hypo,
                        eval_w_hyper,
                        eval_total_w,
                        eval_is_normal,
                        eval_is_trisomy,
                        s_inter_cutoff,
                    )
                    fut_to_range[fut] = (start, end)

                for fut in as_completed(fut_to_range):
                    start, end = fut_to_range[fut]
                    tp, tn, fp, fn, mcc = fut.result()
                    for i, run_index in enumerate(range(start, end)):
                        report_rows.append(
                            {
                                "ref_n": int(N),
                                "run_index": int(run_index),
                                "MCC": float(mcc[i]),
                                "TP": int(tp[i]),
                                "TN": int(tn[i]),
                                "FP": int(fp[i]),
                                "FN": int(fn[i]),
                            }
                        )
                    progress.advance(task, advance=(end - start))

            # Release combos for this N before moving on.
            del combos, chosen_full

    # --------------------------------------------------------------
    # Step 4: write report
    # --------------------------------------------------------------
    report_df = pd.DataFrame(
        report_rows,
        columns=["ref_n", "run_index", "MCC", "TP", "TN", "FP", "FN"],
    )
    report_df = report_df.sort_values(["ref_n", "run_index"]).reset_index(drop=True)
    report_path = out_dir / "report.csv"
    report_df.to_csv(report_path, index=False, float_format="%.6f")
    console.print(
    f"\n[green]OK[/green] Wrote report: {report_path} "
    f"({len(report_df)} rows)"
    )


# Isolated-mode fixed parameters (see module docstring).
_ISOLATED_POOL_SIZE = 50
_ISOLATED_POOL_RUNS = 100
_ISOLATED_REF_RUNS = 50
_ISOLATED_N_MIN = 5
_ISOLATED_N_MAX = 45
_ISOLATED_TOP_N_MIN = 35
_ISOLATED_TOP_N_MAX = 45
_ISOLATED_TOP_LIST_SIZE = 100


def _run_isolated_mode(
    *,
    combo_dir_path: Path,
    meta_csv: str,
    out_dir: Path,
    chr_list: List[str],
    seed: int,
    s_inter_cutoff: float,
    threads: int,
    batch_size: int,
) -> None:
    """Isolated dev/test MCC with random reference pools and top-100 tracking."""
    console.print("\n[bold cyan]Isolated mode parameters[/bold cyan]")
    console.print(f"  Pool size       : {_ISOLATED_POOL_SIZE}")
    console.print(f"  Pool runs       : {_ISOLATED_POOL_RUNS}")
    console.print(f"  Ref runs per N  : {_ISOLATED_REF_RUNS}")
    console.print(
        f"  N range         : [{_ISOLATED_N_MIN}, {_ISOLATED_N_MAX}]"
    )
    console.print(
        f"  Top-100 N range : [{_ISOLATED_TOP_N_MIN}, {_ISOLATED_TOP_N_MAX}]"
    )

    console.print("\n[bold cyan]Step 1: Loading combo TSVs[/bold cyan]")
    samples, arrays = _load_combined(combo_dir_path, chr_list)
    sample_arr = np.asarray(samples, dtype=object)
    console.print(
        f"[green]OK[/green] Loaded {len(samples)} samples x {len(chr_list)} chrs"
    )

    console.print("\n[bold cyan]Step 2: Joining meta and partitioning samples[/bold cyan]")
    meta = _load_meta(meta_csv, samples)
    labels = meta["label"].fillna("").astype(str).to_numpy()
    sets = meta["set"].fillna("").astype(str).to_numpy()

    dev_normal_mask = (sets == "dev") & (labels == "Normal")
    dev_normal_indices = np.flatnonzero(dev_normal_mask)
    dev_normal_size = int(dev_normal_indices.size)
    if dev_normal_size < _ISOLATED_POOL_SIZE:
        raise ValueError(
            f"Dev Normal pool ({dev_normal_size}) is smaller than required "
            f"pool size ({_ISOLATED_POOL_SIZE})."
        )

    all_dev_mask = sets == "dev"
    all_dev_indices = np.flatnonzero(all_dev_mask)
    test_indices = np.flatnonzero(sets == "test")

    test_labels = labels[test_indices]
    test_is_normal, test_is_trisomy = _label_masks(test_labels)
    console.print(
        f"[green]OK[/green] Dev Normal candidates: {dev_normal_size}; "
        f"all dev: {all_dev_indices.size}; test: {test_indices.size}"
    )

    n_values = list(range(_ISOLATED_N_MIN, _ISOLATED_N_MAX + 1))
    top_n_values = list(range(_ISOLATED_TOP_N_MIN, _ISOLATED_TOP_N_MAX + 1))
    top_tracker = _Top100ReferenceTracker(top_n_values, _ISOLATED_TOP_LIST_SIZE)

    rng = np.random.default_rng(seed)
    pool_combos = _generate_unique_combos(
        dev_normal_size,
        _ISOLATED_POOL_SIZE,
        _ISOLATED_POOL_RUNS,
        rng,
    )

    total_runs = _ISOLATED_POOL_RUNS * len(n_values) * _ISOLATED_REF_RUNS
    console.print(
        f"\n[bold cyan]Step 3: Isolated sweep "
        f"({total_runs} planned runs, threads={threads})[/bold cyan]"
    )

    report_rows: List[Dict[str, object]] = []
    progress_columns = (
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    )

    with Progress(*progress_columns, console=console) as progress:
        task = progress.add_task("Runs", total=total_runs)

        for pool_run_index, pool_local in enumerate(pool_combos):
            pool_global = dev_normal_indices[pool_local]
            pool_global_set = set(pool_global.tolist())

            dev_eval_indices = np.array(
                [i for i in all_dev_indices if i not in pool_global_set],
                dtype=np.int64,
            )
            dev_eval_labels = labels[dev_eval_indices]
            dev_is_normal, dev_is_trisomy = _label_masks(dev_eval_labels)

            pool_hypo_z = np.ascontiguousarray(arrays["hypo_z_intra"][pool_global])
            pool_hyper_z = np.ascontiguousarray(arrays["hyper_z_intra"][pool_global])

            dev_hypo_z = np.ascontiguousarray(arrays["hypo_z_intra"][dev_eval_indices])
            dev_hyper_z = np.ascontiguousarray(arrays["hyper_z_intra"][dev_eval_indices])
            dev_w_hypo, dev_w_hyper, dev_total_w = _build_eval_weights(
                arrays["hypo_counts"][dev_eval_indices],
                arrays["hyper_counts"][dev_eval_indices],
            )

            test_hypo_z = np.ascontiguousarray(arrays["hypo_z_intra"][test_indices])
            test_hyper_z = np.ascontiguousarray(arrays["hyper_z_intra"][test_indices])
            test_w_hypo, test_w_hyper, test_total_w = _build_eval_weights(
                arrays["hypo_counts"][test_indices],
                arrays["hyper_counts"][test_indices],
            )

            for N in n_values:
                ref_combos = _generate_unique_combos(
                    _ISOLATED_POOL_SIZE,
                    N,
                    _ISOLATED_REF_RUNS,
                    rng,
                )
                chosen_full = np.stack(ref_combos, axis=0)
                n_runs = chosen_full.shape[0]
                batches: List[Tuple[int, int]] = []
                for start in range(0, n_runs, batch_size):
                    end = min(start + batch_size, n_runs)
                    batches.append((start, end))

                with ThreadPoolExecutor(max_workers=threads) as ex:
                    fut_to_range: Dict = {}
                    for start, end in batches:
                        fut = ex.submit(
                            _evaluate_batch_dual,
                            chosen_full[start:end],
                            pool_hypo_z,
                            pool_hyper_z,
                            dev_hypo_z,
                            dev_hyper_z,
                            dev_w_hypo,
                            dev_w_hyper,
                            dev_total_w,
                            dev_is_normal,
                            dev_is_trisomy,
                            test_hypo_z,
                            test_hyper_z,
                            test_w_hypo,
                            test_w_hyper,
                            test_total_w,
                            test_is_normal,
                            test_is_trisomy,
                            s_inter_cutoff,
                        )
                        fut_to_range[fut] = (start, end)

                    for fut in as_completed(fut_to_range):
                        start, end = fut_to_range[fut]
                        (
                            mcc_dev,
                            mcc_test,
                            tp_dev,
                            tn_dev,
                            fp_dev,
                            fn_dev,
                            tp_test,
                            tn_test,
                            fp_test,
                            fn_test,
                        ) = fut.result()
                        for i, ref_run_index in enumerate(range(start, end)):
                            ref_local = chosen_full[ref_run_index]
                            ref_global = pool_global[ref_local]
                            ref_names = sample_arr[ref_global]
                            ref_list_str = ",".join(ref_names.tolist())
                            mcc_d = float(mcc_dev[i])
                            mcc_t = float(mcc_test[i])
                            report_rows.append(
                                {
                                    "pool_run_index": int(pool_run_index),
                                    "ref_run_index": int(ref_run_index),
                                    "ref_n": int(N),
                                    "MCC_dev": mcc_d,
                                    "MCC_test": mcc_t,
                                    "TP_dev": int(tp_dev[i]),
                                    "TN_dev": int(tn_dev[i]),
                                    "FP_dev": int(fp_dev[i]),
                                    "FN_dev": int(fn_dev[i]),
                                    "TP_test": int(tp_test[i]),
                                    "TN_test": int(tn_test[i]),
                                    "FP_test": int(fp_test[i]),
                                    "FN_test": int(fn_test[i]),
                                }
                            )
                            top_tracker.consider(N, mcc_d, mcc_t, ref_list_str)
                        progress.advance(task, advance=(end - start))

                del ref_combos, chosen_full

    report_df = pd.DataFrame(
        report_rows,
        columns=[
            "pool_run_index",
            "ref_run_index",
            "ref_n",
            "MCC_dev",
            "MCC_test",
            "TP_dev",
            "TN_dev",
            "FP_dev",
            "FN_dev",
            "TP_test",
            "TN_test",
            "FP_test",
            "FN_test",
        ],
    )
    report_df = report_df.sort_values(
        ["pool_run_index", "ref_n", "ref_run_index"],
    ).reset_index(drop=True)
    report_path = out_dir / "report.csv"
    report_df.to_csv(report_path, index=False, float_format="%.6f")
    console.print(
        f"\n[green]OK[/green] Wrote report: {report_path} ({len(report_df)} rows)"
    )

    top_df = top_tracker.to_dataframe()
    top_df = top_df.sort_values(
        ["ref_n", "MCC_dev", "MCC_test"],
        ascending=[True, False, False],
    ).reset_index(drop=True)
    top_path = out_dir / "top_100_reference_list.csv"
    top_df.to_csv(top_path, index=False, float_format="%.6f")
    console.print(
        f"[green]OK[/green] Wrote top-100 list: {top_path} ({len(top_df)} rows)"
    )


# ---------------------------------------------------------------------------
# Small parser for chromosome specs (mirrors beta_to_zscore.py's behaviour)
# ---------------------------------------------------------------------------

def _parse_chr_spec(spec: str) -> List[str]:
    """Parse '1-22' or '1,2,X' or 'chr1,chr2' into ['chr1', 'chr2', ...]."""
    spec = spec.strip()
    tokens = [s.strip() for s in spec.split(",") if s.strip()]
    out: List[str] = []
    for token in tokens:
        if "-" in token and not token.startswith("chr"):
            try:
                start, end = token.split("-")
                start_i = int(start)
                end_i = int(end)
                out.extend([f"chr{i}" for i in range(start_i, end_i + 1)])
                continue
            except Exception:
                pass
        out.append(token if token.startswith("chr") else f"chr{token}")
    return out


if __name__ == "__main__":
    main()
