#!/usr/bin/env python3
"""
Enlarged-reference per-chromosome Z-score recomputation.

Goal:
    Test how the per-chromosome ``s_inter`` values respond to growing the
    reference set with additional dev-set Normal samples. For each enlarged
    reference size ``N`` (and several random draws per ``N``) we recompute
    ``hypo_z_inter`` / ``hyper_z_inter`` / ``s_inter`` for every non-reference
    sample using a fixed combo (threshold=0.5, recall=0.65), then report a
    per-run MCC.

Inputs:
    --combo-dir : directory containing ``_analyze_zscore.tsv.gz`` and
                  ``_reference_zscore.tsv.gz`` produced by ``beta_to_zscore.py``
                  for the fixed (threshold, recall) combo.
    --meta-csv  : sample-level metadata. Must provide
                  ``sample, label, ref_type, set, ff_before_mq``.

The reference candidate pool is ``set == "dev"`` AND ``label == "Normal"``.
For each ``N`` from ``--n-min`` (default 10) to that pool's size, ``--runs``
random draws (default 10) are generated. ``N == pool_size`` always uses a
single run (only one selection is possible).

Outputs (under ``--output-base/<output-subdir>/``, default
``enlarged_reference``):

    ref_n_{N}_run_{run_index}/_reference_zscore.tsv.gz
    ref_n_{N}_run_{run_index}/_analyze_zscore.tsv.gz
    report.csv

``report.csv`` columns: ``ref_n, run_index, MCC, TP, TN, FP, FN``.

The classification rule for the MCC counts is:
    Positive : ``label`` starts with ``T`` (i.e. trisomy)
    Negative : ``label == "Normal"``
    Predicted positive : ``any(chr_s_inter > --s-inter-cutoff)``
    Predicted negative : otherwise
"""

from __future__ import annotations

import math
import os
import sys
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

# Reuse the writer (and the column ordering it implies) from beta_to_zscore.
_THIS_DIR = Path(__file__).resolve().parent
_GRID_SEARCH_DIR = _THIS_DIR.parent / "grid_search"
if str(_GRID_SEARCH_DIR) not in sys.path:
    sys.path.insert(0, str(_GRID_SEARCH_DIR))
from beta_to_zscore import _write_tsv  # noqa: E402

console = Console()


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

_REF_SUFFIXES: Tuple[str, ...] = (
    "hypo_beta",
    "hyper_beta",
    "hypo_z_intra",
    "hyper_z_intra",
    "s_intra",
    "hypo_cpgs_count",
    "hyper_cpgs_count",
)


def _load_combo_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing combo output: {path}")
    return pd.read_csv(path, sep="\t", compression="gzip")


def _stack_per_chr(
    df: pd.DataFrame,
    chr_list: List[str],
    suffix: str,
    *,
    fill: float = np.nan,
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
        if not math.isnan(fill):
            arr = np.where(np.isnan(arr), fill, arr)
    return arr


def _load_combined(
    combo_dir: Path,
    chr_list: List[str],
) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Concat reference + analyze TSVs and return per-metric matrices.

    Returns
    -------
    samples : list[str]
        Sample IDs (reference rows first, then analyze rows; the order is
        otherwise preserved from the source files).
    arrays : dict
        Keys: ``hypo_beta``, ``hyper_beta``, ``hypo_z_intra``, ``hyper_z_intra``,
        ``s_intra``, ``hypo_counts``, ``hyper_counts``. Each value is a
        ``(n_samples, n_chr)`` numpy matrix.
    """
    ref_df = _load_combo_table(combo_dir / "_reference_zscore.tsv.gz")
    an_df = _load_combo_table(combo_dir / "_analyze_zscore.tsv.gz")

    keep = ["sample"] + [
        f"{c}_{suf}" for c in chr_list for suf in _REF_SUFFIXES
    ]
    ref_df = ref_df.loc[:, keep].copy()
    an_df = an_df.loc[:, keep].copy()

    combined = pd.concat([ref_df, an_df], axis=0, ignore_index=True)
    if combined["sample"].duplicated().any():
        dups = combined["sample"][combined["sample"].duplicated()].tolist()
        raise ValueError(
            f"Duplicate sample IDs across reference + analyze: {dups[:5]} "
            f"(total {len(dups)})"
        )

    samples = combined["sample"].astype(str).tolist()

    arrays = {
        "hypo_beta": _stack_per_chr(combined, chr_list, "hypo_beta"),
        "hyper_beta": _stack_per_chr(combined, chr_list, "hyper_beta"),
        "hypo_z_intra": _stack_per_chr(combined, chr_list, "hypo_z_intra"),
        "hyper_z_intra": _stack_per_chr(combined, chr_list, "hyper_z_intra"),
        "s_intra": _stack_per_chr(combined, chr_list, "s_intra"),
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
# Statistics + s_inter (vectorized over samples and chrs at once)
# ---------------------------------------------------------------------------

def _reference_stats(
    hypo_z: np.ndarray,
    hyper_z: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-chr (mean, std) of hypo and hyper z_intra, NaN-safe."""
    with np.errstate(invalid="ignore"):
        hypo_means = np.nanmean(hypo_z, axis=0)
        hypo_stds = np.nanstd(hypo_z, axis=0, ddof=0)
        hyper_means = np.nanmean(hyper_z, axis=0)
        hyper_stds = np.nanstd(hyper_z, axis=0, ddof=0)
    hypo_means = np.where(np.isfinite(hypo_means), hypo_means, 0.0)
    hyper_means = np.where(np.isfinite(hyper_means), hyper_means, 0.0)
    return hypo_means, hypo_stds, hyper_means, hyper_stds


def _s_inter_block(
    hypo_z_intra: np.ndarray,
    hyper_z_intra: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
    hypo_means: np.ndarray,
    hypo_stds: np.ndarray,
    hyper_means: np.ndarray,
    hyper_stds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized version of ``calculate_s_inter_from_stats`` over all samples.

    Inputs are ``(n_samples, n_chr)``. Means/stds are ``(n_chr,)``.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        hypo_std_safe = np.where(hypo_stds > 0, hypo_stds, np.nan)
        hyper_std_safe = np.where(hyper_stds > 0, hyper_stds, np.nan)
        hypo_z_inter = (hypo_z_intra - hypo_means) / hypo_std_safe
        hyper_z_inter = (hyper_z_intra - hyper_means) / hyper_std_safe

    w_hypo = np.sqrt(hypo_counts.astype(np.float64))
    w_hyper = np.sqrt(hyper_counts.astype(np.float64))
    total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        s_inter = (hyper_z_inter * w_hyper - hypo_z_inter * w_hypo) / total_w
    s_inter = np.where(np.isnan(s_inter), 0.0, s_inter)
    return hypo_z_inter, hyper_z_inter, s_inter


# ---------------------------------------------------------------------------
# Output construction
# ---------------------------------------------------------------------------

def _build_output_df(
    samples: List[str],
    chr_list: List[str],
    hypo_beta: np.ndarray,
    hyper_beta: np.ndarray,
    hypo_z_intra: np.ndarray,
    hyper_z_intra: np.ndarray,
    s_intra: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
    inter: Optional[Dict[str, np.ndarray]] = None,
) -> pd.DataFrame:
    """Build a wide DataFrame matching ``beta_to_zscore.py``'s schema."""
    data: Dict[str, np.ndarray] = {"sample": np.asarray(samples, dtype=object)}
    for j, chr_name in enumerate(chr_list):
        data[f"{chr_name}_hypo_beta"] = hypo_beta[:, j]
        data[f"{chr_name}_hyper_beta"] = hyper_beta[:, j]
        data[f"{chr_name}_hypo_z_intra"] = hypo_z_intra[:, j]
        data[f"{chr_name}_hyper_z_intra"] = hyper_z_intra[:, j]
        data[f"{chr_name}_s_intra"] = s_intra[:, j]
        data[f"{chr_name}_hypo_cpgs_count"] = hypo_counts[:, j]
        data[f"{chr_name}_hyper_cpgs_count"] = hyper_counts[:, j]
        if inter is not None:
            data[f"{chr_name}_hypo_z_inter"] = inter["hypo_z_inter"][:, j]
            data[f"{chr_name}_hyper_z_inter"] = inter["hyper_z_inter"][:, j]
            data[f"{chr_name}_s_inter"] = inter["s_inter"][:, j]
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# MCC
# ---------------------------------------------------------------------------

def _classify(
    s_inter: np.ndarray,
    labels: np.ndarray,
    cutoff: float,
) -> Tuple[int, int, int, int]:
    """Return ``(TP, TN, FP, FN)`` over the rows of ``s_inter``.

    ``labels`` is a 1-D string array aligned to ``s_inter``'s rows.
    Rows whose label is missing or not in {"Normal", "T*"} are skipped.
    """
    with np.errstate(invalid="ignore"):
        any_pos = np.nanmax(s_inter, axis=1) > cutoff
    any_pos = np.where(np.isnan(any_pos), False, any_pos).astype(bool)

    is_normal = labels == "Normal"
    is_trisomy = np.array(
        [isinstance(l, str) and l.startswith("T") and l != "Normal" for l in labels],
        dtype=bool,
    )

    tp = int(np.sum(is_trisomy & any_pos))
    fn = int(np.sum(is_trisomy & ~any_pos))
    fp = int(np.sum(is_normal & any_pos))
    tn = int(np.sum(is_normal & ~any_pos))
    return tp, tn, fp, fn


def _mcc(tp: int, tn: int, fp: int, fn: int) -> float:
    num = tp * tn - fp * fn
    denom_sq = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if denom_sq == 0:
        return float("nan")
    return num / math.sqrt(denom_sq)


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
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Directory under which the output subdir is created.",
)
@click.option(
    "--output-subdir",
    default="enlarged_reference",
    type=str,
    help="Subdirectory name under --output-base (default: enlarged_reference).",
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
    default=10,
    type=int,
    help="Number of random draws per N (default 10). N == pool_size always uses 1.",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Base seed; per-run seed is ``seed * 100000 + N * 100 + run_index``.",
)
@click.option(
    "--s-inter-cutoff",
    default=3.0,
    type=float,
    help="s_inter cutoff used for the TP/TN/FP/FN classification (default 3.0).",
)
@click.option(
    "--no-write-tsv",
    is_flag=True,
    default=False,
    help="Skip writing per-run _analyze/_reference TSVs (only emit report.csv).",
)
@click.option(
    "--no-gzip",
    is_flag=True,
    default=False,
    help="Write per-run outputs as plain TSV instead of TSV.gz.",
)
def main(
    combo_dir: str,
    meta_csv: str,
    output_base: str,
    output_subdir: str,
    chr_spec: str,
    n_min: int,
    n_max: Optional[int],
    n_step: int,
    runs: int,
    seed: int,
    s_inter_cutoff: float,
    no_write_tsv: bool,
    no_gzip: bool,
) -> None:
    """Sweep over reference sizes ``N`` and report MCC for each random draw.

    \b
    1. Concat ``_reference_zscore.tsv.gz`` + ``_analyze_zscore.tsv.gz`` from
       --combo-dir (the fixed (threshold, recall) directory).
    2. Join with --meta-csv to attach label / ref_type / set.
    3. For ``N`` in [n_min, n_min+n_step, ..., n_max] and
       ``run_index`` in [0, runs):
         - randomly pick N samples from {set==dev, label==Normal} as the new
           reference;
         - recompute reference per-chr stats and analyze ``s_inter``;
         - write ref_n_{N}_run_{run_index}/_reference_zscore.tsv.gz and
           _analyze_zscore.tsv.gz (unless --no-write-tsv);
         - record TP/TN/FP/FN/MCC over all non-reference samples.
    4. Write report.csv under --output-base/<output-subdir>/.
    """
    console.rule("[bold blue]Enlarged-reference Z-score recomputation")

    output_base_path = Path(output_base)
    out_dir = output_base_path / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    combo_dir_path = Path(combo_dir)
    chr_list = _parse_chr_spec(chr_spec)

    console.print("\n[bold]Input parameters[/bold]")
    console.print(f"  Combo dir       : {combo_dir_path}")
    console.print(f"  Meta CSV        : {meta_csv}")
    console.print(f"  Output dir      : {out_dir}")
    console.print(f"  Chromosomes     : {len(chr_list)}  ({chr_list[0]}..{chr_list[-1]})")
    console.print(f"  N range         : [{n_min}, {n_max if n_max else 'pool_size'}] step {n_step}")
    console.print(f"  Runs per N      : {runs}")
    console.print(f"  Seed (base)     : {seed}")
    console.print(f"  s_inter cutoff  : {s_inter_cutoff}")
    console.print(f"  Write per-run TSVs: {not no_write_tsv}")

    try:
        # --------------------------------------------------------------
        # Step 1: load combo TSVs
        # --------------------------------------------------------------
        console.print("\n[bold cyan]Step 1: Loading combo TSVs[/bold cyan]")
        samples, arrays = _load_combined(combo_dir_path, chr_list)
        n_samples = len(samples)
        console.print(
            f"[green]OK[/green] Loaded {n_samples} samples x {len(chr_list)} chrs "
            f"(reference + analyze concatenated)"
        )

        # --------------------------------------------------------------
        # Step 2: meta + reference candidate pool
        # --------------------------------------------------------------
        console.print("\n[bold cyan]Step 2: Joining meta and building ref pool[/bold cyan]")
        meta = _load_meta(meta_csv, samples)
        labels = meta["label"].astype("string").to_numpy()
        sets = meta["set"].astype("string").to_numpy()

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

        if n_max is None:
            n_max = pool_size
        n_max = min(int(n_max), pool_size)
        if n_max < n_min:
            raise ValueError(f"--n-max ({n_max}) is smaller than --n-min ({n_min}).")

        n_values = list(range(n_min, n_max + 1, max(1, n_step)))
        if n_values[-1] != n_max:
            n_values.append(n_max)
        console.print(f"  N values        : {n_values}")

        # Pre-extract arrays we'll reuse hot.
        hypo_z_intra = arrays["hypo_z_intra"]
        hyper_z_intra = arrays["hyper_z_intra"]
        hypo_counts = arrays["hypo_counts"]
        hyper_counts = arrays["hyper_counts"]
        hypo_beta = arrays["hypo_beta"]
        hyper_beta = arrays["hyper_beta"]
        s_intra = arrays["s_intra"]

        suffix = ".tsv" if no_gzip else ".tsv.gz"

        # --------------------------------------------------------------
        # Step 3: sweep over N and runs
        # --------------------------------------------------------------
        total_runs = sum(1 if N == pool_size else runs for N in n_values)
        console.print(
            f"\n[bold cyan]Step 3: Recomputing s_inter across "
            f"{len(n_values)} N values, {total_runs} total runs[/bold cyan]"
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
                this_runs = 1 if N == pool_size else runs
                for run_index in range(this_runs):
                    rng = np.random.default_rng(seed * 100000 + N * 100 + run_index)
                    chosen = np.sort(rng.choice(pool_indices, size=N, replace=False))

                    is_ref = np.zeros(n_samples, dtype=bool)
                    is_ref[chosen] = True
                    analyze_idx = np.flatnonzero(~is_ref)

                    hypo_means, hypo_stds, hyper_means, hyper_stds = _reference_stats(
                        hypo_z_intra[chosen], hyper_z_intra[chosen],
                    )

                    a_hypo_z_inter, a_hyper_z_inter, a_s_inter = _s_inter_block(
                        hypo_z_intra[analyze_idx], hyper_z_intra[analyze_idx],
                        hypo_counts[analyze_idx], hyper_counts[analyze_idx],
                        hypo_means, hypo_stds, hyper_means, hyper_stds,
                    )

                    if not no_write_tsv:
                        run_dir = out_dir / f"ref_n_{N}_run_{run_index}"
                        run_dir.mkdir(parents=True, exist_ok=True)

                        # Sort by sample for deterministic ordering.
                        ref_order = np.argsort([samples[i] for i in chosen])
                        an_order = np.argsort([samples[i] for i in analyze_idx])

                        ref_samples_sorted = [samples[chosen[i]] for i in ref_order]
                        an_samples_sorted = [samples[analyze_idx[i]] for i in an_order]

                        ref_df = _build_output_df(
                            ref_samples_sorted,
                            chr_list,
                            hypo_beta[chosen][ref_order],
                            hyper_beta[chosen][ref_order],
                            hypo_z_intra[chosen][ref_order],
                            hyper_z_intra[chosen][ref_order],
                            s_intra[chosen][ref_order],
                            hypo_counts[chosen][ref_order],
                            hyper_counts[chosen][ref_order],
                            inter=None,
                        )
                        an_df = _build_output_df(
                            an_samples_sorted,
                            chr_list,
                            hypo_beta[analyze_idx][an_order],
                            hyper_beta[analyze_idx][an_order],
                            hypo_z_intra[analyze_idx][an_order],
                            hyper_z_intra[analyze_idx][an_order],
                            s_intra[analyze_idx][an_order],
                            hypo_counts[analyze_idx][an_order],
                            hyper_counts[analyze_idx][an_order],
                            inter={
                                "hypo_z_inter": a_hypo_z_inter[an_order],
                                "hyper_z_inter": a_hyper_z_inter[an_order],
                                "s_inter": a_s_inter[an_order],
                            },
                        )
                        _write_tsv(ref_df, str(run_dir / f"_reference_zscore{suffix}"))
                        _write_tsv(an_df, str(run_dir / f"_analyze_zscore{suffix}"))

                    tp, tn, fp, fn = _classify(
                        a_s_inter, labels[analyze_idx], s_inter_cutoff,
                    )
                    mcc = _mcc(tp, tn, fp, fn)
                    report_rows.append(
                        {
                            "ref_n": N,
                            "run_index": run_index,
                            "MCC": mcc,
                            "TP": tp,
                            "TN": tn,
                            "FP": fp,
                            "FN": fn,
                        }
                    )
                    progress.advance(task)

        # --------------------------------------------------------------
        # Step 4: report
        # --------------------------------------------------------------
        report_df = pd.DataFrame(
            report_rows, columns=["ref_n", "run_index", "MCC", "TP", "TN", "FP", "FN"],
        )
        report_path = out_dir / "report.csv"
        report_df.to_csv(report_path, index=False, float_format="%.6f")
        console.print(
            f"\n[green]OK[/green] Wrote report: {report_path} "
            f"({len(report_df)} rows)"
        )

        console.rule("[bold green]Enlarged-reference recomputation complete")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


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
