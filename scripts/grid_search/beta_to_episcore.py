#!/usr/bin/env python3
"""
Beta to Z-score Conversion for Multi-sample Trisomy Detection.

This script processes methylation beta values for many analyze and reference
samples to calculate chromosome-level statistics and Z-scores for trisomy
detection.

Pipeline:
    1. Discover ``{sample}_beta_value.tsv.gz`` files in two directories
       (--reference-beta-value-dir and --analyze-beta-value-dir).
    2. For every sample (in parallel) read the beta file, optionally filter by
       depth and a CpG list, then compute per-chromosome aggregate beta values
       and intra-sample Z-scores (s_intra).
    3. Stack reference samples into a reference s_intra matrix and derive
       per-chromosome means and standard deviations of hypo/hyper z_intra.
    4. For every analyze sample, compute inter-sample Z-scores (s_inter)
       against the reference statistics.

Outputs (two TSVs, gzipped if the prefix path does not already end with .tsv):
    - ``{output_prefix}_reference_zscore.tsv.gz`` : reference s_intra matrix
    - ``{output_prefix}_analyze_zscore.tsv.gz``   : analyze s_intra + s_inter

Both outputs include a leading ``sample`` column.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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

# Module-level globals shared with worker processes via Unix fork (copy-on-write).
# They are populated by ``_init_worker_globals`` before the pool is forked, so each
# child inherits them without paying any pickling cost.
_CTX: Dict[str, Any] = {}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SAMPLE_FILE_RE = re.compile(r"^(?P<sample>.+)_beta_value\.tsv\.gz$")


def _parse_chr_list(spec: str) -> List[str]:
    """Parse a chromosome spec like '1-22', '1,2,X,Y', or mixed spec like '1-22,X'."""
    spec = spec.strip()
    tokens = [s.strip() for s in spec.split(",") if s.strip()]
    result = []
    for token in tokens:
        if "-" in token and not token.startswith("chr"):
            # Range token, e.g. "1-22"
            try:
                start, end = token.split("-")
                # Only numeric ranges are expanded as chr{n}
                start_int = int(start)
                end_int = int(end)
                result.extend([f"chr{i}" for i in range(start_int, end_int + 1)])
            except Exception:
                # If it isn't a numeric X-Y, treat as a literal chromosome string
                result.append(token if token.startswith("chr") else f"chr{token}")
        else:
            result.append(token if token.startswith("chr") else f"chr{token}")
    return result


def _discover_samples(directory: str) -> Dict[str, str]:
    """Return ``{sample: filepath}`` for every ``*_beta_value.tsv.gz`` in dir."""
    p = Path(directory)
    if not p.is_dir():
        raise FileNotFoundError(f"Not a directory: {directory}")
    out: Dict[str, str] = {}
    for f in sorted(p.iterdir()):
        m = _SAMPLE_FILE_RE.match(f.name)
        if m:
            out[m.group("sample")] = str(f)
    return out


# ---------------------------------------------------------------------------
# Per-sample numerical kernels
# ---------------------------------------------------------------------------

def read_beta_filtered(
    beta_path: str,
    usecols: List[str],
    chr_list: List[str],
    cpg_filter_df: Optional[pd.DataFrame] = None,
    filter_depth: Optional[int] = None,
    depth_col: Optional[str] = None,
) -> pd.DataFrame:
    """Read a gzipped beta-value TSV and apply chr/depth/CpG filters.

    The CpG filter is applied via an inner merge on ``(chr, start, end)`` which
    avoids the costly string-concatenation pattern used by the original
    single-sample script.
    """
    df = pd.read_csv(
        beta_path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
    )

    chr_set = set(chr_list)
    df = df[df["chr"].isin(chr_set)]

    if filter_depth is not None and depth_col is not None and depth_col in df.columns:
        df = df[df[depth_col] > filter_depth]

    if cpg_filter_df is not None:
        df = df.merge(cpg_filter_df, on=["chr", "start", "end"], how="inner", copy=False)

    return df


def _aggregate_chr(
    df: pd.DataFrame,
    chr_list: List[str],
    meth_col: str,
    unmeth_col: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(beta_per_chr, count_per_chr)`` aligned to ``chr_list``."""
    n = len(chr_list)
    if df.empty:
        return np.full(n, np.nan, dtype=np.float64), np.zeros(n, dtype=np.int64)

    grouped = df.groupby("chr", sort=False)
    meth = grouped[meth_col].sum().reindex(chr_list, fill_value=0).to_numpy(dtype=np.float64)
    unmeth = grouped[unmeth_col].sum().reindex(chr_list, fill_value=0).to_numpy(dtype=np.float64)
    counts = grouped.size().reindex(chr_list, fill_value=0).to_numpy(dtype=np.int64)

    denom = meth + unmeth
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(denom > 0, meth / denom, np.nan)
    return beta, counts


def calculate_chr_level_beta(
    df: pd.DataFrame,
    chr_list: List[str],
    meth_col: str,
    unmeth_col: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split into hypo/hyper by ``meandiff`` and compute chr-level betas+counts."""
    hypo = df[df["meandiff"] < 0]
    hyper = df[df["meandiff"] > 0]

    hypo_beta, hypo_counts = _aggregate_chr(hypo, chr_list, meth_col, unmeth_col)
    hyper_beta, hyper_counts = _aggregate_chr(hyper, chr_list, meth_col, unmeth_col)
    return hypo_beta, hyper_beta, hypo_counts, hyper_counts


def _zscore_across(arr: np.ndarray) -> np.ndarray:
    """Standardize a 1-D array (NaN-safe). Returns NaN-filled if std is 0/NaN."""
    if np.all(np.isnan(arr)):
        return np.full_like(arr, np.nan, dtype=np.float64)
    mean = np.nanmean(arr)
    std = np.nanstd(arr)
    if not np.isfinite(std) or std == 0:
        return np.full_like(arr, np.nan, dtype=np.float64)
    return (arr - mean) / std


def calculate_s_intra(
    hypo_beta: np.ndarray,
    hyper_beta: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute hypo/hyper z_intra and the weighted combined s_intra."""
    hypo_z = _zscore_across(hypo_beta)
    hyper_z = _zscore_across(hyper_beta)

    w_hypo = np.sqrt(hypo_counts.astype(np.float64))
    w_hyper = np.sqrt(hyper_counts.astype(np.float64))
    total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        s_intra = (hyper_z * w_hyper - hypo_z * w_hypo) / total_w
        s_intra = np.where(np.isnan(s_intra), 0.0, s_intra)
    return hypo_z, hyper_z, s_intra


def calculate_s_inter_from_stats(
    hypo_z_intra: np.ndarray,
    hyper_z_intra: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
    hypo_means: np.ndarray,
    hypo_stds: np.ndarray,
    hyper_means: np.ndarray,
    hyper_stds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Normalize z_intra against pre-computed reference stats and combine."""
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
# Worker plumbing
# ---------------------------------------------------------------------------

def _init_worker_globals(ctx: Dict[str, Any]) -> None:
    """Pool initializer. Stores ``ctx`` in the worker's ``_CTX`` global."""
    global _CTX
    _CTX = ctx


def _process_one_sample(item: Tuple[str, str]) -> Dict[str, Any]:
    """Worker entry point: returns per-sample arrays or an error record."""
    sample, file_path = item
    chr_list: List[str] = _CTX["chr_list"]
    beta_columns: List[str] = _CTX["beta_columns"]
    cpg_filter_df: Optional[pd.DataFrame] = _CTX.get("cpg_filter_df")
    depth: Optional[int] = _CTX.get("depth")
    depth_col: Optional[str] = _CTX.get("depth_col")

    try:
        df = read_beta_filtered(
            file_path,
            usecols=beta_columns,
            chr_list=chr_list,
            cpg_filter_df=cpg_filter_df,
            filter_depth=depth,
            depth_col=depth_col,
        )

        if "raw_meth_count" in df.columns and "raw_unmeth_count" in df.columns:
            meth_col, unmeth_col = "raw_meth_count", "raw_unmeth_count"
        elif "target_meth_count" in df.columns and "target_unmeth_count" in df.columns:
            meth_col, unmeth_col = "target_meth_count", "target_unmeth_count"
        else:
            raise ValueError(
                "Cannot find methylation count columns "
                "(expected raw_meth_count/raw_unmeth_count or "
                "target_meth_count/target_unmeth_count)."
            )

        hypo_beta, hyper_beta, hypo_counts, hyper_counts = calculate_chr_level_beta(
            df, chr_list, meth_col, unmeth_col,
        )
        # Free the dataframe ASAP to reduce peak memory of the worker.
        del df

        hypo_z, hyper_z, s_intra = calculate_s_intra(
            hypo_beta, hyper_beta, hypo_counts, hyper_counts,
        )

        return {
            "sample": sample,
            "hypo_beta": hypo_beta.astype(np.float64, copy=False),
            "hyper_beta": hyper_beta.astype(np.float64, copy=False),
            "hypo_z_intra": hypo_z.astype(np.float64, copy=False),
            "hyper_z_intra": hyper_z.astype(np.float64, copy=False),
            "s_intra": s_intra.astype(np.float64, copy=False),
            "hypo_counts": hypo_counts.astype(np.int64, copy=False),
            "hyper_counts": hyper_counts.astype(np.int64, copy=False),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001 - report any failure back to parent
        return {"sample": sample, "error": f"{type(exc).__name__}: {exc}"}


def _run_pool(
    samples: Dict[str, str],
    ncpus: int,
    desc: str,
) -> List[Dict[str, Any]]:
    """Dispatch ``_process_one_sample`` over ``samples``, with progress."""
    items = list(samples.items())
    results: List[Dict[str, Any]] = []
    if not items:
        return results

    progress_columns = (
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    )

    if ncpus <= 1:
        with Progress(*progress_columns, console=console) as progress:
            task = progress.add_task(desc, total=len(items))
            for item in items:
                results.append(_process_one_sample(item))
                progress.advance(task)
        return results

    # Use fork on Linux so workers inherit the parent's already-loaded
    # cpg_filter_df via copy-on-write rather than re-pickling it per worker.
    start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    mp_ctx = mp.get_context(start_method)

    with ProcessPoolExecutor(
        max_workers=ncpus,
        mp_context=mp_ctx,
        initializer=_init_worker_globals,
        initargs=(_CTX,),
    ) as pool:
        futures = {pool.submit(_process_one_sample, item): item[0] for item in items}
        with Progress(*progress_columns, console=console) as progress:
            task = progress.add_task(desc, total=len(items))
            for fut in as_completed(futures):
                results.append(fut.result())
                progress.advance(task)
    return results


# ---------------------------------------------------------------------------
# Output construction
# ---------------------------------------------------------------------------

def _stack(results: List[Dict[str, Any]], key: str) -> np.ndarray:
    """Stack a per-sample 1-D array into a 2-D matrix of shape (n_samples, n_chr)."""
    return np.vstack([r[key] for r in results]) if results else np.empty((0, 0))


def _build_dataframe(
    results: List[Dict[str, Any]],
    chr_list: List[str],
    inter_arrays: Optional[Dict[str, np.ndarray]] = None,
) -> pd.DataFrame:
    """Build a tidy wide DataFrame from per-sample results.

    Columns are constructed once, then populated column-by-column from numpy
    matrices to avoid per-row Python loops.
    """
    n = len(results)
    n_chr = len(chr_list)
    if n == 0:
        cols = ["sample"]
        for chr_name in chr_list:
            cols += [
                f"{chr_name}_hypo_beta", f"{chr_name}_hyper_beta",
                f"{chr_name}_hypo_z_intra", f"{chr_name}_hyper_z_intra",
                f"{chr_name}_s_intra",
                f"{chr_name}_hypo_cpgs_count", f"{chr_name}_hyper_cpgs_count",
            ]
            if inter_arrays is not None:
                cols += [
                    f"{chr_name}_hypo_z_inter", f"{chr_name}_hyper_z_inter",
                    f"{chr_name}_s_inter",
                ]
        return pd.DataFrame(columns=cols)

    samples = [r["sample"] for r in results]
    hypo_beta = _stack(results, "hypo_beta")
    hyper_beta = _stack(results, "hyper_beta")
    hypo_z_intra = _stack(results, "hypo_z_intra")
    hyper_z_intra = _stack(results, "hyper_z_intra")
    s_intra = _stack(results, "s_intra")
    hypo_counts = _stack(results, "hypo_counts")
    hyper_counts = _stack(results, "hyper_counts")

    data: Dict[str, np.ndarray] = {"sample": np.asarray(samples, dtype=object)}
    for idx, chr_name in enumerate(chr_list):
        data[f"{chr_name}_hypo_beta"] = hypo_beta[:, idx]
        data[f"{chr_name}_hyper_beta"] = hyper_beta[:, idx]
        data[f"{chr_name}_hypo_z_intra"] = hypo_z_intra[:, idx]
        data[f"{chr_name}_hyper_z_intra"] = hyper_z_intra[:, idx]
        data[f"{chr_name}_s_intra"] = s_intra[:, idx]
        data[f"{chr_name}_hypo_cpgs_count"] = hypo_counts[:, idx]
        data[f"{chr_name}_hyper_cpgs_count"] = hyper_counts[:, idx]
        if inter_arrays is not None:
            data[f"{chr_name}_hypo_z_inter"] = inter_arrays["hypo_z_inter"][:, idx]
            data[f"{chr_name}_hyper_z_inter"] = inter_arrays["hyper_z_inter"][:, idx]
            data[f"{chr_name}_s_inter"] = inter_arrays["s_inter"][:, idx]

    return pd.DataFrame(data)


def _write_tsv(df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame as TSV; gzip-compress when the path ends in .gz."""
    compression = "gzip" if path.endswith(".gz") else None
    df.to_csv(path, sep="\t", index=False, float_format="%.6f", compression=compression)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--analyze-beta-value-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing analyze sample {sample}_beta_value.tsv.gz files.",
)
@click.option(
    "--reference-beta-value-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing reference sample {sample}_beta_value.tsv.gz files.",
)
@click.option(
    "--output-prefix",
    required=True,
    type=str,
    help="Prefix for the two output matrices (a directory portion is allowed).",
)
@click.option(
    "--depth",
    type=int,
    default=None,
    help="Minimum depth threshold for filtering CpG sites (optional).",
)
@click.option(
    "--cpg-list",
    type=click.Path(exists=True),
    default=None,
    help="TSV file with chr/start/end columns identifying CpG sites to keep.",
)
@click.option(
    "--chr-list",
    default="1-22",
    type=str,
    help='Chromosomes to analyze (e.g., "1-22" or "1,2,3,X,Y").',
)
@click.option(
    "--beta-cols",
    default="chr,start,end,target_meth_count,target_unmeth_count,raw_total_count,meandiff",
    type=str,
    help="Comma-separated list of columns to read from beta files.",
)
@click.option(
    "--depth-col",
    default="raw_total_count",
    type=str,
    help="Column name used for depth filtering.",
)
@click.option(
    "--ncpus",
    default=max(1, (os.cpu_count() or 2) // 2),
    type=int,
    help="Number of worker processes for parallel sample processing.",
)
@click.option(
    "--no-gzip",
    is_flag=True,
    default=False,
    help="Write output as plain TSV instead of TSV.gz.",
)
def main(
    analyze_beta_value_dir: str,
    reference_beta_value_dir: str,
    output_prefix: str,
    depth: Optional[int],
    cpg_list: Optional[str],
    chr_list: str,
    beta_cols: str,
    depth_col: str,
    ncpus: int,
    no_gzip: bool,
) -> None:
    """Multi-sample beta-to-episcore conversion for trisomy detection.

    For each sample in --reference-beta-value-dir and --analyze-beta-value-dir:

    \b
    1. Read the gzipped beta TSV (only required columns).
    2. Apply chr / depth / CpG-list filters.
    3. Aggregate per-chromosome beta (hypo and hyper) and compute s_intra.

    Reference samples are aggregated into a reference s_intra matrix.
    Per-chromosome means and stds of hypo/hyper z_intra across reference
    samples are then used to compute s_inter for each analyze sample.
    """
    console.rule("[bold blue]Multi-sample Beta to Z-score Conversion")
    console.print("\n[bold]Input parameters[/bold]")
    console.print(f"  Reference dir : {reference_beta_value_dir}")
    console.print(f"  Analyze dir   : {analyze_beta_value_dir}")
    console.print(f"  Output prefix : {output_prefix}")
    console.print(f"  Depth filter  : {depth if depth is not None else 'None'}")
    console.print(f"  CpG list      : {cpg_list if cpg_list else 'None'}")
    console.print(f"  Chromosomes   : {chr_list}")
    console.print(f"  Depth column  : {depth_col}")
    console.print(f"  Workers       : {ncpus}")

    try:
        chromosomes = _parse_chr_list(chr_list)
        beta_columns = [c.strip() for c in beta_cols.split(",") if c.strip()]
        # Ensure required columns are present even if user trims aggressively.
        for required in ("chr", "start", "end", "meandiff"):
            if required not in beta_columns:
                beta_columns.append(required)
        console.print(f"  {len(chromosomes)} chromosomes  |  {len(beta_columns)} columns to read")

        # Discover samples
        ref_samples = _discover_samples(reference_beta_value_dir)
        analyze_samples = _discover_samples(analyze_beta_value_dir)
        if not ref_samples:
            raise FileNotFoundError(
                f"No *_beta_value.tsv.gz files found in {reference_beta_value_dir}"
            )
        if not analyze_samples:
            raise FileNotFoundError(
                f"No *_beta_value.tsv.gz files found in {analyze_beta_value_dir}"
            )
        console.print(
            f"\n  Found {len(ref_samples)} reference samples, "
            f"{len(analyze_samples)} analyze samples"
        )

        # Load CpG list once (shared with workers via fork copy-on-write)
        cpg_filter_df: Optional[pd.DataFrame] = None
        if cpg_list:
            console.print("\n[bold cyan]Loading CpG list[/bold cyan]")
            cpg_filter_df = pd.read_csv(
                cpg_list,
                sep="\t",
                usecols=["chr", "start", "end"],
            )
            cpg_filter_df["chr"] = cpg_filter_df["chr"].astype(str)
            cpg_filter_df["start"] = cpg_filter_df["start"].astype(np.int64)
            cpg_filter_df["end"] = cpg_filter_df["end"].astype(np.int64)
            console.print(f"[green]OK[/green] Loaded {len(cpg_filter_df):,} CpG sites")

        # Populate worker context (read by workers after fork)
        global _CTX
        _CTX = {
            "chr_list": chromosomes,
            "beta_columns": beta_columns,
            "cpg_filter_df": cpg_filter_df,
            "depth": depth,
            "depth_col": depth_col if depth is not None else None,
        }

        # 1) Process reference samples in parallel
        console.print("\n[bold cyan]Step 1: Processing reference samples[/bold cyan]")
        ref_results = _run_pool(ref_samples, ncpus, "Reference samples")
        ref_ok = [r for r in ref_results if r.get("error") is None]
        ref_err = [r for r in ref_results if r.get("error") is not None]
        for r in ref_err:
            console.print(f"[bold red]Reference failed[/bold red] {r['sample']}: {r['error']}")
        if not ref_ok:
            raise RuntimeError("All reference samples failed; cannot build reference matrix.")
        console.print(
            f"[green]OK[/green] Reference processed: "
            f"{len(ref_ok)} succeeded, {len(ref_err)} failed"
        )

        # 2) Reference statistics for s_inter normalization
        console.print("\n[bold cyan]Step 2: Computing reference statistics[/bold cyan]")
        hypo_z_stack = _stack(ref_ok, "hypo_z_intra")
        hyper_z_stack = _stack(ref_ok, "hyper_z_intra")
        with np.errstate(invalid="ignore"):
            hypo_means = np.nanmean(hypo_z_stack, axis=0)
            hypo_stds = np.nanstd(hypo_z_stack, axis=0, ddof=0)
            hyper_means = np.nanmean(hyper_z_stack, axis=0)
            hyper_stds = np.nanstd(hyper_z_stack, axis=0, ddof=0)
        # Replace NaN means/stds with safe values so s_inter stays defined.
        hypo_means = np.where(np.isfinite(hypo_means), hypo_means, 0.0)
        hyper_means = np.where(np.isfinite(hyper_means), hyper_means, 0.0)
        console.print("[green]OK[/green] Reference per-chromosome mean/std computed")

        # 3) Process analyze samples in parallel
        console.print("\n[bold cyan]Step 3: Processing analyze samples[/bold cyan]")
        analyze_results = _run_pool(analyze_samples, ncpus, "Analyze samples")
        analyze_ok = [r for r in analyze_results if r.get("error") is None]
        analyze_err = [r for r in analyze_results if r.get("error") is not None]
        for r in analyze_err:
            console.print(f"[bold red]Analyze failed[/bold red] {r['sample']}: {r['error']}")
        console.print(
            f"[green]OK[/green] Analyze processed: "
            f"{len(analyze_ok)} succeeded, {len(analyze_err)} failed"
        )

        # 4) Compute s_inter for analyze samples (vectorized across all of them)
        console.print("\n[bold cyan]Step 4: Computing s_inter for analyze samples[/bold cyan]")
        n_chr = len(chromosomes)
        n_analyze = len(analyze_ok)
        hypo_z_inter_mat = np.empty((n_analyze, n_chr), dtype=np.float64)
        hyper_z_inter_mat = np.empty((n_analyze, n_chr), dtype=np.float64)
        s_inter_mat = np.empty((n_analyze, n_chr), dtype=np.float64)
        for i, r in enumerate(analyze_ok):
            hypo_z_inter, hyper_z_inter, s_inter = calculate_s_inter_from_stats(
                r["hypo_z_intra"], r["hyper_z_intra"],
                r["hypo_counts"], r["hyper_counts"],
                hypo_means, hypo_stds, hyper_means, hyper_stds,
            )
            hypo_z_inter_mat[i] = hypo_z_inter
            hyper_z_inter_mat[i] = hyper_z_inter
            s_inter_mat[i] = s_inter
        inter_arrays = {
            "hypo_z_inter": hypo_z_inter_mat,
            "hyper_z_inter": hyper_z_inter_mat,
            "s_inter": s_inter_mat,
        }
        console.print("[green]OK[/green] s_inter computed for all analyze samples")

        # 5) Sort outputs by sample name for deterministic ordering
        ref_ok.sort(key=lambda r: r["sample"])
        # For analyze, also reorder inter matrices to match the sorted samples.
        order = sorted(range(len(analyze_ok)), key=lambda i: analyze_ok[i]["sample"])
        analyze_ok_sorted = [analyze_ok[i] for i in order]
        inter_arrays_sorted = {k: v[order] for k, v in inter_arrays.items()}

        # 6) Build and write outputs
        console.print("\n[bold cyan]Step 5: Writing output matrices[/bold cyan]")
        ref_df = _build_dataframe(ref_ok, chromosomes, inter_arrays=None)
        analyze_df = _build_dataframe(
            analyze_ok_sorted, chromosomes, inter_arrays=inter_arrays_sorted,
        )

        suffix = ".tsv" if no_gzip else ".tsv.gz"
        ref_path = f"{output_prefix}_reference_zscore{suffix}"
        analyze_path = f"{output_prefix}_analyze_zscore{suffix}"

        out_dir = os.path.dirname(ref_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        _write_tsv(ref_df, ref_path)
        _write_tsv(analyze_df, analyze_path)

        console.print(
            f"[green]OK[/green] Reference matrix: {ref_path}  "
            f"({ref_df.shape[0]} samples x {ref_df.shape[1]} columns)"
        )
        console.print(
            f"[green]OK[/green] Analyze matrix : {analyze_path}  "
            f"({analyze_df.shape[0]} samples x {analyze_df.shape[1]} columns)"
        )

        console.rule("[bold green]Analysis complete")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
