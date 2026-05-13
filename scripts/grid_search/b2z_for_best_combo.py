#!/usr/bin/env python3
"""
Best-combo per-chromosome Z-score recomputation.

For each chromosome, a "best (threshold, recall) combination" is read from a
CSV. This script picks up the per-chromosome aggregated beta values and CpG
counts that ``beta_to_zscore.py`` already wrote for those combinations, then
recomputes per-sample ``hypo_z_intra`` / ``hyper_z_intra`` (across chromosomes
using each chr's best combo), the reference statistics, and finally
``hypo_z_inter`` / ``hyper_z_inter`` / ``s_inter`` for every analyze sample.

Inputs:
    --best-combo-csv : CSV with columns ``chr,threshold,recall`` (extra columns
                        such as ``has_target`` are ignored).
    --output-base    : Directory containing the
                        ``threshold_{t}_recall_{r}/_analyze_zscore.tsv.gz`` and
                        ``_reference_zscore.tsv.gz`` files produced by
                        ``beta_to_zscore.py``.

Outputs (under ``--output-base/best_combo/`` by default):
    _reference_zscore.tsv.gz
    _analyze_zscore.tsv.gz

The output schemas exactly match those of ``beta_to_zscore.py``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

from beta_to_zscore import (
    _build_dataframe,
    _write_tsv,
    calculate_s_inter_from_stats,
    calculate_s_intra,
)

console = Console()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_combo_dirname(threshold: float, recall: float) -> str:
    """Format a (threshold, recall) combo as ``threshold_{t}_recall_{r}``.

    ``%g`` matches the convention used by ``submit_beta_to_zscore.sh`` (it strips
    trailing zeros, e.g. ``0.30`` -> ``0.3``) so the directory names line up
    with the ones already on disk.
    """
    return f"threshold_{threshold:g}_recall_{recall:g}"


def _load_best_combo(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"chr", "threshold", "recall"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"best-combo CSV is missing columns: {sorted(missing)}")
    df = df[["chr", "threshold", "recall"]].copy()
    df["chr"] = df["chr"].astype(str)
    df["threshold"] = df["threshold"].astype(float)
    df["recall"] = df["recall"].astype(float)
    if df["chr"].duplicated().any():
        dups = df["chr"][df["chr"].duplicated()].tolist()
        raise ValueError(f"best-combo CSV has duplicate chr rows: {dups}")
    return df


def _load_combo_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing combo output: {path}")
    return pd.read_csv(path, sep="\t", compression="gzip")


def _compute_zscores_per_sample(
    hypo_beta: np.ndarray,
    hyper_beta: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply ``calculate_s_intra`` row-wise for each sample."""
    n = hypo_beta.shape[0]
    hypo_z = np.empty_like(hypo_beta, dtype=np.float64)
    hyper_z = np.empty_like(hyper_beta, dtype=np.float64)
    s_intra = np.empty_like(hypo_beta, dtype=np.float64)
    for i in range(n):
        hz, hyz, si = calculate_s_intra(
            hypo_beta[i], hyper_beta[i], hypo_counts[i], hyper_counts[i],
        )
        hypo_z[i] = hz
        hyper_z[i] = hyz
        s_intra[i] = si
    return hypo_z, hyper_z, s_intra


def _build_results_list(
    samples: List[str],
    hypo_beta: np.ndarray,
    hyper_beta: np.ndarray,
    hypo_z: np.ndarray,
    hyper_z: np.ndarray,
    s_intra: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
) -> List[Dict[str, np.ndarray]]:
    return [
        {
            "sample": sample,
            "hypo_beta": hypo_beta[i],
            "hyper_beta": hyper_beta[i],
            "hypo_z_intra": hypo_z[i],
            "hyper_z_intra": hyper_z[i],
            "s_intra": s_intra[i],
            "hypo_counts": hypo_counts[i],
            "hyper_counts": hyper_counts[i],
        }
        for i, sample in enumerate(samples)
    ]


# ---------------------------------------------------------------------------
# Core assembly: pull per-chr arrays from each combo's output
# ---------------------------------------------------------------------------

def _assemble_best_combo_arrays(
    output_base: Path,
    combo_df: pd.DataFrame,
) -> Tuple[
    List[str],          # chr_list (in CSV order)
    List[str],          # analyze samples
    List[str],          # reference samples
    np.ndarray,         # analyze hypo_beta (n_analyze, n_chr)
    np.ndarray,         # analyze hyper_beta
    np.ndarray,         # analyze hypo_counts
    np.ndarray,         # analyze hyper_counts
    np.ndarray,         # reference hypo_beta (n_reference, n_chr)
    np.ndarray,         # reference hyper_beta
    np.ndarray,         # reference hypo_counts
    np.ndarray,         # reference hyper_counts
]:
    chr_list = combo_df["chr"].tolist()
    chr_index = {c: i for i, c in enumerate(chr_list)}

    combo_to_chrs: Dict[Tuple[float, float], List[str]] = {}
    for row in combo_df.itertuples(index=False):
        combo_to_chrs.setdefault((row.threshold, row.recall), []).append(row.chr)

    # Anchor the sample order to the first combo encountered; subsequent
    # combos are aligned via ``reindex(samples)``.
    analyze_samples: List[str] = []
    reference_samples: List[str] = []

    n_chr = len(chr_list)
    analyze_hypo_beta = np.full((0, n_chr), np.nan, dtype=np.float64)
    analyze_hyper_beta = np.full((0, n_chr), np.nan, dtype=np.float64)
    analyze_hypo_counts = np.zeros((0, n_chr), dtype=np.int64)
    analyze_hyper_counts = np.zeros((0, n_chr), dtype=np.int64)

    reference_hypo_beta = np.full((0, n_chr), np.nan, dtype=np.float64)
    reference_hyper_beta = np.full((0, n_chr), np.nan, dtype=np.float64)
    reference_hypo_counts = np.zeros((0, n_chr), dtype=np.int64)
    reference_hyper_counts = np.zeros((0, n_chr), dtype=np.int64)

    initialized = False

    for combo_idx, ((thres, recall), chrs) in enumerate(combo_to_chrs.items(), start=1):
        combo_dir = output_base / _format_combo_dirname(thres, recall)
        analyze_path = combo_dir / "_analyze_zscore.tsv.gz"
        reference_path = combo_dir / "_reference_zscore.tsv.gz"

        console.print(
            f"  [{combo_idx}/{len(combo_to_chrs)}] {combo_dir.name}  "
            f"-> {len(chrs)} chr(s): {', '.join(chrs)}"
        )

        a_df = _load_combo_table(analyze_path).set_index("sample")
        r_df = _load_combo_table(reference_path).set_index("sample")

        if not initialized:
            analyze_samples = a_df.index.astype(str).tolist()
            reference_samples = r_df.index.astype(str).tolist()
            n_a = len(analyze_samples)
            n_r = len(reference_samples)

            analyze_hypo_beta = np.full((n_a, n_chr), np.nan, dtype=np.float64)
            analyze_hyper_beta = np.full((n_a, n_chr), np.nan, dtype=np.float64)
            analyze_hypo_counts = np.zeros((n_a, n_chr), dtype=np.int64)
            analyze_hyper_counts = np.zeros((n_a, n_chr), dtype=np.int64)
            reference_hypo_beta = np.full((n_r, n_chr), np.nan, dtype=np.float64)
            reference_hyper_beta = np.full((n_r, n_chr), np.nan, dtype=np.float64)
            reference_hypo_counts = np.zeros((n_r, n_chr), dtype=np.int64)
            reference_hyper_counts = np.zeros((n_r, n_chr), dtype=np.int64)
            initialized = True

        a_df = a_df.reindex(analyze_samples)
        r_df = r_df.reindex(reference_samples)

        for chr_name in chrs:
            j = chr_index[chr_name]
            for col, dtype, into in (
                (f"{chr_name}_hypo_beta", np.float64, analyze_hypo_beta),
                (f"{chr_name}_hyper_beta", np.float64, analyze_hyper_beta),
                (f"{chr_name}_hypo_cpgs_count", np.int64, analyze_hypo_counts),
                (f"{chr_name}_hyper_cpgs_count", np.int64, analyze_hyper_counts),
            ):
                if col not in a_df.columns:
                    raise KeyError(f"Column '{col}' missing in {analyze_path}")
                arr = a_df[col].to_numpy()
                if np.issubdtype(dtype, np.integer):
                    # Counts can come back as float when the row is NaN due to
                    # reindex; replace with 0 so the int cast is safe.
                    arr = np.nan_to_num(arr, nan=0).astype(dtype, copy=False)
                else:
                    arr = arr.astype(dtype, copy=False)
                into[:, j] = arr

            for col, dtype, into in (
                (f"{chr_name}_hypo_beta", np.float64, reference_hypo_beta),
                (f"{chr_name}_hyper_beta", np.float64, reference_hyper_beta),
                (f"{chr_name}_hypo_cpgs_count", np.int64, reference_hypo_counts),
                (f"{chr_name}_hyper_cpgs_count", np.int64, reference_hyper_counts),
            ):
                if col not in r_df.columns:
                    raise KeyError(f"Column '{col}' missing in {reference_path}")
                arr = r_df[col].to_numpy()
                if np.issubdtype(dtype, np.integer):
                    arr = np.nan_to_num(arr, nan=0).astype(dtype, copy=False)
                else:
                    arr = arr.astype(dtype, copy=False)
                into[:, j] = arr

    if not initialized:
        raise RuntimeError("best-combo CSV has no rows; nothing to assemble.")

    return (
        chr_list,
        analyze_samples,
        reference_samples,
        analyze_hypo_beta,
        analyze_hyper_beta,
        analyze_hypo_counts,
        analyze_hyper_counts,
        reference_hypo_beta,
        reference_hyper_beta,
        reference_hypo_counts,
        reference_hyper_counts,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--best-combo-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="CSV with columns chr,threshold,recall identifying the best combo per chr.",
)
@click.option(
    "--output-base",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Directory containing the threshold_{t}_recall_{r}/ subdirectories "
        "produced by beta_to_zscore.py."
    ),
)
@click.option(
    "--output-subdir",
    default="best_combo",
    type=str,
    help="Subdirectory under --output-base where outputs are written.",
)
@click.option(
    "--no-gzip",
    is_flag=True,
    default=False,
    help="Write outputs as plain TSV instead of TSV.gz.",
)
def main(
    best_combo_csv: str,
    output_base: str,
    output_subdir: str,
    no_gzip: bool,
) -> None:
    """Recompute z_intra / z_inter / s_inter using the best (threshold, recall) per chr.

    \b
    1. Read the best-combo CSV.
    2. For each unique (threshold, recall) referenced by the CSV, load the
       analyze and reference TSVs and pull the chr-level beta values and
       CpG counts for the chrs that pick that combo.
    3. Recompute per-sample z_intra (across chrs), then reference statistics,
       then z_inter and s_inter for every analyze sample.
    4. Write best_combo/_reference_zscore.tsv.gz and best_combo/_analyze_zscore.tsv.gz.
    """
    console.rule("[bold blue]Best-combo Z-score recomputation")

    output_base_path = Path(output_base)
    out_dir = output_base_path / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print("\n[bold]Input parameters[/bold]")
    console.print(f"  Best-combo CSV : {best_combo_csv}")
    console.print(f"  Output base    : {output_base_path}")
    console.print(f"  Output dir     : {out_dir}")

    try:
        combo_df = _load_best_combo(best_combo_csv)
        console.print(
            f"  {len(combo_df)} chromosomes, "
            f"{combo_df[['threshold', 'recall']].drop_duplicates().shape[0]} "
            f"unique combos"
        )

        console.print("\n[bold cyan]Step 1: Assembling per-chr arrays from combo outputs[/bold cyan]")
        (
            chr_list,
            analyze_samples,
            reference_samples,
            a_hypo_beta,
            a_hyper_beta,
            a_hypo_counts,
            a_hyper_counts,
            r_hypo_beta,
            r_hyper_beta,
            r_hypo_counts,
            r_hyper_counts,
        ) = _assemble_best_combo_arrays(output_base_path, combo_df)
        console.print(
            f"[green]OK[/green] {len(analyze_samples)} analyze samples, "
            f"{len(reference_samples)} reference samples, "
            f"{len(chr_list)} chromosomes"
        )

        console.print("\n[bold cyan]Step 2: Recomputing reference z_intra & statistics[/bold cyan]")
        r_hypo_z, r_hyper_z, r_s_intra = _compute_zscores_per_sample(
            r_hypo_beta, r_hyper_beta, r_hypo_counts, r_hyper_counts,
        )
        with np.errstate(invalid="ignore"):
            hypo_means = np.nanmean(r_hypo_z, axis=0)
            hypo_stds = np.nanstd(r_hypo_z, axis=0, ddof=0)
            hyper_means = np.nanmean(r_hyper_z, axis=0)
            hyper_stds = np.nanstd(r_hyper_z, axis=0, ddof=0)
        hypo_means = np.where(np.isfinite(hypo_means), hypo_means, 0.0)
        hyper_means = np.where(np.isfinite(hyper_means), hyper_means, 0.0)
        console.print("[green]OK[/green] Reference per-chromosome mean/std computed")

        console.print("\n[bold cyan]Step 3: Recomputing analyze z_intra & s_inter[/bold cyan]")
        a_hypo_z, a_hyper_z, a_s_intra = _compute_zscores_per_sample(
            a_hypo_beta, a_hyper_beta, a_hypo_counts, a_hyper_counts,
        )

        n_analyze = len(analyze_samples)
        n_chr = len(chr_list)
        a_hypo_z_inter = np.empty((n_analyze, n_chr), dtype=np.float64)
        a_hyper_z_inter = np.empty((n_analyze, n_chr), dtype=np.float64)
        a_s_inter = np.empty((n_analyze, n_chr), dtype=np.float64)
        for i in range(n_analyze):
            hz_inter, hyz_inter, si_inter = calculate_s_inter_from_stats(
                a_hypo_z[i], a_hyper_z[i],
                a_hypo_counts[i], a_hyper_counts[i],
                hypo_means, hypo_stds, hyper_means, hyper_stds,
            )
            a_hypo_z_inter[i] = hz_inter
            a_hyper_z_inter[i] = hyz_inter
            a_s_inter[i] = si_inter
        console.print("[green]OK[/green] s_inter computed for all analyze samples")

        console.print("\n[bold cyan]Step 4: Sorting and writing outputs[/bold cyan]")
        # Sort by sample for deterministic output ordering, mirroring beta_to_zscore.py.
        ref_order = sorted(range(len(reference_samples)), key=lambda i: reference_samples[i])
        an_order = sorted(range(n_analyze), key=lambda i: analyze_samples[i])

        ref_results = _build_results_list(
            [reference_samples[i] for i in ref_order],
            r_hypo_beta[ref_order], r_hyper_beta[ref_order],
            r_hypo_z[ref_order], r_hyper_z[ref_order],
            r_s_intra[ref_order],
            r_hypo_counts[ref_order], r_hyper_counts[ref_order],
        )
        an_results = _build_results_list(
            [analyze_samples[i] for i in an_order],
            a_hypo_beta[an_order], a_hyper_beta[an_order],
            a_hypo_z[an_order], a_hyper_z[an_order],
            a_s_intra[an_order],
            a_hypo_counts[an_order], a_hyper_counts[an_order],
        )
        inter_arrays_sorted = {
            "hypo_z_inter": a_hypo_z_inter[an_order],
            "hyper_z_inter": a_hyper_z_inter[an_order],
            "s_inter": a_s_inter[an_order],
        }

        ref_df = _build_dataframe(ref_results, chr_list, inter_arrays=None)
        analyze_df = _build_dataframe(an_results, chr_list, inter_arrays=inter_arrays_sorted)

        suffix = ".tsv" if no_gzip else ".tsv.gz"
        ref_path = str(out_dir / f"_reference_zscore{suffix}")
        analyze_path = str(out_dir / f"_analyze_zscore{suffix}")

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

        console.rule("[bold green]Best-combo recomputation complete")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
