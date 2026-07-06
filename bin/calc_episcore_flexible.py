#!/usr/bin/env python3
"""Flexible per-chromosome episcore for a single sample (AIPT ref-40).

Given a set of per-threshold target bedGraphs (one ``MethylDackel`` extract per
probability threshold) and the grid-search ``best_combo_episcore.csv`` (a
``(threshold, recall)`` combo per chromosome), this script reproduces the
episcore computation used by ``scripts/ref_explore_plus_grid_search`` for a new
sample:

    1. For every distinct ``(threshold, recall)`` combo referenced by the
       best-combo table, take the matching threshold's target bedGraph, filter
       it to the recall CpG list, split hypo/hyper by ``meandiff`` and compute
       per-chromosome beta values, CpG counts and the within-sample
       ``hypo``/``hyper`` ``z_intra`` (across all 22 autosomes for that combo).
    2. For each chromosome, keep the ``z_intra`` / CpG counts coming from *that
       chromosome's* best combo.
    3. Normalise each chromosome against the reference statistics in
       ``best_reference_matrix.tsv`` (``hypo_z_intra_mean/std`` and
       ``hyper_z_intra_mean/std``) and combine hypo/hyper with sqrt-count
       weights to obtain ``s_inter`` == episcore.

The recall CpG lists are nested (a higher recall is a subset of a lower one), so
each recall file is loaded at most once and cached across thresholds.
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_LIST)}

_THR_RE = re.compile(r"thr_([0-9.]+)_target")


def _fmt(value: float) -> str:
    """Format a float the same way the file names are generated (``%g``)."""
    return f"{float(value):g}"


def load_best_combo(path: Path) -> pd.DataFrame:
    """Load the per-chr best ``(threshold, recall)`` combo table."""
    df = pd.read_csv(path)
    for col in ("chr", "threshold", "recall"):
        if col not in df.columns:
            raise ValueError(f"best-combo CSV missing column: {col}")
    df = df[["chr", "threshold", "recall"]].copy()
    df["chr"] = df["chr"].astype(str)
    df["threshold"] = df["threshold"].astype(float)
    df["recall"] = df["recall"].astype(float)
    return df


def map_bedgraphs(bedgraphs: Tuple[str, ...]) -> Dict[str, Path]:
    """Map each threshold (``%g`` string) to its target bedGraph path."""
    out: Dict[str, Path] = {}
    for bg in bedgraphs:
        p = Path(bg)
        m = _THR_RE.search(p.name)
        if not m:
            raise ValueError(
                f"Cannot parse threshold from bedGraph name '{p.name}' "
                "(expected '...thr_<t>_target...')."
            )
        out[_fmt(float(m.group(1)))] = p
    return out


def read_bedgraph(path: Path) -> pd.DataFrame:
    """Read a MethylDackel CpG bedGraph into chr/start/end/meth/unmeth.

    The ``end`` coordinate is shifted by -1 to match the CpG-list convention
    used by ``extract_beta_value.py``.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=["chr", "start", "end", "beta", "meth_count", "unmeth_count"],
        dtype={
            "chr": str,
            "start": np.int64,
            "end": np.int64,
            "beta": float,
            "meth_count": np.int64,
            "unmeth_count": np.int64,
        },
    )
    df["end"] = df["end"] - 1
    return df[["chr", "start", "end", "meth_count", "unmeth_count"]]


def read_recall_list(cpg_recall_dir: Path, recall: float) -> pd.DataFrame:
    """Load the recall CpG list (chr/start/end/meandiff) for ``recall``."""
    path = cpg_recall_dir / f"220k_cpg_recall_{_fmt(recall)}.txt"
    if not path.is_file():
        raise FileNotFoundError(f"Missing recall CpG list: {path}")
    df = pd.read_csv(path, sep="\t", usecols=["chr", "start", "end", "meandiff"])
    df["chr"] = df["chr"].astype(str)
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df["meandiff"] = df["meandiff"].astype(float)
    return df


def _aggregate_chr(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-chr beta and CpG-count arrays aligned to ``CHR_LIST``."""
    if df.empty:
        return np.full(len(CHR_LIST), np.nan), np.zeros(len(CHR_LIST), dtype=np.int64)
    grouped = df.groupby("chr", sort=False)
    meth = grouped["meth_count"].sum().reindex(CHR_LIST, fill_value=0).to_numpy(dtype=np.float64)
    unmeth = grouped["unmeth_count"].sum().reindex(CHR_LIST, fill_value=0).to_numpy(dtype=np.float64)
    counts = grouped.size().reindex(CHR_LIST, fill_value=0).to_numpy(dtype=np.int64)
    denom = meth + unmeth
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(denom > 0, meth / denom, np.nan)
    return beta, counts


def _zscore_across(arr: np.ndarray) -> np.ndarray:
    """NaN-safe standardisation of a 1-D array across chromosomes."""
    if np.all(np.isnan(arr)):
        return np.full_like(arr, np.nan, dtype=np.float64)
    mean = np.nanmean(arr)
    std = np.nanstd(arr)
    if not np.isfinite(std) or std == 0:
        return np.full_like(arr, np.nan, dtype=np.float64)
    return (arr - mean) / std


def compute_combo_z_intra(
    bedgraph: pd.DataFrame,
    recall_df: pd.DataFrame,
    depth: Optional[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-chr hypo/hyper ``z_intra`` and CpG counts for one combo."""
    merged = bedgraph.merge(recall_df, on=["chr", "start", "end"], how="inner")
    if depth is not None:
        merged = merged[(merged["meth_count"] + merged["unmeth_count"]) > depth]
    hypo = merged[merged["meandiff"] < 0]
    hyper = merged[merged["meandiff"] > 0]
    hypo_beta, hypo_counts = _aggregate_chr(hypo)
    hyper_beta, hyper_counts = _aggregate_chr(hyper)
    hypo_z = _zscore_across(hypo_beta)
    hyper_z = _zscore_across(hyper_beta)
    return hypo_z, hyper_z, hypo_counts, hyper_counts


def load_reference_matrix(path: Path) -> pd.DataFrame:
    """Load ``best_reference_matrix.tsv`` indexed by chromosome."""
    df = pd.read_csv(path, sep="\t")
    needed = {
        "chr",
        "hypo_z_intra_mean",
        "hypo_z_intra_std",
        "hyper_z_intra_mean",
        "hyper_z_intra_std",
    }
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Reference matrix missing columns: {sorted(missing)}")
    df["chr"] = df["chr"].astype(str)
    return df.set_index("chr")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--bedgraph", "bedgraphs", multiple=True, required=True,
              type=click.Path(exists=True),
              help="Per-threshold target CpG bedGraph (repeatable); name must contain 'thr_<t>_target'.")
@click.option("--best-combo-episcore", required=True, type=click.Path(exists=True),
              help="best_combo_episcore.csv with columns chr,threshold,recall.")
@click.option("--reference-matrix", required=True, type=click.Path(exists=True),
              help="best_reference_matrix.tsv with per-chr hypo/hyper z_intra mean/std.")
@click.option("--cpg-recall-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory with 220k_cpg_recall_{recall}.txt files.")
@click.option("--output-prefix", required=True, type=str,
              help="Output prefix; writes {prefix}_episcore.tsv.")
@click.option("--depth", type=int, default=None,
              help="Optional minimum (meth+unmeth) depth per CpG.")
def main(
    bedgraphs: Tuple[str, ...],
    best_combo_episcore: str,
    reference_matrix: str,
    cpg_recall_dir: str,
    output_prefix: str,
    depth: Optional[int],
) -> None:
    """Compute flexible per-chromosome episcore for a single sample."""
    try:
        console.rule("[bold blue]Flexible episcore")
        combo_df = load_best_combo(Path(best_combo_episcore))
        ref_mat = load_reference_matrix(Path(reference_matrix))
        bg_map = map_bedgraphs(bedgraphs)
        recall_dir = Path(cpg_recall_dir)

        console.print(f"  bedGraphs (thresholds) : {sorted(bg_map.keys())}")
        n_combos = combo_df[["threshold", "recall"]].drop_duplicates().shape[0]
        console.print(f"  chromosomes / combos   : {len(combo_df)} / {n_combos}")

        # group chromosomes by their (threshold, recall) combo
        combo_to_chrs: Dict[Tuple[float, float], List[str]] = {}
        for row in combo_df.itertuples(index=False):
            combo_to_chrs.setdefault((row.threshold, row.recall), []).append(row.chr)

        # cache bedgraph reads (per threshold) and recall lists (per recall)
        bedgraph_cache: Dict[str, pd.DataFrame] = {}
        recall_cache: Dict[str, pd.DataFrame] = {}

        hypo_z = np.full(len(CHR_LIST), np.nan)
        hyper_z = np.full(len(CHR_LIST), np.nan)
        hypo_counts = np.zeros(len(CHR_LIST), dtype=np.int64)
        hyper_counts = np.zeros(len(CHR_LIST), dtype=np.int64)

        for (thr, rec), chrs in combo_to_chrs.items():
            thr_key = _fmt(thr)
            rec_key = _fmt(rec)
            if thr_key not in bg_map:
                raise ValueError(
                    f"No bedGraph supplied for threshold {thr_key} "
                    f"(needed by chromosomes {chrs})."
                )
            if thr_key not in bedgraph_cache:
                bedgraph_cache[thr_key] = read_bedgraph(bg_map[thr_key])
            if rec_key not in recall_cache:
                recall_cache[rec_key] = read_recall_list(recall_dir, rec)

            combo_hypo_z, combo_hyper_z, combo_hypo_c, combo_hyper_c = compute_combo_z_intra(
                bedgraph_cache[thr_key], recall_cache[rec_key], depth
            )
            for chrom in chrs:
                j = CHR_INDEX[chrom]
                hypo_z[j] = combo_hypo_z[j]
                hyper_z[j] = combo_hyper_z[j]
                hypo_counts[j] = combo_hypo_c[j]
                hyper_counts[j] = combo_hyper_c[j]

        # reference normalisation (s_inter) per chromosome
        ref = ref_mat.reindex(CHR_LIST)
        hypo_mean = ref["hypo_z_intra_mean"].to_numpy(dtype=np.float64)
        hypo_std = ref["hypo_z_intra_std"].to_numpy(dtype=np.float64)
        hyper_mean = ref["hyper_z_intra_mean"].to_numpy(dtype=np.float64)
        hyper_std = ref["hyper_z_intra_std"].to_numpy(dtype=np.float64)

        with np.errstate(divide="ignore", invalid="ignore"):
            hypo_std_safe = np.where(hypo_std > 0, hypo_std, np.nan)
            hyper_std_safe = np.where(hyper_std > 0, hyper_std, np.nan)
            hypo_z_inter = (hypo_z - hypo_mean) / hypo_std_safe
            hyper_z_inter = (hyper_z - hyper_mean) / hyper_std_safe

        w_hypo = np.sqrt(hypo_counts.astype(np.float64))
        w_hyper = np.sqrt(hyper_counts.astype(np.float64))
        total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            s_inter = (hyper_z_inter * w_hyper - hypo_z_inter * w_hypo) / total_w
            s_inter = np.where(np.isnan(s_inter), 0.0, s_inter)

        out_df = pd.DataFrame({"chr": CHR_LIST, "episcore": s_inter})
        out_path = f"{output_prefix}_episcore.tsv"
        out_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
        console.print(f"[green]OK[/green] Wrote {out_path}")
        console.rule("[bold green]Done")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
