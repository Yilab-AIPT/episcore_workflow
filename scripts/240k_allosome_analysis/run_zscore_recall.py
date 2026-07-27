#!/usr/bin/env python3
"""Compute chrX read-count zscore for one recall level (cutoff=0.85).

Uses 20260416 normal females (ref_type=early_ref) to build per-chr percentage
mean/std, then scores analyze samples. Chromosome set includes chrX.

I/O note: deconv parquets are filtered to ``prob_class_1 >= cutoff`` in the
scan before collect to keep memory down.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import click
import numpy as np
import pandas as pd
import polars as pl
from rich.console import Console

console = Console()

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
ALL_CHRS = AUTOSOMES + ["chrX", "chrY"]

# Populated in worker via initializer (avoids re-reading CpG list per sample).
_CPG_POSITIONS: Dict[str, np.ndarray] = {}


def _init_worker(cpg_list: str) -> None:
    global _CPG_POSITIONS
    _CPG_POSITIONS = read_recall_positions(Path(cpg_list))


def _standardize_chr_expr() -> pl.Expr:
    """Map MQ chr codes to UCSC names.

    Upstream deconv parquets store chromosomes as Int64 with ``23`` = X
    and (rarely) ``24`` = Y. Naively prefixing yields ``chr23``, which never
    overlaps CpG lists that use ``chrX``.
    """
    s = pl.col("chr").cast(pl.Utf8)
    return (
        pl.when(s.is_in(["23", "X", "chrX"]))
        .then(pl.lit("chrX"))
        .when(s.is_in(["24", "Y", "chrY"]))
        .then(pl.lit("chrY"))
        .when(s.str.starts_with("chr"))
        .then(s)
        .otherwise(pl.lit("chr") + s)
        .alias("chr")
    )


def read_deconv_paths(paths: str, cutoff: float, mtcount: float) -> pl.DataFrame:
    """Load deconv files, keep only reads with prob_class_1 >= cutoff."""
    cols = ["chr", "start", "end", "text", "prob_class_1", "mTcount"]
    frames: List[pl.DataFrame] = []
    for raw in paths.split(","):
        path = Path(raw.strip())
        if not path.is_file():
            raise FileNotFoundError(f"Missing deconv file: {path}")
        suffix = "".join(path.suffixes).lower()
        if ".parquet" in suffix or ".pq" in suffix:
            lf = pl.scan_parquet(path)
        else:
            lf = pl.scan_csv(path, separator="\t")
        available = lf.collect_schema().names()
        missing = [c for c in cols if c not in available]
        if missing:
            raise ValueError(f"{path} missing columns {missing}")
        frames.append(
            lf.select(cols)
            .filter(
                (pl.col("prob_class_1") >= cutoff) & (pl.col("mTcount") >= mtcount)
            )
            .with_columns(
                pl.col("start").cast(pl.Int64),
                pl.col("end").cast(pl.Int64),
                pl.col("text").cast(pl.Utf8),
                pl.col("prob_class_1").cast(pl.Float64),
                pl.col("mTcount").cast(pl.Float64),
            )
            .with_columns(_standardize_chr_expr())
            .drop_nulls(["chr", "start", "end", "prob_class_1"])
            .collect()
        )
    df = pl.concat(frames, how="vertical_relaxed")
    return df.unique(subset=["chr", "start", "end", "text"], keep="first")


def read_recall_positions(cpg_list: Path) -> Dict[str, np.ndarray]:
    cdf = pl.read_csv(cpg_list, separator="\t", columns=["chr", "start"]).with_columns(
        pl.col("chr").cast(pl.Utf8), pl.col("start").cast(pl.Int64)
    )
    cdf = cdf.with_columns(
        pl.when(pl.col("chr").str.starts_with("chr"))
        .then(pl.col("chr"))
        .otherwise(pl.lit("chr") + pl.col("chr"))
        .alias("chr")
    )
    out: Dict[str, np.ndarray] = {}
    for chrom, grp in cdf.group_by("chr"):
        name = str(chrom[0] if isinstance(chrom, tuple) else chrom)
        out[name] = np.sort(grp.get_column("start").to_numpy())
    return out


def overlap_mask(starts: np.ndarray, ends: np.ndarray, cpg_starts: np.ndarray) -> np.ndarray:
    if cpg_starts.size == 0:
        return np.zeros(starts.shape[0], dtype=bool)
    left = np.searchsorted(cpg_starts, starts, side="left")
    right = np.searchsorted(cpg_starts, ends + 1, side="left")
    return right > left


def chromosome_percentages(
    df: pl.DataFrame,
    cpg_positions: Dict[str, np.ndarray],
) -> Dict[str, float]:
    counts = {c: 0 for c in ALL_CHRS}
    if df.height == 0:
        return {c: float("nan") for c in ALL_CHRS}

    for chrom, grp in df.group_by("chr"):
        name = str(chrom[0] if isinstance(chrom, tuple) else chrom)
        if name not in counts:
            continue
        starts = grp.get_column("start").to_numpy()
        ends = grp.get_column("end").to_numpy()
        cpg = cpg_positions.get(name, np.empty(0, dtype=np.int64))
        counts[name] = int(overlap_mask(starts, ends, cpg).sum())

    sum_auto = sum(counts[c] for c in AUTOSOMES)
    if sum_auto <= 0:
        return {c: float("nan") for c in ALL_CHRS}
    return {c: counts[c] / float(sum_auto) for c in ALL_CHRS}


def _worker(args: Tuple[str, str, float]) -> Tuple[str, Dict[str, float]]:
    sample, deconv_paths, cutoff = args
    df = read_deconv_paths(deconv_paths, cutoff=cutoff, mtcount=1.0)
    pct = chromosome_percentages(df, _CPG_POSITIONS)
    return sample, pct


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--samples-meta", required=True, type=click.Path(exists=True))
@click.option("--cpg-list", required=True, type=click.Path(exists=True))
@click.option("--output-prefix", required=True, type=str)
@click.option("--cutoff", type=float, default=0.85, show_default=True)
@click.option("--ncpus", type=int, default=4, show_default=True)
def main(
    samples_meta: str,
    cpg_list: str,
    output_prefix: str,
    cutoff: float,
    ncpus: int,
) -> None:
    """Write reference percentages + analyze chrX zscores for one recall."""
    meta = pd.read_csv(samples_meta)
    needed = {"sample", "deconv_paths", "ref_type"}
    missing = needed - set(meta.columns)
    if missing:
        raise click.ClickException(f"samples-meta missing: {sorted(missing)}")

    ref = meta.loc[meta["ref_type"] == "early_ref"].copy()
    analyze = meta.loc[meta["ref_type"] == "analyze"].copy()
    if ref.empty:
        raise click.ClickException("No early_ref samples in samples-meta")
    if analyze.empty:
        raise click.ClickException("No analyze samples in samples-meta")

    score_samples = pd.concat([ref, analyze], ignore_index=True)

    console.print(f"Recall CpG list : {cpg_list}")
    console.print(f"Cutoff          : {cutoff}")
    console.print(f"Reference n     : {len(ref)}")
    console.print(f"Analyze n       : {len(analyze)} (+{len(ref)} ref scored)")
    console.print(f"Workers         : {ncpus}")

    tasks = [
        (str(r.sample), str(r.deconv_paths), float(cutoff))
        for r in score_samples.itertuples(index=False)
    ]

    results: Dict[str, Dict[str, float]] = {}
    with ProcessPoolExecutor(
        max_workers=ncpus,
        initializer=_init_worker,
        initargs=(str(cpg_list),),
    ) as pool:
        futures = {pool.submit(_worker, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            sample, pct = fut.result()
            results[sample] = pct
            console.print(f"  done {sample}")

    ref_rows = []
    for sample in ref["sample"]:
        pct = results[sample]
        row = {"sample": sample, **{f"{c}_percentage": pct[c] for c in ALL_CHRS}}
        ref_rows.append(row)
    ref_df = pd.DataFrame(ref_rows)

    stats = {}
    for c in ALL_CHRS:
        col = f"{c}_percentage"
        vals = ref_df[col].to_numpy(dtype=float)
        stats[c] = (float(np.nanmean(vals)), float(np.nanstd(vals, ddof=1)))

    analyze_rows = []
    for sample in score_samples["sample"]:
        pct = results[sample]
        row = {"sample": sample}
        for c in ALL_CHRS:
            mean, std = stats[c]
            p = pct[c]
            row[f"{c}_percentage"] = p
            if std > 0 and np.isfinite(p):
                row[f"{c}_zscore"] = (p - mean) / std
            else:
                row[f"{c}_zscore"] = 0.0
        analyze_rows.append(row)
    analyze_df = pd.DataFrame(analyze_rows)

    out_dir = Path(output_prefix.rstrip("/"))
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_path = out_dir / "_reference_percentage.tsv.gz"
    an_path = out_dir / "_analyze_zscore.tsv.gz"
    ref_df.to_csv(ref_path, sep="\t", index=False, float_format="%.6f", compression="gzip")
    analyze_df.to_csv(an_path, sep="\t", index=False, float_format="%.6f", compression="gzip")
    console.print(f"[green]Wrote[/green] {ref_path}")
    console.print(f"[green]Wrote[/green] {an_path}")


if __name__ == "__main__":
    main()
