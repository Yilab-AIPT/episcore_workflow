#!/usr/bin/env python3
"""Flexible per-chromosome read-count zscore for a single sample (AIPT ref-40).

Reproduces the read-count z-score path of ``workflow/zscore`` for a new sample,
driven by the grid-search ``best_combo_zscore.csv``. For each distinct
``(threshold, recall)`` combo:

    1. Keep deconv reads with ``prob_class_1 >= threshold`` and
       ``mTcount >= --mtcount`` (deduplicated on chr/start/end/text, matching
       ``mq_zscore_analyzer.py``).
    2. Keep reads overlapping a CpG in the recall list (the recall ``bed`` filter).
    3. Count reads per chromosome and compute the percentage
       ``readscount / sum(autosome readscounts)``.

Each chromosome then takes its percentage at *its* best combo and is normalised
against the reference statistics in ``best_reference_matrix.tsv``
(``percentage_mean`` / ``percentage_std``) to produce the zscore.

The deconv table is read once; recall CpG positions are cached per recall.
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import click
import numpy as np
import polars as pl
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_LIST)}


def _fmt(value: float) -> str:
    return f"{float(value):g}"


def load_best_combo(path: Path) -> "pl.DataFrame":
    import pandas as pd

    df = pd.read_csv(path)
    for col in ("chr", "threshold", "recall"):
        if col not in df.columns:
            raise ValueError(f"best-combo CSV missing column: {col}")
    df = df[["chr", "threshold", "recall"]].copy()
    df["chr"] = df["chr"].astype(str)
    df["threshold"] = df["threshold"].astype(float)
    df["recall"] = df["recall"].astype(float)
    return df


def load_reference_matrix(path: Path):
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    needed = {"chr", "percentage_mean", "percentage_std"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Reference matrix missing columns: {sorted(missing)}")
    df["chr"] = df["chr"].astype(str)
    return df.set_index("chr")


def read_deconv(path: Path, mtcount: float) -> pl.DataFrame:
    """Load deconv reads, standardise chr, dedup on chr/start/end/text."""
    cols = ["chr", "start", "end", "text", "prob_class_1", "mTcount"]
    suffix = "".join(path.suffixes).lower()
    if ".parquet" in suffix or ".pq" in suffix:
        lf = pl.scan_parquet(path)
    else:
        lf = pl.scan_csv(path, separator="\t")
    available = lf.collect_schema().names()
    missing = [c for c in cols if c not in available]
    if missing:
        raise ValueError(f"Deconv file missing column(s) {missing}; available: {available}")
    df = (
        lf.select(cols)
        .with_columns(
            pl.col("chr").cast(pl.Utf8),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64),
            pl.col("text").cast(pl.Utf8),
            pl.col("prob_class_1").cast(pl.Float64),
            pl.col("mTcount").cast(pl.Float64),
        )
        .drop_nulls(["chr", "start", "end", "prob_class_1"])
        .collect()
    )
    df = df.unique(subset=["chr", "start", "end", "text"], keep="first")
    # standardise chromosome naming -> chr-prefixed
    df = df.with_columns(
        pl.when(pl.col("chr").str.starts_with("chr"))
        .then(pl.col("chr"))
        .otherwise(pl.lit("chr") + pl.col("chr"))
        .alias("chr")
    )
    return df


def read_recall_positions(cpg_recall_dir: Path, recall: float) -> Dict[str, np.ndarray]:
    """Return per-chr sorted CpG start positions for the recall list."""
    path = cpg_recall_dir / f"220k_cpg_recall_{_fmt(recall)}.txt"
    if not path.is_file():
        raise FileNotFoundError(f"Missing recall CpG list: {path}")
    cdf = pl.read_csv(path, separator="\t", columns=["chr", "start"]).with_columns(
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
        name = chrom[0] if isinstance(chrom, tuple) else chrom
        out[str(name)] = np.sort(grp.get_column("start").to_numpy())
    return out


def overlap_mask(starts: np.ndarray, ends: np.ndarray, cpg_starts: np.ndarray) -> np.ndarray:
    """Boolean mask: read [start, end] contains at least one CpG start."""
    if cpg_starts.size == 0:
        return np.zeros(starts.shape[0], dtype=bool)
    left = np.searchsorted(cpg_starts, starts, side="left")
    right = np.searchsorted(cpg_starts, ends + 1, side="left")
    return right > left


def combo_percentages(
    df: pl.DataFrame,
    threshold: float,
    cpg_positions: Dict[str, np.ndarray],
    mtcount: float,
) -> np.ndarray:
    """Per-chr read-count percentage for one (threshold, recall) combo."""
    filt = df.filter(
        (pl.col("prob_class_1") >= threshold) & (pl.col("mTcount") >= mtcount)
    )
    counts = np.zeros(len(CHR_LIST), dtype=np.int64)
    if filt.height == 0:
        return np.full(len(CHR_LIST), np.nan)

    for chrom, grp in filt.group_by("chr"):
        name = str(chrom[0] if isinstance(chrom, tuple) else chrom)
        if name not in CHR_INDEX:
            continue
        starts = grp.get_column("start").to_numpy()
        ends = grp.get_column("end").to_numpy()
        cpg = cpg_positions.get(name, np.empty(0, dtype=np.int64))
        mask = overlap_mask(starts, ends, cpg)
        counts[CHR_INDEX[name]] = int(mask.sum())

    sum_auto = counts.sum()
    if sum_auto <= 0:
        return np.full(len(CHR_LIST), np.nan)
    return counts.astype(np.float64) / float(sum_auto)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--deconv-res", required=True, type=click.Path(exists=True),
              help="Deconv result (parquet/tsv) with chr,start,end,text,prob_class_1,mTcount.")
@click.option("--best-combo-zscore", required=True, type=click.Path(exists=True),
              help="best_combo_zscore.csv with columns chr,threshold,recall.")
@click.option("--reference-matrix", required=True, type=click.Path(exists=True),
              help="best_reference_matrix.tsv with per-chr percentage_mean/std.")
@click.option("--cpg-recall-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory with 220k_cpg_recall_{recall}.txt files.")
@click.option("--output-prefix", required=True, type=str,
              help="Output prefix; writes {prefix}_zscore.tsv.")
@click.option("--mtcount", type=float, default=1.0, show_default=True,
              help="Minimum mTcount per read (matches workflow/zscore default).")
def main(
    deconv_res: str,
    best_combo_zscore: str,
    reference_matrix: str,
    cpg_recall_dir: str,
    output_prefix: str,
    mtcount: float,
) -> None:
    """Compute flexible per-chromosome read-count zscore for a single sample."""
    try:
        console.rule("[bold blue]Flexible zscore")
        combo_df = load_best_combo(Path(best_combo_zscore))
        ref_mat = load_reference_matrix(Path(reference_matrix))
        recall_dir = Path(cpg_recall_dir)

        df = read_deconv(Path(deconv_res), mtcount)
        console.print(f"  deconv reads (deduped) : {df.height:,}")
        n_combos = combo_df[["threshold", "recall"]].drop_duplicates().shape[0]
        console.print(f"  chromosomes / combos   : {len(combo_df)} / {n_combos}")

        combo_to_chrs: Dict[Tuple[float, float], List[str]] = {}
        for row in combo_df.itertuples(index=False):
            combo_to_chrs.setdefault((row.threshold, row.recall), []).append(row.chr)

        recall_cache: Dict[str, Dict[str, np.ndarray]] = {}
        percentage = np.full(len(CHR_LIST), np.nan)

        for (thr, rec), chrs in combo_to_chrs.items():
            rec_key = _fmt(rec)
            if rec_key not in recall_cache:
                recall_cache[rec_key] = read_recall_positions(recall_dir, rec)
            combo_pct = combo_percentages(df, thr, recall_cache[rec_key], mtcount)
            for chrom in chrs:
                percentage[CHR_INDEX[chrom]] = combo_pct[CHR_INDEX[chrom]]

        ref = ref_mat.reindex(CHR_LIST)
        pct_mean = ref["percentage_mean"].to_numpy(dtype=np.float64)
        pct_std = ref["percentage_std"].to_numpy(dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            std_safe = np.where(pct_std > 0, pct_std, np.nan)
            zscore = (percentage - pct_mean) / std_safe
            zscore = np.where(np.isfinite(zscore), zscore, 0.0)

        import pandas as pd

        out_df = pd.DataFrame({"chr": CHR_LIST, "zscore": zscore})
        out_path = f"{output_prefix}_zscore.tsv"
        out_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
        console.print(f"[green]OK[/green] Wrote {out_path}")
        console.rule("[bold green]Done")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
