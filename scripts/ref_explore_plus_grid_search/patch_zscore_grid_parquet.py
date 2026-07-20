#!/usr/bin/env python3
"""Append newly generated zscore CSVs into zscore_grid_search.parquet."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import List

import click
import pandas as pd
from rich.console import Console

from grid_coverage import fmt_float

console = Console()

CSV_USECOLS = ["sample", "chr", "percentage"]
CSV_DTYPES = {"sample": "string", "chr": "string", "percentage": "float64"}
OUTPUT_COLS = ["sample", "chr", "threshold", "recall", "percentage"]

DEFAULT_ZSCORE_ROOT = (
    "/lustre1/cqyi/AIPT_2.0/results/zscore_output/20260513-grid_search/zscore_results"
)
DEFAULT_PARQUET = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/"
    "20260621-ref_40_rebuild_consider_lib_ng/zscore_grid_search.parquet"
)


def _zscore_csv_path(
    zscore_root: Path, sample: str, threshold: float, recall: float
) -> Path:
    thr_s = fmt_float(threshold)
    rec_s = fmt_float(recall)
    combo_dir = zscore_root / f"recall.{rec_s}_cutoff.{thr_s}"
    return (
        combo_dir
        / f"{sample}.{thr_s}.1.0.220k_cpg_recall_{rec_s}.NoLen.zscore.csv"
    )


def _load_csv(path: Path, threshold: float, recall: float) -> pd.DataFrame:
    frame = pd.read_csv(path, usecols=CSV_USECOLS, dtype=CSV_DTYPES)
    if len(frame) < 22:
        raise click.ClickException(f"{path} has only {len(frame)} rows (need >=22)")
    frame["threshold"] = float(threshold)
    frame["recall"] = float(recall)
    return frame[OUTPUT_COLS]


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--missing-tsv",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="missing_coverage.tsv from check_grid_coverage.py",
)
@click.option(
    "--zscore-root",
    default=DEFAULT_ZSCORE_ROOT,
    show_default=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--parquet-path",
    default=DEFAULT_PARQUET,
    show_default=True,
    type=click.Path(path_type=Path),
)
@click.option("--backup/--no-backup", default=True, show_default=True)
def main(
    missing_tsv: Path,
    zscore_root: Path,
    parquet_path: Path,
    backup: bool,
) -> None:
    """Merge filled zscore CSV rows into the grid-search parquet."""
    missing = pd.read_csv(missing_tsv, sep="\t")
    if "score_type" in missing.columns:
        missing = missing[missing["score_type"].astype(str) == "zscore"].copy()
    if missing.empty:
        console.print("[yellow]No zscore missing rows to patch[/yellow]")
        return

    frames: List[pd.DataFrame] = []
    still_missing: List[str] = []
    for _, row in missing.iterrows():
        sample = str(row["sample"])
        thr = float(row["threshold"])
        rec = float(row["recall"])
        path = _zscore_csv_path(zscore_root, sample, thr, rec)
        if not path.is_file():
            still_missing.append(str(path))
            continue
        frames.append(_load_csv(path, thr, rec))

    if still_missing:
        preview = "\n".join(f"  - {p}" for p in still_missing[:20])
        extra = (
            f"\n  ... and {len(still_missing) - 20} more"
            if len(still_missing) > 20
            else ""
        )
        raise click.ClickException(
            f"{len(still_missing)} zscore CSVs still missing:\n{preview}{extra}"
        )

    new_rows = pd.concat(frames, ignore_index=True)
    console.print(f"  new rows to append : {len(new_rows)}")

    if not parquet_path.is_file():
        raise click.ClickException(f"Parquet not found: {parquet_path}")

    if backup:
        bak = parquet_path.with_suffix(parquet_path.suffix + ".bak_before_patch")
        if not bak.is_file():
            shutil.copy2(parquet_path, bak)
            console.print(f"  backup            : {bak}")

    console.print("[cyan]Loading existing parquet ...[/cyan]")
    old = pd.read_parquet(parquet_path)
    # Drop any pre-existing rows for these keys (should be none)
    keys = missing[["sample", "threshold", "recall"]].drop_duplicates()
    keys["sample"] = keys["sample"].astype(str)
    keys["threshold"] = keys["threshold"].astype(float)
    keys["recall"] = keys["recall"].astype(float)
    old["sample"] = old["sample"].astype(str)
    old["threshold"] = old["threshold"].astype(float)
    old["recall"] = old["recall"].astype(float)
    merged_keys = old.merge(keys, on=["sample", "threshold", "recall"], how="left", indicator=True)
    keep = old.loc[merged_keys["_merge"] == "left_only"].copy()
    n_dropped = len(old) - len(keep)
    out = pd.concat([keep, new_rows], ignore_index=True)

    tmp = parquet_path.with_suffix(".parquet.tmp")
    out.to_parquet(tmp, index=False, compression="snappy")
    tmp.replace(parquet_path)
    console.print(f"[green]Done[/green] Patched {parquet_path}")
    console.print(f"  dropped prior rows : {n_dropped}")
    console.print(f"  final rows         : {len(out):,}")


if __name__ == "__main__":
    main()
