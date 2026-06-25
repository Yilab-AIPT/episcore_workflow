#!/usr/bin/env python3
"""Build long-format zscore grid-search table from per-sample CSV outputs.

Reads ``recall.{recall}_cutoff.{threshold}/{sample}*.csv`` files, writes one
parquet shard per combo directory, then streams shards into a single output
parquet without loading the full table into memory.
"""

from __future__ import annotations

import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

import click
import pandas as pd
import pyarrow.parquet as pq
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

console = Console()

ZSCORE_COMBO_RE = re.compile(r"^recall\.(?P<recall>[\d.]+)_cutoff\.(?P<threshold>[\d.]+)$")
OUTPUT_COLS = ["sample", "chr", "threshold", "recall", "percentage"]
CSV_USECOLS = ["sample", "chr", "percentage"]
CSV_DTYPES = {"sample": "string", "chr": "string", "percentage": "float64"}


def _resolve_threads(threads: int) -> int:
    if threads > 0:
        return threads
    from_env = int(os.environ.get("SLURM_CPUS_PER_TASK", "0") or 0)
    return from_env or 32


def _list_combo_dirs(root: Path) -> List[Path]:
    combos: List[Path] = []
    with os.scandir(root) as entries:
        for entry in entries:
            if entry.is_dir() and ZSCORE_COMBO_RE.match(entry.name):
                combos.append(Path(entry.path))
    combos.sort(key=lambda path: path.name)
    return combos


def _list_csv_paths(combo_dir: Path) -> List[Path]:
    csv_paths: List[Path] = []
    with os.scandir(combo_dir) as entries:
        for entry in entries:
            if entry.is_file(follow_symlinks=False) and entry.name.endswith(".csv"):
                csv_paths.append(Path(entry.path))
    csv_paths.sort()
    return csv_paths


def _process_combo(combo_dir: Path, shard_dir: Path) -> Optional[Path]:
    match = ZSCORE_COMBO_RE.match(combo_dir.name)
    if not match:
        return None

    threshold = float(match.group("threshold"))
    recall = float(match.group("recall"))
    shard_path = shard_dir / f"{combo_dir.name}.parquet"
    if shard_path.is_file():
        return shard_path

    csv_paths = _list_csv_paths(combo_dir)
    if not csv_paths:
        return None

    frames = [
        pd.read_csv(path, usecols=CSV_USECOLS, dtype=CSV_DTYPES)
        for path in csv_paths
    ]
    out = pd.concat(frames, ignore_index=True, copy=False)
    out["threshold"] = threshold
    out["recall"] = recall
    out = out[OUTPUT_COLS]
    out.to_parquet(shard_path, index=False, compression="snappy")
    return shard_path


def _merge_shards(shard_dir: Path, output_path: Path) -> int:
    shard_paths = sorted(shard_dir.glob("*.parquet"))
    if not shard_paths:
        raise FileNotFoundError(f"No parquet shards found under {shard_dir}")

    writer: Optional[pq.ParquetWriter] = None
    total_rows = 0
    try:
        for shard_path in shard_paths:
            table = pq.read_table(shard_path)
            total_rows += table.num_rows
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema, compression="snappy")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    return total_rows


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--zscore-grid-search-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Root directory containing recall.{recall}_cutoff.{threshold}/ subdirs",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    help="Directory for zscore_grid_search.parquet",
)
@click.option(
    "--threads",
    default=0,
    show_default="SLURM_CPUS_PER_TASK or 32",
    type=int,
    help="Worker threads for parallel combo parsing",
)
@click.option(
    "--keep-shards",
    is_flag=True,
    default=False,
    help="Keep intermediate per-combo parquet shards after merge",
)
def main(
    zscore_grid_search_dir: Path,
    output_dir: Path,
    threads: int,
    keep_shards: bool,
) -> None:
    """Parse zscore grid-search CSVs into one long-format parquet table."""
    threads = _resolve_threads(threads)
    output_dir.mkdir(parents=True, exist_ok=True)
    shard_dir = output_dir / "_zscore_grid_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "zscore_grid_search.parquet"

    combo_dirs = _list_combo_dirs(zscore_grid_search_dir)
    if not combo_dirs:
        raise click.ClickException(
            f"No recall.*_cutoff.* combo directories found under {zscore_grid_search_dir}"
        )

    console.rule("[bold blue]Build zscore grid-search table")
    console.print(f"  input root : {zscore_grid_search_dir}")
    console.print(f"  output     : {output_path}")
    console.print(f"  shard dir  : {shard_dir}")
    console.print(f"  combos     : {len(combo_dirs)}")
    console.print(f"  threads    : {threads}")

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        console=console,
    )
    shard_paths: List[Path] = []
    with progress:
        task = progress.add_task("Parsing combo directories", total=len(combo_dirs))
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {
                pool.submit(_process_combo, combo_dir, shard_dir): combo_dir
                for combo_dir in combo_dirs
            }
            for future in as_completed(futures):
                combo_dir = futures[future]
                try:
                    shard_path = future.result()
                except Exception as exc:
                    raise click.ClickException(
                        f"Failed parsing combo directory {combo_dir}: {exc}"
                    ) from exc
                if shard_path is not None:
                    shard_paths.append(shard_path)
                progress.advance(task)

    console.print(f"  shards written : {len(shard_paths)}")
    total_rows = _merge_shards(shard_dir, output_path)
    console.print(f"  merged rows    : {total_rows:,}")
    console.print(f"  wrote          : {output_path}")

    if not keep_shards:
        shutil.rmtree(shard_dir)


if __name__ == "__main__":
    main()
