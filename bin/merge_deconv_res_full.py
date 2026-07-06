#!/usr/bin/env python3
"""Merge multiple deconvolution-result files for a single sample (AIPT ref-40).

Unlike ``merge_deconv_res.py`` (which keeps only ``name``/``prob_class_1``/
``insert_size`` for read-name based BAM splitting), the AIPT ref-40 pipeline
also computes a read-count z-score directly from the deconv table, which needs
the genomic coordinates and ``mTcount``. This script therefore concatenates the
input files while preserving every column required downstream:

    name, chr, start, end, text, prob_class_1, mTcount, insert_size

The episcore read-splitting (``split_reads_by_thresholds.py``) deduplicates on
``name`` and the z-score (``calc_zscore_flexible.py``) deduplicates on
``chr/start/end/text``, so this step performs a plain vertical concatenation and
leaves deduplication to each consumer. Input files may be txt/tsv (optionally
gzipped) or parquet.
"""

import gc
import gzip
import io
import sys
from pathlib import Path
from typing import List

import click
import polars as pl
from rich.console import Console

console = Console()

REQUIRED_COLUMNS = [
    "name",
    "chr",
    "start",
    "end",
    "text",
    "prob_class_1",
    "mTcount",
    "insert_size",
]

CAST_MAP = {
    "name": pl.String,
    "chr": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "text": pl.String,
    "prob_class_1": pl.Float64,
    "mTcount": pl.Float64,
    "insert_size": pl.Int64,
}


def is_gzipped(file_path: Path) -> bool:
    """Check if a file is gzip-compressed by reading its magic bytes."""
    with open(file_path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def detect_file_format(file_path: Path) -> str:
    """Detect whether the input file is txt (TSV) or parquet format."""
    suffixes = [s.lower() for s in file_path.suffixes]
    if suffixes and suffixes[-1] == ".gz":
        base_suffix = suffixes[-2] if len(suffixes) >= 2 else ""
    else:
        base_suffix = suffixes[-1] if suffixes else ""

    if base_suffix in (".txt", ".tsv", ".csv"):
        return "txt"
    if base_suffix in (".parquet", ".pq"):
        return "parquet"

    try:
        with open(file_path, "rb") as f:
            if f.read(4) == b"PAR1":
                return "parquet"
    except Exception:
        pass

    return "txt"


def read_file_lazy(file_path: Path) -> pl.LazyFrame:
    """Read one deconv file into a LazyFrame with the required columns cast."""
    file_format = detect_file_format(file_path)

    if file_format == "parquet":
        lazy_df = pl.scan_parquet(file_path)
    elif is_gzipped(file_path):
        with gzip.open(file_path, "rb") as f:
            raw = f.read()
        lazy_df = pl.read_csv(
            io.BytesIO(raw),
            separator="\t",
            has_header=True,
            null_values=["NA", "na", "N/A", ""],
        ).lazy()
        del raw
    else:
        lazy_df = pl.scan_csv(
            file_path,
            separator="\t",
            has_header=True,
            null_values=["NA", "na", "N/A", ""],
            low_memory=True,
        )

    schema_columns = lazy_df.collect_schema().names()
    missing_columns = [c for c in REQUIRED_COLUMNS if c not in schema_columns]
    if missing_columns:
        raise ValueError(
            f"File {file_path} is missing required columns: {', '.join(missing_columns)}"
        )

    lazy_df = lazy_df.select(REQUIRED_COLUMNS).with_columns(
        [pl.col(col).cast(dtype) for col, dtype in CAST_MAP.items()]
    )
    return lazy_df


def merge_files(input_files: List[Path], output_path: Path) -> None:
    """Concatenate the deconv files preserving all required columns.

    The deconv files can be very large (tens of GB of TSV), so the merged result
    is streamed straight to a parquet file via ``sink_parquet`` to avoid
    materialising the whole table in memory. A standard collect+write is used as
    a fallback if streaming is unavailable.
    """
    lazy_frames: List[pl.LazyFrame] = []
    for file_path in input_files:
        console.print(f"[cyan]Preparing[/cyan] {file_path.name}")
        lazy_frames.append(read_file_lazy(file_path))

    if not lazy_frames:
        raise ValueError("No valid input files could be read")

    merged_lazy = pl.concat(lazy_frames, how="vertical")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        merged_lazy.sink_parquet(output_path)
    except Exception as exc:
        console.print(
            f"[yellow]sink_parquet failed ({exc}); falling back to in-memory collect.[/yellow]"
        )
        output_df = merged_lazy.collect(engine="streaming")
        if output_df.height == 0:
            raise ValueError("No data remains after merging. Check your input files.")
        output_df.write_parquet(output_path, use_pyarrow=True)
        del output_df

    gc.collect()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--inputs",
    required=True,
    type=str,
    help="Space-separated list of input deconv files (txt/tsv/gz or parquet).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output parquet path for the merged deconv table.",
)
def main(inputs: str, output: str) -> None:
    """Merge multiple deconv files for one sample, preserving all needed columns."""
    console.rule("[bold blue]Deconv merge (full columns)")
    try:
        input_file_paths = [Path(f.strip()) for f in inputs.split() if f.strip()]
        if not input_file_paths:
            raise ValueError("No input files provided")

        for fp in input_file_paths:
            if not fp.is_file():
                raise FileNotFoundError(f"Input file not found: {fp}")

        console.print(f"  input files : {len(input_file_paths)}")
        console.print(f"  output      : {output}")

        merge_files(input_file_paths, Path(output))
        console.print(f"[green]OK[/green] Wrote {output}")
        console.rule("[bold green]Done")
    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
