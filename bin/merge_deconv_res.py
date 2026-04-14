#!/usr/bin/env python3
"""
Merge multiple files containing genomic deconvolution results.

This script merges multiple files (txt/tsv or parquet format) with columns 'name',
'prob_class_1', 'mTcount', and 'insert_size', optionally filters by mTcount threshold,
and outputs a merged parquet file with 'name', 'prob_class_1', and 'insert_size'.
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
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

console = Console()


def validate_input_files(input_files: List[Path]) -> None:
    """
    Validate that all input files exist and are readable.

    Args:
        input_files: List of Path objects pointing to input files.

    Raises:
        FileNotFoundError: If any input file does not exist.
        PermissionError: If any input file is not readable.
    """
    for file_path in input_files:
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")
        if not file_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")


def is_gzipped(file_path: Path) -> bool:
    """Check if a file is gzip-compressed by reading its magic bytes."""
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def detect_file_format(file_path: Path) -> str:
    """Detect whether the input file is txt (TSV) or parquet format.

    Checks file extension first, then falls back to magic byte inspection.
    Handles gzip-compressed files (e.g. .txt.gz).

    Returns:
        'txt' or 'parquet'.
    """
    suffixes = [s.lower() for s in file_path.suffixes]

    if suffixes and suffixes[-1] == '.gz':
        base_suffix = suffixes[-2] if len(suffixes) >= 2 else ''
    else:
        base_suffix = suffixes[-1] if suffixes else ''

    if base_suffix in ('.txt', '.tsv', '.csv'):
        return 'txt'
    if base_suffix in ('.parquet', '.pq'):
        return 'parquet'

    # Fall back to magic bytes
    try:
        with open(file_path, 'rb') as f:
            magic = f.read(4)
            if magic == b'PAR1':
                return 'parquet'
    except Exception:
        pass

    if is_gzipped(file_path):
        return 'txt'

    console.print(
        f"[yellow]Warning: Cannot determine format of '{file_path.name}', "
        f"assuming txt/tsv[/yellow]"
    )
    return 'txt'


def read_file_lazy(file_path: Path, required_columns: List[str], ncpgs: int = 0) -> pl.LazyFrame:
    """
    Read a single file (txt/tsv or parquet) into a Polars LazyFrame.

    Args:
        file_path: Path to the input file.
        required_columns: List of column names that must be present.
        ncpgs: Minimum mTcount threshold for filtering (applied lazily).

    Returns:
        A Polars LazyFrame with filters and column selection applied.

    Raises:
        ValueError: If required columns are missing from the file.
    """
    try:
        file_format = detect_file_format(file_path)

        if file_format == 'parquet':
            lazy_df = pl.scan_parquet(file_path)
        elif is_gzipped(file_path):
            with gzip.open(file_path, 'rb') as f:
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
        missing_columns = [col for col in required_columns if col not in schema_columns]
        if missing_columns:
            raise ValueError(
                f"File {file_path} is missing required columns: {', '.join(missing_columns)}"
            )

        if ncpgs > 0:
            lazy_df = lazy_df.filter(pl.col("mTcount") >= ncpgs)

        lazy_df = lazy_df.select(required_columns)

        lazy_df = lazy_df.with_columns([
            pl.col("name").cast(pl.String),
            pl.col("prob_class_1").cast(pl.Float64),
            pl.col("mTcount").cast(pl.Float64),
            pl.col("insert_size").cast(pl.Int64),
        ])

        return lazy_df

    except Exception as e:
        console.print(f"[bold red]Error reading file {file_path}:[/bold red] {str(e)}")
        raise


def merge_files(
    input_files: List[Path],
    ncpgs: int = 0
) -> pl.DataFrame:
    """
    Merge multiple files and apply filtering using lazy evaluation.

    Args:
        input_files: List of Path objects pointing to input files.
        ncpgs: Minimum mTcount threshold for filtering (default: 0, no filtering).

    Returns:
        A merged and filtered Polars DataFrame with 'name', 'prob_class_1', and 'insert_size' columns.

    Raises:
        ValueError: If no valid data remains after filtering.
    """
    required_columns = ["name", "prob_class_1", "mTcount", "insert_size"]
    lazy_frames = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console
    ) as progress:
        task = progress.add_task(
            "[cyan]Preparing input files...",
            total=len(input_files)
        )

        for file_path in input_files:
            try:
                lazy_df = read_file_lazy(file_path, required_columns, ncpgs)
                lazy_frames.append(lazy_df)
                progress.update(
                    task,
                    advance=1,
                    description=f"[cyan]Preparing input files... ({file_path.name})"
                )
            except Exception as e:
                console.print(f"[bold yellow]Warning:[/bold yellow] Skipping {file_path}: {str(e)}")
                progress.update(task, advance=1)

    if not lazy_frames:
        raise ValueError("No valid input files could be read")

    console.print(f"[green]✓[/green] Successfully prepared {len(lazy_frames)} file(s)")

    console.print("[cyan]Merging and processing data (this may take a while for large files)...[/cyan]")

    merged_lazy = pl.concat(lazy_frames, how="vertical")

    output_lazy = merged_lazy.select(["name", "prob_class_1", "insert_size"])

    try:
        output_df = output_lazy.collect(engine="streaming")
    except Exception:
        console.print("[yellow]⚠[/yellow] Streaming mode failed, falling back to standard collection...")
        output_df = output_lazy.collect()

    gc.collect()

    output_rows = len(output_df)
    console.print(f"[green]✓[/green] Processed data contains {output_rows:,} rows")

    if ncpgs > 0:
        console.print(f"[green]✓[/green] Filter applied: mTcount >= {ncpgs}")
    else:
        console.print("[yellow]ℹ[/yellow] No filtering applied (ncpgs = 0)")

    if output_rows == 0:
        raise ValueError(
            f"No data remains after processing. "
            f"Check your input files and filter threshold (ncpgs={ncpgs})."
        )

    return output_df


def write_output(df: pl.DataFrame, output_path: Path) -> None:
    """
    Write the merged DataFrame to an output file.

    Supports parquet (default), gzipped TSV, and plain TSV based on file extension.

    Args:
        df: Polars DataFrame to write.
        output_path: Path where the output file will be written.
    """
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print("[cyan]Writing output file...[/cyan]")

        suffix = output_path.suffix.lower()
        if suffix in ('.parquet', '.pq'):
            df.write_parquet(output_path, use_pyarrow=True)
        elif str(output_path).endswith('.gz'):
            with gzip.open(output_path, 'wb') as f:
                df.write_csv(f, separator="\t", include_header=True)
        else:
            df.write_csv(output_path, separator="\t", include_header=True)

        console.print(f"[green]✓[/green] Output written to: {output_path}")
        console.print(f"[green]✓[/green] Output contains {len(df):,} rows")

        gc.collect()

    except Exception as e:
        console.print(f"[bold red]Error writing output file:[/bold red] {str(e)}")
        raise


@click.command()
@click.option(
    "--inputs",
    required=True,
    type=str,
    help="Space-separated list of input files to merge (txt/tsv or parquet)."
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Output file path for the merged data (use .parquet extension for parquet output)."
)
@click.option(
    "--ncpgs",
    default=0,
    type=int,
    help="Minimum mTcount threshold for filtering (default: 0, no filtering)."
)
def main(inputs: str, output: str, ncpgs: int) -> None:
    """
    Merge multiple files containing genomic deconvolution results.

    Reads multiple files (txt/tsv or parquet) with columns 'name', 'prob_class_1',
    'mTcount', and 'insert_size', merges them, applies optional mTcount filtering,
    and outputs a file with 'name', 'prob_class_1', and 'insert_size' columns.

    Examples:
        merge_deconv_res.py --inputs "file1.parquet file2.parquet" --output merged.parquet

        merge_deconv_res.py --inputs "file1.txt file2.txt" --output merged.parquet --ncpgs 5
    """
    console.print("[bold cyan]Deconvolution Results Merger[/bold cyan]")
    console.print("=" * 60)

    try:
        input_file_paths = [Path(f.strip()) for f in inputs.split() if f.strip()]

        if not input_file_paths:
            raise ValueError("No input files provided")

        console.print(f"[cyan]Input files:[/cyan] {len(input_file_paths)}")
        console.print(f"[cyan]Output file:[/cyan] {output}")
        console.print(f"[cyan]mTcount threshold:[/cyan] {ncpgs}")
        console.print("=" * 60)

        console.print("[cyan]Validating input files...[/cyan]")
        validate_input_files(input_file_paths)
        console.print(f"[green]✓[/green] All {len(input_file_paths)} input file(s) validated")

        merged_df = merge_files(input_file_paths, ncpgs)

        output_path = Path(output)
        write_output(merged_df, output_path)

        del merged_df
        gc.collect()

        console.print("=" * 60)
        console.print("[bold green]✓ Process completed successfully![/bold green]")

    except Exception as e:
        console.print("=" * 60)
        console.print(f"[bold red]✗ Error:[/bold red] {str(e)}")
        console.print("[bold red]Process failed![/bold red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
