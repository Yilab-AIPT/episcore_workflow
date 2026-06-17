#!/usr/bin/env python3
"""
Replace original read probabilities with perturbed read probabilities.

This script takes a merged *original* deconvolution result (columns 'name',
'prob_class_1', 'insert_size') and a *perturbed* deconvolution result (reads
whose methylation status was randomly perturbed and then re-inferred; columns
'name', 'prob_class_1', ...). For every read present in the perturbed file the
original 'prob_class_1' is overwritten by the perturbed value; reads absent from
the perturbed file keep their original probability. The output is a deconvolution
result (columns 'name', 'prob_class_1', 'insert_size') ready for split_bam.
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

    try:
        with open(file_path, 'rb') as f:
            if f.read(4) == b'PAR1':
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


def scan_file(file_path: Path, required_columns: List[str]) -> pl.LazyFrame:
    """Read a single file (txt/tsv or parquet) into a Polars LazyFrame.

    Args:
        file_path: Path to the input file.
        required_columns: Columns that must be present; only these are selected.

    Returns:
        A Polars LazyFrame restricted to ``required_columns``.

    Raises:
        ValueError: If required columns are missing from the file.
    """
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
    missing_columns = [c for c in required_columns if c not in schema_columns]
    if missing_columns:
        raise ValueError(
            f"File {file_path} is missing required columns: {', '.join(missing_columns)}"
        )

    return lazy_df.select(required_columns)


def collect_lazy(lazy_df: pl.LazyFrame) -> pl.DataFrame:
    """Collect a LazyFrame, preferring the streaming engine for large inputs."""
    try:
        return lazy_df.collect(engine="streaming")
    except Exception:
        console.print(
            "[yellow]⚠[/yellow] Streaming mode failed, falling back to standard collection..."
        )
        return lazy_df.collect()


def replace_probabilities(
    original_path: Path,
    perturbed_path: Path,
) -> pl.DataFrame:
    """Overwrite original read probabilities with perturbed ones.

    Args:
        original_path: Merged original deconvolution result (name, prob_class_1, insert_size).
        perturbed_path: Perturbed deconvolution result (name, prob_class_1, ...).

    Returns:
        A Polars DataFrame with columns 'name', 'prob_class_1', 'insert_size' in
        which 'prob_class_1' is taken from the perturbed file where the read name
        matches, and from the original file otherwise.

    Raises:
        ValueError: If no rows remain after the replacement.
    """
    original_lazy = scan_file(original_path, ["name", "prob_class_1", "insert_size"]).with_columns([
        pl.col("name").cast(pl.String),
        pl.col("prob_class_1").cast(pl.Float64),
        pl.col("insert_size").cast(pl.Int64),
    ])

    # Perturbed reads are re-inferred per read; keep one probability per name.
    perturbed_lazy = (
        scan_file(perturbed_path, ["name", "prob_class_1"])
        .with_columns([
            pl.col("name").cast(pl.String),
            pl.col("prob_class_1").cast(pl.Float64).alias("perturbed_prob_class_1"),
        ])
        .select(["name", "perturbed_prob_class_1"])
        .unique(subset=["name"], keep="first")
    )

    console.print("[cyan]Joining original and perturbed probabilities (left join on read name)...[/cyan]")

    merged_lazy = (
        original_lazy
        .join(perturbed_lazy, on="name", how="left")
        .with_columns(
            pl.coalesce(["perturbed_prob_class_1", "prob_class_1"]).alias("prob_class_1")
        )
        .select(["name", "prob_class_1", "insert_size"])
    )

    output_df = collect_lazy(merged_lazy)
    gc.collect()

    if output_df.height == 0:
        raise ValueError(
            "No rows remain after replacing probabilities. "
            "Check that the original and perturbed files share read names."
        )

    return output_df


def write_output(df: pl.DataFrame, output_path: Path) -> None:
    """Write the deconvolution result, dispatching by output extension."""
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
    console.print(f"[green]✓[/green] Output contains {df.height:,} rows")
    gc.collect()


@click.command()
@click.option(
    "--original",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Merged original deconvolution result (txt/tsv or parquet) with columns "
         "'name', 'prob_class_1', 'insert_size'.",
)
@click.option(
    "--perturbed",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Perturbed deconvolution result (txt/tsv or parquet) with columns "
         "'name' and 'prob_class_1'.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output file path (use .parquet extension for parquet output).",
)
def main(original: Path, perturbed: Path, output: Path) -> None:
    """Replace original read probabilities with perturbed read probabilities.

    Examples:
        replace_deconv_prob.py --original merged.parquet \\
            --perturbed perturbed.parquet --output replaced.parquet
    """
    console.print("[bold cyan]Deconvolution Probability Replacer[/bold cyan]")
    console.print("=" * 60)
    console.print(f"[cyan]Original:[/cyan]  {original}")
    console.print(f"[cyan]Perturbed:[/cyan] {perturbed}")
    console.print(f"[cyan]Output:[/cyan]    {output}")
    console.print("=" * 60)

    try:
        result_df = replace_probabilities(original, perturbed)
        write_output(result_df, output)

        del result_df
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
