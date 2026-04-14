#!/usr/bin/env python3
"""Filter deconvolution results and split reads by probability threshold.

This script processes deconvolution result files (txt/parquet format) and splits
read names into target and background sets based on a probability threshold.
The output files can be used with samtools to split BAM files accordingly.
"""

import gzip
import io
import sys
from pathlib import Path

import click
import polars as pl
from rich.console import Console

console = Console()


def is_gzipped(file_path: Path) -> bool:
    """Check if a file is gzip-compressed by reading its magic bytes."""
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def detect_file_format(input_path: Path) -> str:
    """Detect whether the input file is txt (TSV) or parquet format.

    Handles gzip-compressed files (e.g. .txt.gz) by inspecting the
    second-to-last suffix.

    Args:
        input_path: Path to the input file.

    Returns:
        File format as string: 'txt' or 'parquet'.

    Raises:
        ValueError: If file format cannot be determined or is unsupported.
    """
    suffixes = [s.lower() for s in input_path.suffixes]

    # Strip trailing .gz to examine the real extension
    if suffixes and suffixes[-1] == '.gz':
        base_suffix = suffixes[-2] if len(suffixes) >= 2 else ''
    else:
        base_suffix = suffixes[-1] if suffixes else ''

    if base_suffix in ['.txt', '.tsv', '.csv']:
        return 'txt'
    elif base_suffix in ['.parquet', '.pq']:
        return 'parquet'
    else:
        try:
            with open(input_path, 'rb') as f:
                magic = f.read(4)
                if magic == b'PAR1':
                    return 'parquet'
        except Exception:
            pass

        # If gzip-compressed with no recognisable inner extension, assume txt
        if is_gzipped(input_path):
            return 'txt'

        console.print(
            f"[yellow]Warning: Cannot determine file format from extension "
            f"'{input_path.suffix}', assuming txt/tsv format[/yellow]"
        )
        return 'txt'


def scan_input(
    input_path: Path,
    file_format: str,
    required_cols: list[str],
) -> pl.LazyFrame:
    """Create a lazy frame from the input file, selecting only required columns.

    Args:
        input_path: Path to the input file.
        file_format: 'txt' or 'parquet'.
        required_cols: List of columns to select.

    Returns:
        LazyFrame with only the required columns.

    Raises:
        ValueError: If any required column is missing.
    """
    if file_format == 'txt':
        if is_gzipped(input_path):
            with gzip.open(input_path, 'rb') as f:
                raw = f.read()
            lf = pl.read_csv(io.BytesIO(raw), separator='\t').lazy()
            del raw
        else:
            lf = pl.scan_csv(input_path, separator='\t')
    else:
        lf = pl.scan_parquet(input_path)

    available = lf.collect_schema().names()
    missing = [c for c in required_cols if c not in available]
    if missing:
        raise ValueError(
            f"Required column(s) {missing} not found in input file. "
            f"Available columns: {available}"
        )

    return lf.select(required_cols)


def write_read_names(output_path: Path, names: pl.Series) -> int:
    """Write read names to a file, one per line, sorted.

    Args:
        output_path: Path to the output file.
        names: Series of read names to write.

    Returns:
        Number of names written.
    """
    sorted_names = names.sort()
    with open(output_path, 'w') as f:
        for name in sorted_names:
            f.write(f"{name}\n")
    return len(sorted_names)


@click.command()
@click.option(
    '--input',
    'input_file',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Input file path (txt/tsv or parquet format) with columns: name, prob_class_1',
)
@click.option(
    '--threshold',
    type=float,
    default=None,
    help='Probability threshold to filter reads (prob_class_1 > threshold -> target)',
)
@click.option(
    '--output-dir',
    type=click.Path(path_type=Path),
    default=Path.cwd(),
    help='Output directory for result files (default: current directory)',
)
@click.option(
    '--size',
    type=int,
    default=None,
    help=(
        'Insert size threshold. When set, target reads must also have '
        'insert_size < size (requires insert_size column in input).'
    ),
)
def main(
    input_file: Path,
    threshold: float | None,
    output_dir: Path,
    size: int | None,
):
    """Filter deconvolution results and split reads into target/background.

    Reads a deconvolution result file (txt/tsv or parquet) and splits read names
    based on the provided filtering criteria. At least one of --threshold or
    --size must be specified.

    \b
    Filtering modes:
    - --threshold only: target reads have prob_class_1 > threshold
    - --size only: target reads have insert_size < size
    - Both: target reads satisfy both conditions

    \b
    Output files:
    - target_reads.txt: Read names passing the filter(s)
    - background_reads.txt: Read names failing the filter(s)
    - classified_reads.txt: All read names

    The output files contain unique, sorted read names (one per line) suitable
    for use with 'samtools view -N' to split BAM files.

    Examples:
        filter_deconv_res.py --input deconv.txt --threshold 0.5

        filter_deconv_res.py --input deconv.txt --size 200

        filter_deconv_res.py --input deconv.txt --threshold 0.5 --size 200
    """
    try:
        if threshold is None and size is None:
            console.print(
                "[bold red]Error:[/bold red] At least one of --threshold or "
                "--size must be specified."
            )
            sys.exit(1)

        console.print("\n[bold cyan]Deconvolution Results Filter[/bold cyan]")
        console.print(f"Input file: {input_file}")
        if threshold is not None:
            console.print(f"Threshold: {threshold}")
        if size is not None:
            console.print(f"Insert size filter: < {size}")
        console.print(f"Output directory: {output_dir}\n")

        output_dir.mkdir(parents=True, exist_ok=True)

        file_format = detect_file_format(input_file)
        console.print(f"[green]Detected file format: {file_format}[/green]")

        required_cols = ['name']
        cast_exprs: list[pl.Expr] = [pl.col('name').cast(pl.Utf8)]
        if threshold is not None:
            required_cols.append('prob_class_1')
            cast_exprs.append(pl.col('prob_class_1').cast(pl.Float64))
        if size is not None:
            required_cols.append('insert_size')
            cast_exprs.append(pl.col('insert_size').cast(pl.Int64))

        console.print("[cyan]Reading and filtering data...[/cyan]")

        try:
            lf = scan_input(input_file, file_format, required_cols)
            df = lf.with_columns(cast_exprs).drop_nulls().collect()

            console.print(f"  Loaded {len(df):,} rows after dropping nulls")

            conditions: list[pl.Expr] = []
            if threshold is not None:
                conditions.append(pl.col('prob_class_1') >= threshold)
            if size is not None:
                conditions.append(pl.col('insert_size') <= size)

            target_filter = conditions[0]
            for cond in conditions[1:]:
                target_filter = target_filter & cond

            target_names = df.filter(target_filter).get_column('name').unique()
            background_names = df.filter(~target_filter).get_column('name').unique()
            classified_names = df.get_column('name').unique()

        except ValueError as e:
            console.print(f"[bold red]Error:[/bold red] {e}")
            sys.exit(1)
        except Exception as e:
            console.print(
                f"[bold red]Unexpected error during processing:[/bold red] {e}"
            )
            sys.exit(1)

        console.print("\n[cyan]Writing output files...[/cyan]")

        output_data = {
            'target_reads.txt': target_names,
            'background_reads.txt': background_names,
            'classified_reads.txt': classified_names,
        }

        for filename, names in output_data.items():
            output_path = output_dir / filename
            count = write_read_names(output_path, names)
            console.print(f"  [green]✓[/green] {filename}: {count:,} unique reads")

        n_classified = len(classified_names)
        n_target = len(target_names)
        n_background = len(background_names)

        filter_parts = []
        if threshold is not None:
            filter_parts.append(f"prob >= {threshold}")
        if size is not None:
            filter_parts.append(f"insert_size <= {size}")
        filter_label = ", ".join(filter_parts)

        console.print("\n[bold green]Summary:[/bold green]")
        console.print(f"  Total rows processed: {len(df):,}")
        console.print(f"  Unique classified reads: {n_classified:,}")
        console.print(
            f"  Unique target reads ({filter_label}): "
            f"{n_target:,} ({n_target / n_classified * 100:.2f}%)"
        )
        console.print(
            f"  Unique background reads: {n_background:,} "
            f"({n_background / n_classified * 100:.2f}%)"
        )

        console.print(
            "\n[bold green]✓ Processing completed successfully![/bold green]\n"
        )

    except KeyboardInterrupt:
        console.print("\n[yellow]Operation cancelled by user[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[bold red]Fatal error:[/bold red] {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
