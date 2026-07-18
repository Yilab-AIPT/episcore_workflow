#!/usr/bin/env python3
"""Split deconvolution reads into per-threshold target sets in a single pass.

In the AIPT ref-40 flexible-combo setting each chromosome may pick a different
probability threshold for episcore. Rather than re-reading the (very large)
deconvolution table once per threshold, this script reads it exactly once and
writes one target read-name list per requested threshold. Because the target
sets are nested (``prob_class_1 >= t`` shrinks as ``t`` grows) the caller can
turn each list into a ``{sample}__thr_{t}_target.bam`` with a single
``samtools view -N`` per threshold.
"""

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
    """Return True when ``file_path`` starts with the gzip magic bytes."""
    with open(file_path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def detect_file_format(input_path: Path) -> str:
    """Detect whether the input file is txt (TSV) or parquet."""
    suffixes = [s.lower() for s in input_path.suffixes]
    if suffixes and suffixes[-1] == ".gz":
        base_suffix = suffixes[-2] if len(suffixes) >= 2 else ""
    else:
        base_suffix = suffixes[-1] if suffixes else ""

    if base_suffix in [".txt", ".tsv", ".csv"]:
        return "txt"
    if base_suffix in [".parquet", ".pq"]:
        return "parquet"

    try:
        with open(input_path, "rb") as handle:
            if handle.read(4) == b"PAR1":
                return "parquet"
    except Exception:
        pass
    return "txt"


def scan_input(input_path: Path, file_format: str, required_cols: List[str]) -> pl.LazyFrame:
    """Lazily read the deconv file, keeping only ``required_cols``."""
    if file_format == "txt":
        if is_gzipped(input_path):
            with gzip.open(input_path, "rb") as handle:
                raw = handle.read()
            lf = pl.read_csv(io.BytesIO(raw), separator="\t").lazy()
            del raw
        else:
            lf = pl.scan_csv(input_path, separator="\t")
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


def parse_thresholds(thresholds: str) -> List[float]:
    """Parse and sort a comma-separated list of unique probability thresholds."""
    values = sorted({float(t.strip()) for t in thresholds.split(",") if t.strip()})
    if not values:
        raise ValueError("No thresholds provided")
    return values


def write_read_names(output_path: Path, names: pl.Series) -> int:
    """Write sorted, unique read names one per line; return the count."""
    sorted_names = names.unique().sort()
    with open(output_path, "w") as handle:
        for name in sorted_names:
            handle.write(f"{name}\n")
    return len(sorted_names)


@click.command()
@click.option(
    "--input",
    "input_file",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Deconvolution result file (txt/tsv or parquet) with columns name, prob_class_1.",
)
@click.option(
    "--thresholds",
    type=str,
    required=True,
    help="Comma-separated probability thresholds, e.g. '0.1,0.33,0.5,0.67,0.9'.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path.cwd(),
    help="Output directory for the per-threshold read-name files.",
)
@click.option(
    "--background-threshold",
    type=float,
    default=None,
    help=(
        "Also write background_thr_{t}.txt for this threshold "
        "(reads with prob_class_1 < t). Used for raw_total depth filtering."
    ),
)
def main(
    input_file: Path,
    thresholds: str,
    output_dir: Path,
    background_threshold: float | None,
) -> None:
    """Write one ``target_thr_{t}.txt`` per threshold (reads with prob_class_1 >= t).

    The deconv table is read a single time; thresholds are applied in memory.
    When ``--background-threshold`` is set, also emit the complementary
    background read-name list for that threshold.
    """
    try:
        threshold_values = parse_thresholds(thresholds)
        if background_threshold is not None:
            bg_t = float(background_threshold)
            if bg_t not in threshold_values:
                threshold_values = sorted(set(threshold_values) | {bg_t})
        output_dir.mkdir(parents=True, exist_ok=True)

        console.print("\n[bold cyan]Split reads by thresholds[/bold cyan]")
        console.print(f"Input file : {input_file}")
        console.print(f"Thresholds : {threshold_values}")
        if background_threshold is not None:
            console.print(f"Background : {float(background_threshold):g}")
        console.print(f"Output dir : {output_dir}\n")

        file_format = detect_file_format(input_file)
        console.print(f"[green]Detected file format: {file_format}[/green]")

        lf = scan_input(input_file, file_format, ["name", "prob_class_1"])
        df = (
            lf.with_columns(
                pl.col("name").cast(pl.Utf8),
                pl.col("prob_class_1").cast(pl.Float64),
            )
            .drop_nulls()
            .collect()
        )
        console.print(f"  Loaded {len(df):,} rows after dropping nulls")

        for t in threshold_values:
            target_names = df.filter(pl.col("prob_class_1") >= t).get_column("name")
            out_path = output_dir / f"target_thr_{t:g}.txt"
            count = write_read_names(out_path, target_names)
            console.print(
                f"  [green]\u2713[/green] {out_path.name}: {count:,} target reads "
                f"(prob_class_1 >= {t:g})"
            )

        if background_threshold is not None:
            bg_t = float(background_threshold)
            bg_names = df.filter(pl.col("prob_class_1") < bg_t).get_column("name")
            bg_path = output_dir / f"background_thr_{bg_t:g}.txt"
            bg_count = write_read_names(bg_path, bg_names)
            console.print(
                f"  [green]\u2713[/green] {bg_path.name}: {bg_count:,} background reads "
                f"(prob_class_1 < {bg_t:g})"
            )

        console.print("\n[bold green]\u2713 Done[/bold green]\n")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)


if __name__ == "__main__":
    main()
