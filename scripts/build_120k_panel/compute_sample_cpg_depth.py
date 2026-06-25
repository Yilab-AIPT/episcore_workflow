#!/usr/bin/env python3
"""Compute raw read depth at CpG sites from a bisulfite BAM.

For each site in the recall-0.6 CpG list, counts pileup depth (pysam) with
MAPQ >= min_mapq and base quality >= min_bq. Writes one TSV per sample.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import click
import polars as pl
import pysam
from rich.console import Console

from bam_utils import bam_index_path

console = Console()

DEFAULT_SITES = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list"
    "/220k_cpg_recall_0.6.txt"
)


def open_bam(bam_path: Path) -> pysam.AlignmentFile:
    """Open BAM with an adjacent .bai (required for pileup)."""
    real = Path(os.path.realpath(bam_path))
    bai = bam_index_path(bam_path)
    if bai is None:
        raise FileNotFoundError(f"No BAM index (.bai) found for {bam_path}")
    return pysam.AlignmentFile(str(real), "rb", index_filename=str(bai))


def load_sites(sites_path: Path) -> pl.DataFrame:
    return pl.read_csv(
        sites_path,
        separator="\t",
        columns=["chr", "start", "end"],
    )


def compute_depths(
    bam_path: Path,
    sites: pl.DataFrame,
    min_mapq: int,
    min_bq: int,
) -> pl.DataFrame:
    by_chr: dict[str, set[int]] = defaultdict(set)
    for row in sites.iter_rows(named=True):
        by_chr[row["chr"]].add(row["start"])

    depth_map: dict[tuple[str, int, int], int] = {}
    with open_bam(bam_path) as bam:
        ref_names = set(bam.references)
        for chrom, pos_set in by_chr.items():
            if chrom not in ref_names:
                continue
            positions = sorted(pos_set)
            region_start, region_end = positions[0], positions[-1] + 1
            for col in bam.pileup(
                chrom,
                region_start,
                region_end,
                truncate=True,
                stepper="all",
                min_mapping_quality=min_mapq,
                min_base_quality=min_bq,
                max_depth=1_000_000,
            ):
                pos = col.reference_pos
                if pos in pos_set:
                    depth_map[(chrom, pos, pos + 1)] = col.n

    if not depth_map:
        return sites.with_columns(pl.lit(0).alias("raw_depth").cast(pl.Int32))

    depth_df = pl.DataFrame(
        {
            "chr": [k[0] for k in depth_map],
            "start": [k[1] for k in depth_map],
            "end": [k[2] for k in depth_map],
            "raw_depth": list(depth_map.values()),
        }
    )
    return (
        sites.join(depth_df, on=["chr", "start", "end"], how="left")
        .with_columns(pl.col("raw_depth").fill_null(0).cast(pl.Int32))
    )


@click.command()
@click.option("--bam", required=True, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--sample-id",
    required=True,
    help="Sample identifier used in the output filename.",
)
@click.option(
    "--sites",
    type=click.Path(exists=True, path_type=Path),
    default=DEFAULT_SITES,
    show_default=True,
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory for {sample_id}_cpg_depth.tsv.gz",
)
@click.option("--min-mapq", default=20, show_default=True)
@click.option("--min-bq", default=13, show_default=True)
def main(
    bam: Path,
    sample_id: str,
    sites: Path,
    output_dir: Path,
    min_mapq: int,
    min_bq: int,
) -> None:
    """Compute per-site raw depth for one sample."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{sample_id}_cpg_depth.tsv.gz"

    console.print(f"[bold]Sample[/bold] {sample_id}")
    console.print(f"  BAM   : {bam}")
    console.print(f"  index : {bam_index_path(bam)}")
    console.print(f"  sites : {sites}")
    console.print(f"  out   : {out_path}")

    site_df = load_sites(sites)
    depth_df = compute_depths(bam, site_df, min_mapq, min_bq)
    depth_df.write_csv(out_path, separator="\t")

    covered = depth_df.filter(pl.col("raw_depth") > 0).height
    mean_depth = depth_df["raw_depth"].mean()
    console.print(
        f"[green]Wrote {out_path}[/green] "
        f"({covered:,}/{depth_df.height:,} sites with depth > 0, "
        f"mean raw depth {mean_depth:.2f})"
    )


if __name__ == "__main__":
    main()
