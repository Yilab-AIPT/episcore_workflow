#!/usr/bin/env python3
"""
Build per-recall SNP site lists from 220k CpG recall filtering.

For a target CpG recall, keep 220k probe regions that overlap at least one CpG in
``220k_cpg_recall_<recall>.txt``, then retain SNP sites from the full 220k panel
SNP list whose positions fall inside those probe regions.

Output TSVs match the VCF-like format of ``snps_in_220kprobes_v2.tsv`` (no header).
"""

from __future__ import annotations

from pathlib import Path

import click
import polars as pl
from rich.console import Console

console = Console()

DEFAULT_RECALL_DIR = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list"
)
DEFAULT_PANEL_BED = (
    "/lustre1/cqyi/syfan/snp_nipt/data/240k_panel/panel_info/220k.probe.bed"
)
DEFAULT_SNP_FILE = (
    "/lustre1/cqyi/syfan/snp_nipt/data/240k_panel/panel_info/snps_in_220kprobes_v2.tsv"
)
DEFAULT_OUT_DIR = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260604-snp_subset"
)


def panel_regions_overlapping_cpgs(panel: pl.DataFrame, cpgs: pl.DataFrame) -> pl.DataFrame:
    """Keep probe regions overlapping at least one CpG (join_asof per chromosome)."""
    probes = (
        panel.filter(pl.col("chr").is_in(cpgs["chr"].unique()))
        .rename({"start": "probe_start", "end": "probe_end"})
        .sort("chr", "probe_start")
    )
    sites = cpgs.rename({"start": "cpg_start", "end": "cpg_end"}).sort("chr", "cpg_start")
    return (
        sites.join_asof(
            probes,
            left_on="cpg_start",
            right_on="probe_start",
            by="chr",
            strategy="backward",
        )
        .filter(
            (pl.col("cpg_start") < pl.col("probe_end"))
            & (pl.col("cpg_end") > pl.col("probe_start"))
        )
        .select("chr", pl.col("probe_start").alias("start"), pl.col("probe_end").alias("end"))
        .unique()
    )


def filter_snps_by_probes(snps: pl.DataFrame, probes: pl.DataFrame) -> pl.DataFrame:
    """Keep SNP rows whose 1-based ``pos`` overlaps at least one probe interval."""
    sites = (
        snps.select("chr", "pos")
        .with_columns(
            (pl.col("pos") - 1).alias("snp_start"),
            pl.col("pos").alias("snp_end"),
        )
        .sort("chr", "snp_start")
    )
    probe_tbl = (
        probes.filter(pl.col("chr").is_in(sites["chr"].unique()))
        .rename({"start": "probe_start", "end": "probe_end"})
        .sort("chr", "probe_start")
    )
    hits = (
        sites.join_asof(
            probe_tbl,
            left_on="snp_start",
            right_on="probe_start",
            by="chr",
            strategy="backward",
        )
        .filter(
            (pl.col("snp_start") < pl.col("probe_end"))
            & (pl.col("snp_end") > pl.col("probe_start"))
        )
        .select("chr", "pos")
        .unique()
    )
    return snps.join(hits, on=["chr", "pos"], how="inner")


def format_recall(recall: float) -> str:
    """Match recall filenames (e.g. 0.10 -> 0.1)."""
    return f"{recall:g}"


def build_one_recall(
    recall: float,
    recall_dir: Path,
    panel_bed: Path,
    snp_file: Path,
    out_dir: Path,
) -> Path:
    recall_str = format_recall(recall)
    cpg_list = recall_dir / f"220k_cpg_recall_{recall_str}.txt"
    out_path = out_dir / f"snp_for_recall_{recall_str}.tsv"

    if not cpg_list.is_file():
        raise FileNotFoundError(f"CpG recall list not found: {cpg_list}")
    if not panel_bed.is_file():
        raise FileNotFoundError(f"Panel BED not found: {panel_bed}")
    if not snp_file.is_file():
        raise FileNotFoundError(f"SNP list not found: {snp_file}")

    out_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold]Recall {recall_str}[/bold]")
    console.print(f"  CpG list : {cpg_list}")
    console.print(f"  panel    : {panel_bed}")
    console.print(f"  SNPs     : {snp_file}")
    console.print(f"  output   : {out_path}")

    autosome_cpgs = pl.read_csv(
        cpg_list,
        separator="\t",
        columns=["chr", "start", "end"],
    )
    panel_220k = pl.read_csv(
        panel_bed,
        separator="\t",
        has_header=False,
        new_columns=["chr", "start", "end"],
    )
    snp_full = pl.read_csv(snp_file, separator="\t", has_header=False)
    chr_col, pos_col = snp_full.columns[0], snp_full.columns[1]
    snp_sites = snp_full.select(
        pl.col(chr_col).alias("chr"),
        pl.col(pos_col).alias("pos"),
    )

    filtered_probes = panel_regions_overlapping_cpgs(panel_220k, autosome_cpgs)
    filtered_snps = filter_snps_by_probes(snp_sites, filtered_probes)
    snp_out = snp_full.join(
        filtered_snps.select("chr", "pos"),
        left_on=[chr_col, pos_col],
        right_on=["chr", "pos"],
        how="semi",
    )

    n_panel = panel_220k.height
    n_probes = filtered_probes.height
    n_snps_in = snp_sites.height
    n_snps_out = filtered_snps.height

    console.print(
        f"  probes retained : {n_probes:,} / {n_panel:,} "
        f"({n_probes / n_panel:.4f})"
    )
    console.print(
        f"  SNPs retained   : {n_snps_out:,} / {n_snps_in:,} "
        f"({n_snps_out / n_snps_in:.4f})"
    )

    snp_out.write_csv(out_path, separator="\t", include_header=False)
    console.print(f"[green]Wrote {out_path}[/green] ({n_snps_out:,} sites)")
    return out_path


@click.command()
@click.option(
    "--recall",
    type=float,
    default=None,
    help="Single recall value (e.g. 0.65). Omit when using --recall-min/--recall-max.",
)
@click.option("--recall-min", type=float, default=0.05, show_default=True)
@click.option("--recall-max", type=float, default=0.65, show_default=True)
@click.option("--recall-step", type=float, default=0.05, show_default=True)
@click.option("--recall-dir", type=click.Path(path_type=Path), default=DEFAULT_RECALL_DIR)
@click.option("--panel-bed", type=click.Path(path_type=Path), default=DEFAULT_PANEL_BED)
@click.option("--snp-file", type=click.Path(path_type=Path), default=DEFAULT_SNP_FILE)
@click.option("--out-dir", type=click.Path(path_type=Path), default=DEFAULT_OUT_DIR)
def main(
    recall: float | None,
    recall_min: float,
    recall_max: float,
    recall_step: float,
    recall_dir: Path,
    panel_bed: Path,
    snp_file: Path,
    out_dir: Path,
) -> None:
    """Build SNP subset TSVs for one or more CpG recall levels."""
    if recall is not None:
        recalls = [recall]
    else:
        n = round((recall_max - recall_min) / recall_step)
        recalls = [round(recall_min + i * recall_step, 10) for i in range(n + 1)]

    console.print(f"Processing {len(recalls)} recall value(s)")
    for r in recalls:
        build_one_recall(r, recall_dir, panel_bed, snp_file, out_dir)


if __name__ == "__main__":
    main()
