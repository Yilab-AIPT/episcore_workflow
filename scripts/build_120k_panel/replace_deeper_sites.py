#!/usr/bin/env python3
"""Select samples, summarize normalized CpG depth, and replace shallow recall-0.65 sites.

Workflow:
  1. Randomly select N samples with indexed BAMs from a samplesheet (one BAM per sample).
  2. Load per-sample depth tables and compute mean normalized depth across samples
     (per sample: site depth / mean site depth on recall 0.6).
  3. Set A = recall 0.6 minus recall 0.65; set B = recall 0.65.
  4. Rank-match weakest B sites with strongest A sites; replace when A depth > B depth.
  5. Write replaced_deeper_recall_0.65_sites.bed and a summary log.
"""

from __future__ import annotations

import random
from datetime import datetime, timezone
from pathlib import Path

import click
import polars as pl
from rich.console import Console

from bam_utils import bam_has_index

console = Console()

DEFAULT_RECALL_DIR = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260508-grid_search/recall_list"
)
DEFAULT_SAMPLESHEET = (
    "/lustre1/cqyi/syfan/nipt_article_plot/dev_and_test_mqres_samplesheet.csv"
)
DEFAULT_OUT_DIR = "/lustre1/cqyi/AIPT_2.0/results/small_panel"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)]


def load_recall_sites(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t").select(
        pl.col("chr").cast(pl.Utf8),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.all().exclude(["chr", "start", "end"]),
    )


def sort_sites(df: pl.DataFrame) -> pl.DataFrame:
    chr_rank = {c: i for i, c in enumerate(CHR_ORDER)}
    return df.with_columns(
        pl.col("chr").replace_strict(chr_rank, default=len(CHR_ORDER)).alias("_chr_rank")
    ).sort(["_chr_rank", "start"]).drop("_chr_rank")


def select_samples(
    samplesheet: Path,
    n_samples: int,
    seed: int,
    out_csv: Path,
) -> pl.DataFrame:
    df = pl.read_csv(samplesheet)
    if not {"sample", "clean_bam"}.issubset(df.columns):
        raise ValueError("Samplesheet must contain columns: sample, clean_bam")

    unique = df.group_by("sample").agg(pl.col("clean_bam").first().alias("clean_bam"))
    indexed = unique.filter(
        pl.col("clean_bam").map_elements(
            lambda p: bam_has_index(Path(p)),
            return_dtype=pl.Boolean,
        )
    )
    excluded = unique.filter(
        pl.col("clean_bam").map_elements(
            lambda p: not bam_has_index(Path(p)),
            return_dtype=pl.Boolean,
        )
    )

    console.print(
        f"Indexed BAMs: {indexed.height:,} / {unique.height:,} unique samples "
        f"({excluded.height:,} excluded without .bai)"
    )
    if indexed.height < n_samples:
        raise ValueError(
            f"Only {indexed.height} unique samples with indexed BAMs available, "
            f"requested {n_samples}"
        )

    rng = random.Random(seed)
    indices = rng.sample(range(indexed.height), n_samples)
    selected = indexed[sorted(indices)]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    selected.write_csv(out_csv)
    excluded_path = out_csv.parent / "excluded_samples_no_index.csv"
    excluded.select("sample", "clean_bam").write_csv(excluded_path)
    console.print(f"[green]Wrote {out_csv}[/green] ({selected.height} samples, seed={seed})")
    console.print(f"[yellow]Wrote {excluded_path}[/yellow] ({excluded.height} excluded)")
    return selected


def summarize_normalized_depth(
    selected: pl.DataFrame,
    depth_dir: Path,
) -> pl.DataFrame:
    """Return chr/start/end/mean_norm_depth averaged across selected samples."""
    depth_files = []
    missing = []
    for row in selected.iter_rows(named=True):
        sample_id = row["sample"]
        depth_path = depth_dir / f"{sample_id}_cpg_depth.tsv.gz"
        if not depth_path.is_file():
            missing.append(f"{sample_id}: {depth_path}")
            continue
        depth_files.append(
            pl.read_csv(depth_path, separator="\t").with_columns(
                pl.lit(sample_id).alias("sample")
            )
        )

    if missing:
        msg = "Missing depth files:\n  " + "\n  ".join(missing)
        raise FileNotFoundError(msg)
    if not depth_files:
        raise ValueError("No depth files loaded")

    all_depths = pl.concat(depth_files, how="vertical")
    norm = all_depths.with_columns(
        (pl.col("raw_depth") / pl.col("raw_depth").mean().over("sample")).alias(
            "norm_depth"
        )
    )
    summary = norm.group_by("chr", "start", "end").agg(
        pl.col("norm_depth").mean().alias("mean_norm_depth")
    )
    return summary


def replace_sites(
    set_a: pl.DataFrame,
    set_b: pl.DataFrame,
    depth_summary: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, float | int]]:
    key_cols = ["chr", "start", "end"]
    meta_cols = [c for c in set_b.columns if c not in key_cols]

    a = set_a.join(depth_summary, on=key_cols, how="left").with_columns(
        pl.col("mean_norm_depth").fill_null(0.0)
    )
    b = set_b.join(depth_summary, on=key_cols, how="left").with_columns(
        pl.col("mean_norm_depth").fill_null(0.0)
    )

    n_pairs = min(a.height, b.height)
    a_pairs = a.sort("mean_norm_depth", descending=True).head(n_pairs).with_row_index("_idx")
    b_pairs = b.sort("mean_norm_depth", descending=False).head(n_pairs).with_row_index("_idx")

    paired = b_pairs.join(
        a_pairs.select(
            "_idx",
            *[pl.col(c).alias(f"{c}_donor") for c in key_cols + meta_cols],
            pl.col("mean_norm_depth").alias("mean_norm_depth_donor"),
        ),
        on="_idx",
        how="inner",
    )
    replace_pairs = paired.filter(
        pl.col("mean_norm_depth_donor") > pl.col("mean_norm_depth")
    )
    n_replace = replace_pairs.height

    replacement_map = replace_pairs.select(
        pl.col("chr").alias("_orig_chr"),
        pl.col("start").alias("_orig_start"),
        pl.col("end").alias("_orig_end"),
        *[pl.col(f"{c}_donor").alias(c) for c in key_cols + meta_cols],
    )

    joined = set_b.join(
        replacement_map,
        left_on=key_cols,
        right_on=["_orig_chr", "_orig_start", "_orig_end"],
        how="left",
    )
    out_cols = key_cols + meta_cols
    b_out = joined.select(
        [
            pl.coalesce([pl.col(f"{c}_right"), pl.col(c)]).alias(c)
            for c in out_cols
        ]
    )

    pairs_detail = replace_pairs.select(
        pl.col("chr").alias("b_chr"),
        pl.col("start").alias("b_start"),
        pl.col("end").alias("b_end"),
        pl.col("mean_norm_depth").alias("b_mean_norm_depth"),
        pl.col("chr_donor").alias("a_chr"),
        pl.col("start_donor").alias("a_start"),
        pl.col("end_donor").alias("a_end"),
        pl.col("mean_norm_depth_donor").alias("a_mean_norm_depth"),
    )

    stats: dict[str, float | int] = {
        "n_set_a": set_a.height,
        "n_set_b": set_b.height,
        "n_pairs_compared": n_pairs,
        "n_replaced": n_replace,
        "pct_replaced_in_b": 100.0 * n_replace / set_b.height if set_b.height else 0.0,
    }
    return sort_sites(b_out), pairs_detail, stats


def write_log(
    log_path: Path,
    stats: dict[str, float | int],
    selected: pl.DataFrame,
    seed: int,
    recall_low: float,
    recall_high: float,
    pairs_detail: pl.DataFrame,
) -> None:
    lines = [
        f"# replace_deeper_sites log",
        f"timestamp: {datetime.now(timezone.utc).isoformat()}",
        f"random_seed: {seed}",
        f"sample_selection: indexed BAMs only (.bai required)",
        f"recall_low: {recall_low}",
        f"recall_high: {recall_high}",
        f"n_samples: {selected.height}",
        f"samples: {', '.join(selected['sample'].to_list())}",
        "",
        f"n_set_a (recall {recall_low:g} only): {stats['n_set_a']:,}",
        f"n_set_b (recall {recall_high:g}): {stats['n_set_b']:,}",
        f"n_rank_pairs_compared: {stats['n_pairs_compared']:,}",
        f"n_sites_in_b_replaced_by_a: {stats['n_replaced']:,}",
        f"pct_of_set_b_replaced: {stats['pct_replaced_in_b']:.4f}%",
        "",
        "replacement_strategy: sort set B by mean_norm_depth ascending, "
        "set A descending; pair by rank; replace B with A when A depth > B depth.",
        "",
    ]
    if pairs_detail.height > 0:
        lines.append("top_10_replacements_by_depth_gain:")
        top = pairs_detail.with_columns(
            (pl.col("a_mean_norm_depth") - pl.col("b_mean_norm_depth")).alias("depth_gain")
        ).sort("depth_gain", descending=True).head(10)
        for row in top.iter_rows(named=True):
            lines.append(
                f"  B {row['b_chr']}:{row['b_start']}-{row['b_end']} "
                f"(depth={row['b_mean_norm_depth']:.4f}) -> "
                f"A {row['a_chr']}:{row['a_start']}-{row['a_end']} "
                f"(depth={row['a_mean_norm_depth']:.4f}, gain={row['depth_gain']:.4f})"
            )

    log_path.write_text("\n".join(lines) + "\n")


@click.command()
@click.option("--samplesheet", type=click.Path(exists=True, path_type=Path), default=DEFAULT_SAMPLESHEET)
@click.option("--recall-dir", type=click.Path(exists=True, path_type=Path), default=DEFAULT_RECALL_DIR)
@click.option("--recall-low", type=float, default=0.6, show_default=True)
@click.option("--recall-high", type=float, default=0.65, show_default=True)
@click.option("--n-samples", type=int, default=50, show_default=True)
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--depth-dir", type=click.Path(path_type=Path), default=None)
@click.option("--out-dir", type=click.Path(path_type=Path), default=DEFAULT_OUT_DIR)
@click.option(
    "--selected-samples",
    type=click.Path(path_type=Path),
    default=None,
    help="CSV with sample,clean_bam. If missing, samples are selected and written here.",
)
@click.option(
    "--select-only",
    is_flag=True,
    default=False,
    help="Only select samples and exit (no depth merge / replacement).",
)
def main(
    samplesheet: Path,
    recall_dir: Path,
    recall_low: float,
    recall_high: float,
    n_samples: int,
    seed: int,
    depth_dir: Path | None,
    out_dir: Path,
    selected_samples: Path | None,
    select_only: bool,
) -> None:
    """Merge per-sample depths and build the replaced recall-0.65 site list."""
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if selected_samples is None:
        selected_samples = out_dir / "selected_samples.csv"
    else:
        selected_samples = selected_samples.resolve()

    if depth_dir is None:
        depth_dir = out_dir / "depth_per_sample"
    else:
        depth_dir = depth_dir.resolve()

    if select_only or not selected_samples.is_file():
        selected = select_samples(samplesheet, n_samples, seed, selected_samples)
    else:
        selected = pl.read_csv(selected_samples)
        console.print(f"Using existing sample list: {selected_samples}")

    if select_only:
        return

    recall_low_path = recall_dir / f"220k_cpg_recall_{recall_low:g}.txt"
    recall_high_path = recall_dir / f"220k_cpg_recall_{recall_high:g}.txt"
    recall_low_sites = load_recall_sites(recall_low_path)
    recall_high_sites = load_recall_sites(recall_high_path)

    set_b = recall_high_sites
    set_a = recall_low_sites.join(
        recall_high_sites.select("chr", "start", "end"),
        on=["chr", "start", "end"],
        how="anti",
    )

    console.print(f"Set A (recall {recall_low:g} only): {set_a.height:,} sites")
    console.print(f"Set B (recall {recall_high:g}): {set_b.height:,} sites")

    depth_summary = summarize_normalized_depth(selected, depth_dir)
    depth_summary_path = out_dir / "cpg_mean_normalized_depth.tsv.gz"
    depth_summary.write_csv(depth_summary_path, separator="\t")
    console.print(f"[green]Wrote {depth_summary_path}[/green]")

    replaced, pairs_detail, stats = replace_sites(set_a, set_b, depth_summary)

    bed_path = out_dir / "replaced_deeper_recall_0.65_sites.bed"
    replaced.select("chr", "start", "end").write_csv(
        bed_path, separator="\t", include_header=False
    )
    console.print(f"[green]Wrote {bed_path}[/green] ({replaced.height:,} sites)")

    pairs_path = out_dir / "replacement_pairs.tsv"
    pairs_detail.write_csv(pairs_path, separator="\t")
    console.print(f"[green]Wrote {pairs_path}[/green] ({pairs_detail.height:,} pairs)")

    log_path = out_dir / "replace_deeper_sites.log"
    write_log(
        log_path,
        stats,
        selected,
        seed,
        recall_low,
        recall_high,
        pairs_detail,
    )
    console.print(f"[green]Wrote {log_path}[/green]")
    console.print(
        f"Replaced {stats['n_replaced']:,} / {stats['n_set_b']:,} set-B sites "
        f"({stats['pct_replaced_in_b']:.4f}%)"
    )


if __name__ == "__main__":
    main()
