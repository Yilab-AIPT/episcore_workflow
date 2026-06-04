#!/usr/bin/env python3
"""Run XO allosome episcore calculation for one recall level."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent


def _find_episcore_dir() -> Path:
    candidates = []
    for parent in [_SCRIPT_DIR, *_SCRIPT_DIR.parents]:
        candidates.extend(
            [
                parent / "episcore",
                parent / "scripts" / "episcore",
            ]
        )
    for candidate in candidates:
        if (candidate / "episcore_fast_calculator.py").is_file():
            return candidate
    raise ImportError("Could not locate episcore_fast_calculator.py")


_EPISCORE_DIR = _find_episcore_dir()
if str(_EPISCORE_DIR) not in sys.path:
    sys.path.insert(0, str(_EPISCORE_DIR))

from episcore_fast_calculator import episcore_fast_calculator  # noqa: E402


def _write_tsv(df: pd.DataFrame, path: str) -> None:
    compression = "gzip" if path.endswith(".gz") else None
    df.to_csv(path, sep="\t", index=False, float_format="%.6f", compression=compression)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--samples-meta",
    required=True,
    type=click.Path(exists=True),
    help="CSV with columns sample, beta_path, ref_type.",
)
@click.option(
    "--cpg-list",
    required=True,
    type=click.Path(exists=True),
    help="CpG list TSV with chr/start/end columns.",
)
@click.option(
    "--output-prefix",
    required=True,
    type=str,
    help="Output prefix; writes _reference_zscore.tsv.gz and _analyze_zscore.tsv.gz.",
)
@click.option(
    "--chr-list",
    default="1-22,X",
    show_default=True,
    help="Chromosomes to analyze.",
)
@click.option(
    "--depth",
    type=int,
    default=None,
    help="Optional minimum depth filter.",
)
@click.option(
    "--depth-col",
    default="raw_total_count",
    show_default=True,
    help="Column used for depth filtering.",
)
@click.option(
    "--ncpus",
    type=int,
    default=None,
    help="Worker count for parallel sample processing.",
)
def main(
    samples_meta: str,
    cpg_list: str,
    output_prefix: str,
    chr_list: str,
    depth: int | None,
    depth_col: str,
    ncpus: int | None,
) -> None:
    """Compute reference/analyze z-scores for one recall CpG list."""
    samples_df = pd.read_csv(samples_meta)
    reference_df, analyze_df = episcore_fast_calculator(
        samples_df,
        cpg_list,
        chr_list=chr_list,
        depth=depth,
        depth_col=depth_col,
        ncpus=ncpus,
    )

    out_dir = Path(output_prefix.rstrip("/"))
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_path = f"{out_dir}/_reference_zscore.tsv.gz"
    analyze_path = f"{out_dir}/_analyze_zscore.tsv.gz"
    _write_tsv(reference_df, ref_path)
    _write_tsv(analyze_df, analyze_path)

    click.echo(f"Wrote reference matrix: {ref_path} ({reference_df.shape[0]} samples)")
    click.echo(f"Wrote analyze matrix  : {analyze_path} ({analyze_df.shape[0]} samples)")


if __name__ == "__main__":
    main()
