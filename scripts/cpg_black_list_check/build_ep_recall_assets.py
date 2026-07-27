#!/usr/bin/env python3
"""Build fixed-combo grid_search_assets for a chosen episcore recall.

Copies zscore/ezscore assets from an existing original-combo asset dir and
rewrites ``best_combo_episcore.csv`` to ``threshold=0.5, recall=<ep-recall>``.
Zscore combo stays 0.85 / 0.95.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import click
import pandas as pd
from rich.console import Console

console = Console()
CHR_LIST = [f"chr{i}" for i in range(1, 23)]
COPY_FILES = [
    "best_combo_zscore.csv",
    "best_reference_matrix.tsv",
    "best_reference_samples.txt",
    "best_ezscore_ref_20_matrix.tsv",
    "best_ezscore_ref_20_samples.txt",
    "best_sample_scores_recalc_ezscore.tsv",
]


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--source-assets",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Existing original-combo assets (e.g. grid_search_assets).",
)
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--ep-threshold", default=0.5, show_default=True, type=float)
@click.option("--ep-recall", required=True, type=float)
@click.option("--z-threshold", default=0.85, show_default=True, type=float)
@click.option("--z-recall", default=0.95, show_default=True, type=float)
def main(
    source_assets: str,
    output_dir: str,
    ep_threshold: float,
    ep_recall: float,
    z_threshold: float,
    z_recall: float,
) -> None:
    src = Path(source_assets)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for name in COPY_FILES:
        s = src / name
        if not s.is_file():
            raise click.ClickException(f"Missing {s}")
        shutil.copy2(s, out / name)

    # Rewrite episcore combo
    ep = pd.DataFrame(
        {
            "chr": CHR_LIST,
            "threshold": ep_threshold,
            "recall": ep_recall,
            "has_target": False,
            "min_recall": ep_recall,
        }
    )
    ep.to_csv(out / "best_combo_episcore.csv", index=False)

    # Ensure zscore combo matches requested (usually already 0.85/0.95)
    z = pd.DataFrame(
        {
            "chr": CHR_LIST,
            "threshold": z_threshold,
            "recall": z_recall,
            "has_target": False,
            "min_recall": z_recall,
        }
    )
    z.to_csv(out / "best_combo_zscore.csv", index=False)

    (out / "ASSETS_NOTES.txt").write_text(
        f"ep combo: {ep_threshold}/{ep_recall}\n"
        f"z combo: {z_threshold}/{z_recall}\n"
        f"source: {src}\n"
        "reference matrices copied from source (early production refs).\n"
    )
    console.print(
        f"[green]OK[/green] {out}  ep={ep_threshold}/{ep_recall}  z={z_threshold}/{z_recall}"
    )


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
