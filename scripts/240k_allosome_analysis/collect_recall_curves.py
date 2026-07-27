#!/usr/bin/env python3
"""Collect chrX episcore / zscore vs recall tables from recall_* job outputs."""

from __future__ import annotations

import re
from pathlib import Path

import click
import pandas as pd
from rich.console import Console

import config as cfg

console = Console()
RECALL_DIR_RE = re.compile(r"^recall_([\d.]+)$")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """Write chrX_episcore_vs_recall.tsv and chrX_zscore_vs_recall.tsv."""
    labels = pd.read_csv(cfg.COHORT_LABELS).set_index("sample")

    # --- episcore ---
    epi_rows = []
    for out_dir in sorted(cfg.EPISCORE_RECALL_DIR.glob("recall_*")):
        m = RECALL_DIR_RE.match(out_dir.name)
        if not m:
            continue
        recall = float(m.group(1))
        path = out_dir / "_analyze_zscore.tsv.gz"
        if not path.is_file():
            console.print(f"[yellow]Skip[/yellow] missing {path}")
            continue
        df = pd.read_csv(path, sep="\t", usecols=["sample", "chrX_s_inter"])
        for _, row in df.iterrows():
            sample = row["sample"]
            lab = labels.loc[sample] if sample in labels.index else None
            epi_rows.append(
                {
                    "sample": sample,
                    "recall": recall,
                    "chrX_episcore": row["chrX_s_inter"],
                    "label": None if lab is None else lab["label"],
                    "cohort": None if lab is None else lab["cohort"],
                }
            )
    epi_df = pd.DataFrame(epi_rows)
    cfg.EPISCORE_COLLECTED.parent.mkdir(parents=True, exist_ok=True)
    epi_df.to_csv(cfg.EPISCORE_COLLECTED, sep="\t", index=False)
    console.print(
        f"[green]Wrote[/green] {cfg.EPISCORE_COLLECTED}  "
        f"(recalls={epi_df['recall'].nunique() if len(epi_df) else 0}, "
        f"rows={len(epi_df)})"
    )

    # --- zscore ---
    z_rows = []
    for out_dir in sorted(cfg.ZSCORE_RECALL_DIR.glob("recall_*")):
        m = RECALL_DIR_RE.match(out_dir.name)
        if not m:
            continue
        recall = float(m.group(1))
        path = out_dir / "_analyze_zscore.tsv.gz"
        if not path.is_file():
            console.print(f"[yellow]Skip[/yellow] missing {path}")
            continue
        df = pd.read_csv(path, sep="\t")
        col = "chrX_zscore" if "chrX_zscore" in df.columns else None
        if col is None:
            console.print(f"[yellow]Skip[/yellow] no chrX_zscore in {path}")
            continue
        for _, row in df.iterrows():
            sample = row["sample"]
            lab = labels.loc[sample] if sample in labels.index else None
            z_rows.append(
                {
                    "sample": sample,
                    "recall": recall,
                    "chrX_zscore": row[col],
                    "label": None if lab is None else lab["label"],
                    "cohort": None if lab is None else lab["cohort"],
                }
            )
    z_df = pd.DataFrame(z_rows)
    z_df.to_csv(cfg.ZSCORE_COLLECTED, sep="\t", index=False)
    console.print(
        f"[green]Wrote[/green] {cfg.ZSCORE_COLLECTED}  "
        f"(recalls={z_df['recall'].nunique() if len(z_df) else 0}, "
        f"rows={len(z_df)})"
    )


if __name__ == "__main__":
    main()
