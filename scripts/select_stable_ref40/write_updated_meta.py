#!/usr/bin/env python3
"""
Write an updated meta samplesheet with ref_40 as early_ref and recalculated zscores.

- ref_type: 'early_ref' for ref_40 samples; previous early_ref that are not in
  ref_40 become 'analyze'; other ref_types left unchanged unless they are in ref_40
- beta_zscores / rc_zscores / final_zscores from ref40_score.tsv
- pred_label from ezscore (cutoff 4.5); set to empty for early_ref rows
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from rich.console import Console

console = Console()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--meta-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Original meta samplesheet",
)
@click.option(
    "--ref40-samples",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="ref40_samples.txt or .tsv with a sample column / one ID per line",
)
@click.option(
    "--score-tsv",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="ref40_score.tsv from select_ref40.py",
)
@click.option(
    "--output-csv",
    required=True,
    type=click.Path(dir_okay=False),
    help="Updated meta samplesheet path",
)
@click.option(
    "--clear-ref-pred-label/--keep-ref-pred-label",
    default=True,
    show_default=True,
    help="Clear pred_label for early_ref rows (matches original convention)",
)
def main(
    meta_csv: str,
    ref40_samples: str,
    score_tsv: str,
    output_csv: str,
    clear_ref_pred_label: bool,
) -> None:
    """Build updated meta samplesheet for ref_40 replacement."""
    meta = pd.read_csv(meta_csv)
    meta["sample"] = meta["sample"].astype(str)

    ref_path = Path(ref40_samples)
    if ref_path.suffix.lower() in {".tsv", ".csv"}:
        ref_df = pd.read_csv(ref_path, sep="\t" if ref_path.suffix.lower() == ".tsv" else ",")
        if "sample" not in ref_df.columns:
            raise click.ClickException(f"No sample column in {ref_path}")
        ref_set = set(ref_df["sample"].astype(str))
    else:
        ref_set = {
            line.strip()
            for line in ref_path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        }

    scores = pd.read_csv(score_tsv, sep="\t")
    scores["sample"] = scores["sample"].astype(str)
    score_cols = ["sample", "beta_zscores", "rc_zscores", "final_zscores", "pred_label"]
    missing = [c for c in score_cols if c not in scores.columns]
    if missing:
        raise click.ClickException(f"score tsv missing columns: {missing}")

    out = meta.merge(scores[score_cols], on="sample", how="left", suffixes=("", "_new"))
    for col in ("beta_zscores", "rc_zscores", "final_zscores", "pred_label"):
        new_col = f"{col}_new"
        if new_col in out.columns:
            out[col] = out[new_col].where(out[new_col].notna(), out[col])
            out = out.drop(columns=[new_col])

    # Update ref_type
    old_early = out["ref_type"].astype(str) == "early_ref"
    in_ref40 = out["sample"].isin(ref_set)
    out.loc[old_early & ~in_ref40, "ref_type"] = "analyze"
    out.loc[in_ref40, "ref_type"] = "early_ref"

    if clear_ref_pred_label:
        out.loc[out["ref_type"].astype(str) == "early_ref", "pred_label"] = pd.NA

    out_path = Path(output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    console.print(f"[green]OK[/green] Wrote {out_path}")
    console.print(
        f"  early_ref now: {(out['ref_type'] == 'early_ref').sum()} "
        f"(expected {len(ref_set)})"
    )
    console.print(
        f"  samples with updated scores: {out['sample'].isin(scores['sample']).sum()}"
    )


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
