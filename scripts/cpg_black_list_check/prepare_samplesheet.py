#!/usr/bin/env python3
"""Build aipt_ref_40 samplesheet from meta.csv ∩ mqres.csv.

Keeps ``analyze`` samples that have meta ground-truth scores
(``beta_zscores``, ``rc_zscores``, ``final_zscores``). Emits one row per
mqres (sample, clean_bam, deconv_res) entry so PREPARE_INPUTS can merge.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
import pandas as pd
from rich.console import Console

console = Console()

META_COLS = ("beta_zscores", "rc_zscores", "final_zscores")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option(
    "--ref-type",
    default="analyze",
    show_default=True,
    help="meta.ref_type filter; use 'all' to keep every ref_type with scores.",
)
def main(input_dir: str, output_dir: str, ref_type: str) -> None:
    """Write samplesheet.csv + sample_list.txt + prepare_config.json."""
    inp = Path(input_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(inp / "meta.csv")
    meta["sample"] = meta["sample"].astype(str)
    mq = pd.read_csv(inp / "mqres.csv")
    mq["sample"] = mq["sample"].astype(str)

    scored = meta.loc[meta[list(META_COLS)].notna().all(axis=1)].copy()
    if ref_type != "all":
        scored = scored.loc[scored["ref_type"].astype(str) == ref_type]

    keep = set(scored["sample"])
    sheet = mq.loc[mq["sample"].isin(keep)].copy()
    missing = sorted(keep - set(sheet["sample"]))
    if missing:
        raise click.ClickException(
            f"{len(missing)} meta samples missing from mqres.csv "
            f"(examples: {missing[:10]})"
        )

    for col in ("clean_bam", "deconv_res"):
        bad = sheet.loc[~sheet[col].map(lambda p: Path(str(p)).is_file()), "sample"]
        if len(bad):
            raise click.ClickException(
                f"{bad.nunique()} samples have missing {col} "
                f"(examples: {bad.unique()[:5].tolist()})"
            )

    sheet = sheet[["sample", "clean_bam", "deconv_res"]].sort_values(
        ["sample", "clean_bam"]
    )
    samples = sorted(sheet["sample"].unique())
    sheet.to_csv(out / "samplesheet.csv", index=False)
    (out / "sample_list.txt").write_text("\n".join(samples) + "\n")

    cfg = {
        "n_samples": len(samples),
        "n_samplesheet_rows": int(len(sheet)),
        "ref_type": ref_type,
        "meta_score_cols": list(META_COLS),
        "label_counts": scored.set_index("sample")
        .loc[samples, "label"]
        .astype(str)
        .value_counts()
        .to_dict(),
    }
    (out / "prepare_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    console.print(
        f"[green]OK[/green] {len(samples)} samples / {len(sheet)} rows → {out / 'samplesheet.csv'}"
    )


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
