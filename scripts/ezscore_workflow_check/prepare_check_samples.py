#!/usr/bin/env python3
"""Prepare a filtered samplesheet of check samples for aipt_ref_40 vs meta."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Sequence

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

REQUIRED = "PTAY1183P"
META_COLS = ("beta_zscores", "rc_zscores", "final_zscores")


def _select(
    candidates: Sequence[str],
    *,
    n: int,
    required: str,
    seed: int,
) -> List[str]:
    cand = sorted(set(candidates))
    if required not in cand:
        raise click.ClickException(f"Required sample {required!r} not in candidates")
    others = [s for s in cand if s != required]
    rng = np.random.default_rng(seed)
    extra = (
        [others[i] for i in sorted(rng.choice(len(others), size=n - 1, replace=False).tolist())]
        if n > 1
        else []
    )
    return sorted([required, *extra])


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--n-check", default=10, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--required-sample", default=REQUIRED, show_default=True)
def main(
    input_dir: str,
    output_dir: str,
    n_check: int,
    seed: int,
    required_sample: str,
) -> None:
    """Write check_samples.txt + check_samplesheet.csv."""
    inp = Path(input_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(inp / "meta.csv")
    meta["sample"] = meta["sample"].astype(str)
    sheet = pd.read_csv(inp / "samplesheet.csv")
    sheet["sample"] = sheet["sample"].astype(str)
    sheet_ids = set(sheet["sample"])

    analyze = meta.loc[meta["ref_type"].astype(str) == "analyze"].set_index("sample")
    candidates = [
        s
        for s in analyze.index.astype(str)
        if s in sheet_ids and all(pd.notna(analyze.loc[s, c]) for c in META_COLS)
    ]
    check = _select(candidates, n=n_check, required=required_sample, seed=seed)

    (out / "check_samples.txt").write_text("\n".join(check) + "\n")
    filtered = sheet[sheet["sample"].isin(check)].sort_values(["sample", "clean_bam"])
    filtered.to_csv(out / "check_samplesheet.csv", index=False)

    cfg = {
        "n_check": len(check),
        "check_samples": check,
        "seed": seed,
        "n_samplesheet_rows": int(len(filtered)),
    }
    (out / "prepare_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    console.print(f"[green]OK[/green] check samples: {check}")
    console.print(f"[green]OK[/green] samplesheet rows: {len(filtered)}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
