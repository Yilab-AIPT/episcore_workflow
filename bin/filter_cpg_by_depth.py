#!/usr/bin/env python3
"""Build a CpG site list filtered by raw_total depth (target + background).

Matches the NIPT depth filter used by ``extract_beta_value`` + ``beta_to_episcore``:
merge MethylDackel target and background bedGraphs, compute
``raw_total = target_meth+target_unmeth+background_meth+background_unmeth``,
and keep sites with ``raw_total > --depth``.

Output is a TSV of ``chr/start/end`` (end already shifted -1 to match CpG-list
convention) for downstream ``calc_episcore_flexible.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


def read_bedgraph(path: Path, prefix: str) -> pd.DataFrame:
    """Read MethylDackel CpG bedGraph; shift end by -1 to match CpG lists."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=["chr", "start", "end", "beta", f"{prefix}_meth", f"{prefix}_unmeth"],
        dtype={
            "chr": str,
            "start": np.int64,
            "end": np.int64,
            "beta": float,
            f"{prefix}_meth": np.int64,
            f"{prefix}_unmeth": np.int64,
        },
    )
    df["end"] = df["end"] - 1
    return df[["chr", "start", "end", f"{prefix}_meth", f"{prefix}_unmeth"]]


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--target-bedgraph", required=True, type=click.Path(exists=True))
@click.option("--background-bedgraph", required=True, type=click.Path(exists=True))
@click.option("--depth", required=True, type=int, help="Keep sites with raw_total > depth.")
@click.option("--output-prefix", required=True, type=str)
def main(
    target_bedgraph: str,
    background_bedgraph: str,
    depth: int,
    output_prefix: str,
) -> None:
    """Write ``{prefix}_depth_filtered_cpgs.tsv`` (chr/start/end)."""
    try:
        console.rule("[bold blue]Filter CpGs by raw_total depth")
        target = read_bedgraph(Path(target_bedgraph), "target")
        background = read_bedgraph(Path(background_bedgraph), "background")
        console.print(f"  target sites      : {len(target):,}")
        console.print(f"  background sites  : {len(background):,}")
        console.print(f"  depth threshold   : {depth}")

        merged = target.merge(background, on=["chr", "start", "end"], how="outer").fillna(0)
        for col in ("target_meth", "target_unmeth", "background_meth", "background_unmeth"):
            merged[col] = merged[col].astype(np.int64)
        merged["raw_total"] = (
            merged["target_meth"]
            + merged["target_unmeth"]
            + merged["background_meth"]
            + merged["background_unmeth"]
        )
        filtered = merged.loc[merged["raw_total"] > depth, ["chr", "start", "end"]].copy()
        filtered = filtered.sort_values(["chr", "start", "end"]).reset_index(drop=True)

        out_path = f"{output_prefix}_depth_filtered_cpgs.tsv"
        filtered.to_csv(out_path, sep="\t", index=False)
        console.print(
            f"[green]OK[/green] Wrote {out_path} "
            f"({len(filtered):,} / {len(merged):,} sites pass depth>{depth})"
        )
        console.rule("[bold green]Done")
    except Exception as exc:  # noqa: BLE001
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
