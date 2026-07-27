#!/usr/bin/env python3
"""Build a cpg_recall_dir whose episcore 0.65 list is deeper (replaced) sites.

Default bed: ``replaced_deeper_recall_0.65_sites.bed`` — same site count as
production recall 0.65, but ~15% of sites swapped for deeper (higher-depth)
sites from the recall-0.6-only set. All sites ⊂ production recall 0.6, so
``meandiff`` is joined from ``220k_cpg_recall_0.6.txt``.

Other recall files (incl. zscore 0.95) are symlinked from production.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from rich.console import Console

console = Console()

DEFAULT_DEEPER_CANDIDATES = [
    "replaced_deeper_recall_0.65_sites.bed",
    "deeper_recall_0.65_sites.bed",
]
FALLBACK_SMALL_PANEL = Path(
    "/lustre1/cqyi/AIPT_2.0/results/small_panel/replaced_deeper_recall_0.65_sites.bed"
)


def _load_bed_sites(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    df["chr"] = df["chr"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df.drop_duplicates(["chr", "start", "end"])


def _resolve_deeper_bed(input_dir: Path, deeper_bed: Optional[str]) -> Path:
    if deeper_bed:
        p = Path(deeper_bed)
        if not p.is_file():
            raise click.ClickException(f"Missing --deeper-bed {p}")
        return p
    for name in DEFAULT_DEEPER_CANDIDATES:
        p = input_dir / name
        if p.is_file() or p.is_symlink():
            # Prefer replaced_* when both exist
            if name.startswith("replaced_") or not (input_dir / DEFAULT_DEEPER_CANDIDATES[0]).exists():
                if name.startswith("replaced_"):
                    return p
    # explicit preference order
    for name in DEFAULT_DEEPER_CANDIDATES:
        p = input_dir / name
        if p.is_file() or p.is_symlink():
            resolved = p.resolve() if p.is_symlink() else p
            # if only raw deeper is in input, prefer small_panel replaced
            if name == "deeper_recall_0.65_sites.bed" and FALLBACK_SMALL_PANEL.is_file():
                console.print(
                    f"[yellow]NOTE[/yellow] input has raw deeper; using replaced panel:\n"
                    f"  {FALLBACK_SMALL_PANEL}"
                )
                return FALLBACK_SMALL_PANEL
            return resolved
    if FALLBACK_SMALL_PANEL.is_file():
        return FALLBACK_SMALL_PANEL
    raise click.ClickException(
        "No deeper bed found (looked for replaced_deeper_recall_0.65_sites.bed / "
        "deeper_recall_0.65_sites.bed under input-dir and small_panel)."
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False),
              help="Recall dir root (writes recall_dirs/<scheme>/).")
@click.option("--scheme-name", default="baseline-deeper-0.65", show_default=True)
@click.option(
    "--deeper-bed",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Override deeper sites BED (default: replaced_deeper_recall_0.65).",
)
def main(input_dir: str, output_dir: str, scheme_name: str, deeper_bed: Optional[str]) -> None:
    inp = Path(input_dir)
    out_root = Path(output_dir)
    recall_src = inp / "recall_list"
    r06_path = recall_src / "220k_cpg_recall_0.6.txt"
    r065_path = recall_src / "220k_cpg_recall_0.65.txt"
    deeper_path = _resolve_deeper_bed(inp, deeper_bed)
    for p in (r06_path, r065_path):
        if not p.is_file():
            raise click.ClickException(f"Missing {p}")

    deeper = _load_bed_sites(deeper_path)
    r06 = pd.read_csv(r06_path, sep="\t")
    r065 = pd.read_csv(r065_path, sep="\t")
    keys = ["chr", "start", "end"]
    for df in (r06, r065):
        df["chr"] = df["chr"].astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

    merged = deeper.merge(r06, on=keys, how="left", validate="one_to_one")
    n_miss = int(merged["meandiff"].isna().sum())
    if n_miss:
        raise click.ClickException(
            f"{n_miss}/{len(merged)} deeper sites missing from production recall 0.6. "
            "Use replaced_deeper_recall_0.65_sites.bed (depth-swapped panel)."
        )

    out_df = merged[r06.columns]
    dst = out_root / "recall_dirs" / scheme_name
    dst.mkdir(parents=True, exist_ok=True)

    for path in sorted(recall_src.glob("220k_cpg_recall_*.txt")):
        link = dst / path.name
        if link.exists() or link.is_symlink():
            link.unlink()
        if path.name == "220k_cpg_recall_0.65.txt":
            continue
        link.symlink_to(path.resolve())

    out_path = dst / "220k_cpg_recall_0.65.txt"
    out_df.to_csv(out_path, sep="\t", index=False)

    r065_sites = set(zip(r065["chr"], r065["start"], r065["end"]))
    d_sites = set(zip(deeper["chr"], deeper["start"], deeper["end"]))
    stats = {
        "scheme": scheme_name,
        "deeper_bed": str(deeper_path),
        "n_deeper": len(deeper),
        "n_production_0.65": len(r065),
        "n_production_0.6": len(r06),
        "n_written_0.65": len(out_df),
        "same_n_as_production_0.65": len(deeper) == len(r065),
        "symmetric_diff_vs_production_0.65": len(d_sites ^ r065_sites),
        "deeper_only": len(d_sites - r065_sites),
        "production_0.65_only": len(r065_sites - d_sites),
        "output_dir": str(dst),
    }
    (dst / "DEEPER_RECALL_NOTES.json").write_text(json.dumps(stats, indent=2) + "\n")
    console.print(
        f"[green]OK[/green] {dst}\n"
        f"  bed={deeper_path.name}\n"
        f"  deeper n={stats['n_deeper']}  prod 0.65 n={stats['n_production_0.65']}  "
        f"|symdiff|={stats['symmetric_diff_vs_production_0.65']}  "
        f"same_n={stats['same_n_as_production_0.65']}"
    )


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
