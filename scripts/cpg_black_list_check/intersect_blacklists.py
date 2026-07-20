#!/usr/bin/env python3
"""Intersect probe-design blacklists with submitted / production CpG panels.

Writes:
  blacklist_source_report.tsv  — per-blacklist overlap stats + inferred source
  blacklist_source_report.json — same as structured JSON
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Set, Tuple

import click
import pandas as pd
from rich.console import Console

console = Console()

Site = Tuple[str, int, int]


def load_bed(path: Path) -> Set[Site]:
    """Load (chr, start, end) from tab- or comma-separated BED (no header)."""
    sites: Set[Site] = set()
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith("chr\t") or line.startswith("chr,"):
            continue
        sep = "," if ("," in line and "\t" not in line) else "\t"
        parts = line.split(sep)
        sites.add((parts[0], int(float(parts[1])), int(float(parts[2]))))
    return sites


def load_recall_txt(path: Path) -> Set[Site]:
    df = pd.read_csv(path, sep="\t", usecols=["chr", "start", "end"])
    return set(
        zip(df["chr"].astype(str), df["start"].astype(int), df["end"].astype(int))
    )


def _frac(n: int, d: int) -> float:
    return float(n) / d if d else 0.0


def infer_source(frac_deeper: float, frac_r06: float) -> str:
    """Heuristic source label from overlap fractions."""
    if frac_deeper >= 0.99 and frac_r06 >= 0.99:
        return "deeper_recall_0.65 (subset of recall_0.6)"
    if frac_r06 >= 0.99 and frac_deeper < 0.99:
        return "recall_0.6 (extends beyond deeper_0.65)"
    if frac_deeper >= 0.99:
        return "deeper_recall_0.65"
    if frac_r06 >= 0.95 and frac_deeper >= 0.90:
        return "mostly deeper_recall_0.65 / recall_0.6"
    return "ambiguous/partial"


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(input_dir: str, output_dir: str) -> None:
    """Report which submitted panel each blacklist most likely came from."""
    inp = Path(input_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    deeper = load_bed(inp / "deeper_recall_0.65_sites.bed")
    recall06 = load_bed(inp / "recall_0.6_sites.bed")
    r065 = load_recall_txt(inp / "recall_list" / "220k_cpg_recall_0.65.txt")
    r06 = load_recall_txt(inp / "recall_list" / "220k_cpg_recall_0.6.txt")
    r095 = load_recall_txt(inp / "recall_list" / "220k_cpg_recall_0.95.txt")

    panel_summary = {
        "deeper_recall_0.65_sites": len(deeper),
        "recall_0.6_sites": len(recall06),
        "recall_list_0.65": len(r065),
        "recall_list_0.6": len(r06),
        "recall_list_0.95": len(r095),
        "deeper_subset_of_recall_0.6_sites": deeper <= recall06,
        "r065_subset_of_r06": r065 <= r06,
        "r095_subset_of_r065": r095 <= r065,
        "deeper_vs_recall_list_0.65_symmetric_diff": len(deeper ^ r065),
        "recall_0.6_sites_vs_recall_list_0.6_symmetric_diff": len(recall06 ^ r06),
    }

    rows = []
    bl_dir = inp / "cpg_black_list"
    for bed in sorted(bl_dir.glob("*.bed")):
        bl = load_bed(bed)
        n = len(bl)
        inter_deeper = len(bl & deeper)
        inter_r06_sites = len(bl & recall06)
        inter_r065 = len(bl & r065)
        inter_r06 = len(bl & r06)
        inter_r095 = len(bl & r095)
        outside_submitted = len(bl - deeper - recall06)
        in_06_not_deeper = len((bl & recall06) - deeper)
        frac_deeper = _frac(inter_deeper, n)
        frac_r06_sites = _frac(inter_r06_sites, n)
        source = infer_source(frac_deeper, frac_r06_sites)
        row = {
            "blacklist": bed.name,
            "scheme": bed.stem,
            "n_sites": n,
            "inter_deeper_0.65": inter_deeper,
            "frac_deeper_0.65": round(frac_deeper, 4),
            "inter_recall_0.6_sites": inter_r06_sites,
            "frac_recall_0.6_sites": round(frac_r06_sites, 4),
            "outside_both_submitted": outside_submitted,
            "in_0.6_sites_not_deeper": in_06_not_deeper,
            "inter_recall_list_0.65": inter_r065,
            "frac_recall_list_0.65": round(_frac(inter_r065, n), 4),
            "n_removed_from_episcore_0.65": inter_r065,
            "inter_recall_list_0.6": inter_r06,
            "frac_recall_list_0.6": round(_frac(inter_r06, n), 4),
            "inter_recall_list_0.95": inter_r095,
            "frac_recall_list_0.95": round(_frac(inter_r095, n), 4),
            "n_removed_from_zscore_0.95": inter_r095,
            "inferred_source": source,
        }
        rows.append(row)
        console.print(
            f"[cyan]{bed.name}[/cyan] n={n} "
            f"∩deeper={inter_deeper} ({frac_deeper:.1%}) "
            f"∩0.6sites={inter_r06_sites} ({frac_r06_sites:.1%}) "
            f"→ remove ep={inter_r065} z={inter_r095} | {source}"
        )

    report = pd.DataFrame(rows)
    report.to_csv(out / "blacklist_source_report.tsv", sep="\t", index=False)
    payload = {"panel_summary": panel_summary, "blacklists": rows}
    (out / "blacklist_source_report.json").write_text(json.dumps(payload, indent=2) + "\n")
    console.print(f"[green]OK[/green] wrote {out / 'blacklist_source_report.tsv'}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
