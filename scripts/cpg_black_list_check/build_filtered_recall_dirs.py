#!/usr/bin/env python3
"""Build per-blacklist cpg_recall_dir trees for aipt_ref_40.

Shared upstream (MethylDackel / depth filter) does not depend on recall lists.
Only scoring reads ``220k_cpg_recall_{recall}.txt``. For each blacklist scheme
we write a recall dir where:

  * episcore recall 0.65  = production 0.65 minus blacklist
  * zscore   recall 0.95  = production 0.95 minus blacklist

A ``baseline`` dir keeps the unmodified production lists. Other recall files
from the source recall_list directory are symlinked so the dir stays complete.
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
EP_RECALL = "0.65"
Z_RECALL = "0.95"


def load_bed_sites(path: Path) -> Set[Site]:
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


def scheme_name(bed: Path) -> str:
    """Short stable directory name from blacklist filename."""
    # e.g. 15-site1-BisCap-20260624-10-1x-J.3k.bed → 15-site1-J
    stem = bed.stem  # drop .bed; may still end with .3k
    stem = stem.replace(".3k", "")
    parts = stem.split("-")
    # keep leading size-siteN and trailing J/M tag
    head = "-".join(parts[:2]) if len(parts) >= 2 else stem
    tag = parts[-1] if parts else "bl"
    return f"{head}-{tag}"


def write_filtered(
    src: Path,
    dst: Path,
    blacklist: Set[Site],
) -> dict:
    df = pd.read_csv(src, sep="\t")
    keys = list(zip(df["chr"].astype(str), df["start"].astype(int), df["end"].astype(int)))
    keep = [k not in blacklist for k in keys]
    out = df.loc[keep]
    dst.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dst, sep="\t", index=False)
    n_src = len(df)
    n_out = len(out)
    return {
        "source": str(src),
        "n_before": n_src,
        "n_after": n_out,
        "n_removed": n_src - n_out,
    }


def link_or_copy_rest(src_dir: Path, dst_dir: Path, overwrite_names: Set[str]) -> None:
    """Symlink every recall file not explicitly overwritten."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(src_dir.glob("220k_cpg_recall_*.txt")):
        if path.name in overwrite_names:
            continue
        link = dst_dir / path.name
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(path.resolve())


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(input_dir: str, output_dir: str) -> None:
    """Create baseline + per-blacklist recall directories."""
    inp = Path(input_dir)
    out = Path(output_dir)
    recall_src = inp / "recall_list"
    bl_dir = inp / "cpg_black_list"
    schemes_root = out / "recall_dirs"
    schemes_root.mkdir(parents=True, exist_ok=True)

    ep_src = recall_src / f"220k_cpg_recall_{EP_RECALL}.txt"
    z_src = recall_src / f"220k_cpg_recall_{Z_RECALL}.txt"
    for p in (ep_src, z_src):
        if not p.is_file():
            raise click.ClickException(f"Missing {p}")

    stats = []
    n_ep = len(pd.read_csv(ep_src, sep="\t", usecols=["chr"]))
    n_z = len(pd.read_csv(z_src, sep="\t", usecols=["chr"]))

    # baseline: unmodified production lists
    baseline = schemes_root / "baseline"
    link_or_copy_rest(recall_src, baseline, set())
    stats.append(
        {
            "scheme": "baseline",
            "blacklist": None,
            "episcore": {"n_before": n_ep, "n_after": n_ep, "n_removed": 0},
            "zscore": {"n_before": n_z, "n_after": n_z, "n_removed": 0},
        }
    )
    console.print(f"[green]OK[/green] baseline → {baseline}")

    for bed in sorted(bl_dir.glob("*.bed")):
        scheme = scheme_name(bed)
        dst_dir = schemes_root / scheme
        bl = load_bed_sites(bed)
        overwrite = {
            f"220k_cpg_recall_{EP_RECALL}.txt",
            f"220k_cpg_recall_{Z_RECALL}.txt",
        }
        link_or_copy_rest(recall_src, dst_dir, overwrite)
        ep_stat = write_filtered(ep_src, dst_dir / f"220k_cpg_recall_{EP_RECALL}.txt", bl)
        z_stat = write_filtered(z_src, dst_dir / f"220k_cpg_recall_{Z_RECALL}.txt", bl)
        # also write filtered 0.6 for offline beta-subset bookkeeping (not used by fixed combo)
        r06 = recall_src / "220k_cpg_recall_0.6.txt"
        if r06.is_file():
            write_filtered(r06, dst_dir / "220k_cpg_recall_0.6.filtered_for_audit.txt", bl)
        stats.append(
            {
                "scheme": scheme,
                "blacklist": bed.name,
                "n_blacklist_sites": len(bl),
                "episcore": ep_stat,
                "zscore": z_stat,
            }
        )
        console.print(
            f"[green]OK[/green] {scheme}: ep remove {ep_stat['n_removed']}, "
            f"z remove {z_stat['n_removed']} → {dst_dir}"
        )

    (out / "recall_filter_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    rows = []
    for s in stats:
        rows.append(
            {
                "scheme": s["scheme"],
                "blacklist": s.get("blacklist"),
                "n_blacklist_sites": s.get("n_blacklist_sites", 0),
                "ep_n_before": s["episcore"]["n_before"],
                "ep_n_after": s["episcore"]["n_after"],
                "ep_n_removed": s["episcore"]["n_removed"],
                "z_n_before": s["zscore"]["n_before"],
                "z_n_after": s["zscore"]["n_after"],
                "z_n_removed": s["zscore"]["n_removed"],
            }
        )
    pd.DataFrame(rows).to_csv(out / "recall_filter_stats.tsv", sep="\t", index=False)

    schemes = [s["scheme"] for s in stats]
    (out / "schemes.txt").write_text("\n".join(schemes) + "\n")
    console.print(f"[green]OK[/green] schemes: {schemes}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
