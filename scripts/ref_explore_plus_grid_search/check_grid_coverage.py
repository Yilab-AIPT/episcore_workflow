#!/usr/bin/env python3
"""Audit episcore / zscore grid-search parquet coverage and list missing cells.

Writes:
    missing_coverage.tsv   score_type, sample, threshold, recall
    coverage_summary.tsv   per-score combo counts / missing totals
"""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
from rich.console import Console

from grid_coverage import (
    N_AUTOSOMES,
    find_missing_coverage,
    majority_combos,
)

console = Console()

DEFAULT_INPUT = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng"
)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input-dir",
    default=DEFAULT_INPUT,
    show_default=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output-dir",
    default=None,
    type=click.Path(file_okay=False, path_type=Path),
    help="Default: <input-dir>/coverage_check",
)
@click.option("--majority-frac", default=0.95, show_default=True, type=float)
@click.option(
    "--fail-on-missing/--no-fail-on-missing",
    default=False,
    show_default=True,
    help="Exit non-zero when any majority-combo cell is missing",
)
def main(
    input_dir: Path,
    output_dir: Path | None,
    majority_frac: float,
    fail_on_missing: bool,
) -> None:
    """Check episcore/zscore parquet coverage against majority combos."""
    out = output_dir or (input_dir / "coverage_check")
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(input_dir / "meta.csv").drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)

    console.print("[cyan]Loading parquets ...[/cyan]")
    ep = pd.read_parquet(input_dir / "episcore_grid_search.parquet")
    z = pd.read_parquet(input_dir / "zscore_grid_search.parquet")

    ep_samples = set(ep["sample"].astype(str).unique())
    z_samples = set(z["sample"].astype(str).unique())
    universe = sorted(set(meta["sample"]) & ep_samples & z_samples)
    if not universe:
        raise click.ClickException("No shared samples across meta ∩ episcore ∩ zscore")

    summary_rows = []
    missing_frames = []

    for score_type, df in (("episcore", ep), ("zscore", z)):
        combos = majority_combos(df, universe, majority_frac=majority_frac)
        missing = find_missing_coverage(df, universe, combos)
        missing = missing.assign(score_type=score_type)
        missing_frames.append(missing)
        summary_rows.append(
            {
                "score_type": score_type,
                "n_universe": len(universe),
                "n_majority_combos": len(combos),
                "n_expected_cells": len(universe) * len(combos),
                "n_missing_cells": len(missing),
                "n_missing_samples": int(missing["sample"].nunique()) if len(missing) else 0,
            }
        )
        console.print(
            f"  {score_type}: majority_combos={len(combos)} "
            f"missing_cells={len(missing)} "
            f"missing_samples={summary_rows[-1]['n_missing_samples']}"
        )

    missing_all = pd.concat(missing_frames, ignore_index=True)
    if not missing_all.empty:
        missing_all = missing_all[["score_type", "sample", "threshold", "recall"]]
        missing_all = missing_all.sort_values(
            ["score_type", "sample", "threshold", "recall"]
        )
    else:
        missing_all = pd.DataFrame(
            columns=["score_type", "sample", "threshold", "recall"]
        )

    missing_path = out / "missing_coverage.tsv"
    summary_path = out / "coverage_summary.tsv"
    missing_all.to_csv(missing_path, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    (out / "universe_samples.txt").write_text("\n".join(universe) + "\n")
    console.print(f"[green]OK[/green] Wrote {missing_path}")
    console.print(f"[green]OK[/green] Wrote {summary_path}")
    console.print(f"  n_chr required : {N_AUTOSOMES}")
    console.print(f"  universe       : {len(universe)}")

    if fail_on_missing and len(missing_all):
        raise click.ClickException(
            f"Coverage gaps remain: {len(missing_all)} missing cells "
            f"(see {missing_path})"
        )


if __name__ == "__main__":
    main()
