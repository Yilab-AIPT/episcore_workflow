#!/usr/bin/env python3
"""Build a ref-40 asset dir that keeps grid-search combos but uses fixed 25-sample ezscore ref."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from grid_search_ref40 import (  # noqa: E402
    CHR_LIST,
    _build_dense,
    _read_sample_list,
    build_ezscore_ref_matrix,
    compute_episcore,
    compute_zscore,
)
from minimize_episcore_panel import load_combo_dict  # noqa: E402

console = Console()

NEXTFLOW_LINKS = (
    "best_combo_zscore.csv",
    "best_reference_matrix.tsv",
    "best_sample_scores_recalc_ezscore.tsv",
)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--source-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Grid-search input dir with parquets + ezscore_ref_samples.txt")
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--min-ff", default=0.01, show_default=True, type=float)
def main(source_dir: str, input_dir: str, output_dir: str, min_ff: float) -> None:
    """Emit asset bundle with fixed ezscore_ref_samples.txt normalization matrix."""
    src = Path(source_dir)
    inp = Path(input_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    ep_best, _, _ = load_combo_dict(src / "best_combo_episcore.csv")
    z_best, _, _ = load_combo_dict(src / "best_combo_zscore.csv")
    ref_samples = _read_sample_list(src / "best_reference_samples.txt")
    ez_samples = _read_sample_list(inp / "ezscore_ref_samples.txt")

    meta = pd.read_csv(inp / "meta.csv").drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")
    ep_df = pd.read_parquet(inp / "episcore_grid_search.parquet")
    z_df = pd.read_parquet(inp / "zscore_grid_search.parquet")

    ep_samples = set(ep_df["sample"].astype(str))
    z_samples = set(z_df["sample"].astype(str))
    meta_samples = set(meta["sample"])
    ff_pass = set(meta.set_index("sample").index[meta["ff_before_mq"] > min_ff].astype(str))
    ref_set, ez_set = set(ref_samples), set(ez_samples)
    universe = sorted(
        (meta_samples & ep_samples & z_samples & ff_pass)
        | (ref_set & meta_samples & ep_samples & z_samples)
        | (ez_set & meta_samples & ep_samples & z_samples)
    )
    sample_index = {s: i for i, s in enumerate(universe)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}

    ep_df = ep_df[ep_df["sample"].astype(str).isin(universe)]
    z_df = z_df[z_df["sample"].astype(str).isin(universe)]

    ref_idx = np.array([sample_index[s] for s in ref_samples if s in sample_index], dtype=np.int64)
    ez_idx = np.array([sample_index[s] for s in ez_samples if s in sample_index], dtype=np.int64)
    if ref_idx.size != len(ref_samples):
        raise click.ClickException("Some reference samples missing from parquet universe")
    if ez_idx.size == 0:
        raise click.ClickException("No ezscore reference samples found in parquet universe")
    if ez_idx.size != len(ez_samples):
        console.print(
            f"[yellow]Warning[/yellow] using {ez_idx.size}/{len(ez_samples)} ezscore ref samples "
            "present in parquet universe"
        )
    used_ez_samples = [s for s in ez_samples if s in sample_index]

    ep_combos, ep_arrays = _build_dense(
        ep_df,
        ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
        sample_index,
        chr_index,
    )
    z_combos, z_arrays = _build_dense(z_df, ["percentage"], sample_index, chr_index)
    episcore_all = compute_episcore(
        ep_arrays[0], ep_arrays[1], ep_arrays[2], ep_arrays[3], ref_idx,
    )
    zscore_all = compute_zscore(z_arrays[0], ref_idx)

    ez_matrix = build_ezscore_ref_matrix(
        episcore_all, zscore_all, ep_combos, z_combos, ep_best, z_best, ez_idx,
    )

    shutil.copyfile(src / "best_combo_episcore.csv", out / "best_combo_episcore.csv")
    for name in NEXTFLOW_LINKS:
        dst = out / name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src / name)

    ez_matrix.to_csv(out / "best_ezscore_ref_20_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    (out / "best_ezscore_ref_20_samples.txt").write_text("\n".join(used_ez_samples) + "\n")

    console.print(f"[green]OK[/green] Wrote {out / 'best_combo_episcore.csv'}")
    console.print(
        f"[green]OK[/green] Wrote {out / 'best_ezscore_ref_20_matrix.tsv'} "
        f"({len(used_ez_samples)} ezscore ref samples)"
    )
    console.print(f"[green]OK[/green] Wrote {out / 'best_ezscore_ref_20_samples.txt'}")
    console.print(f"[green]OK[/green] Symlinked {len(NEXTFLOW_LINKS)} files from {src}")


if __name__ == "__main__":
    main()
