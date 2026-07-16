#!/usr/bin/env python3
"""Aggregate ref_free_ezscore slice outputs into per-sample signal ratios."""

from __future__ import annotations

import json
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-base", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--total-repeats", default=None, type=int)
def main(output_base: str, total_repeats: int | None) -> None:
    out_root = Path(output_base) / "ref_free_ezscore"
    eval_path = out_root / "eval_samples.tsv"
    config_path = out_root / "run_config.json"
    if not eval_path.is_file() or not config_path.is_file():
        raise click.ClickException(f"Missing outputs under {out_root}; run ref_free_ezscore.py first")

    eval_info = pd.read_csv(eval_path, sep="\t")
    config = json.loads(config_path.read_text())
    n_eval = len(eval_info)

    repeats = total_repeats if total_repeats is not None else int(config["total_repeats"])
    ep_denom = float(int(config["n_ep_combos"]) * repeats)
    z_denom = float(int(config["n_z_combos"]) * repeats)
    ez_denom = float(int(config["n_ez_combos"]) * repeats)

    slice_files = sorted(out_root.glob("abnormality_counts_*.tsv"))
    if not slice_files:
        raise click.ClickException(f"No abnormality_counts_*.tsv under {out_root}")

    ep_counts = np.zeros(n_eval, dtype=np.int64)
    z_counts = np.zeros(n_eval, dtype=np.int64)
    ez_counts = np.zeros(n_eval, dtype=np.int64)
    for path in slice_files:
        df = pd.read_csv(path, sep="\t")
        pos = df["eval_pos"].to_numpy(dtype=np.int64)
        ep_counts[pos] += df["episcore_abnormal_count"].to_numpy(dtype=np.int64)
        z_counts[pos] += df["zscore_abnormal_count"].to_numpy(dtype=np.int64)
        ez_counts[pos] += df["ezscore_abnormal_count"].to_numpy(dtype=np.int64)

    result = eval_info.copy()
    result["episcore_abnormal_count"] = ep_counts
    result["episcore_signal_ratio"] = ep_counts / ep_denom
    result["zscore_abnormal_count"] = z_counts
    result["zscore_signal_ratio"] = z_counts / z_denom
    result["ezscore_abnormal_count"] = ez_counts
    result["ezscore_signal_ratio"] = ez_counts / ez_denom

    out_path = out_root / "abnormality_signal_ratio.tsv"
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")

    summary = {
        "total_repeats": repeats,
        "n_ep_combos": int(config["n_ep_combos"]),
        "n_z_combos": int(config["n_z_combos"]),
        "n_ez_combos": int(config["n_ez_combos"]),
        "episcore_denominator": ep_denom,
        "zscore_denominator": z_denom,
        "ezscore_denominator": ez_denom,
        "n_slice_files": len(slice_files),
        "n_eval_samples": n_eval,
        "mean_episcore_signal_ratio": float(result["episcore_signal_ratio"].mean()),
        "mean_zscore_signal_ratio": float(result["zscore_signal_ratio"].mean()),
        "mean_ezscore_signal_ratio": float(result["ezscore_signal_ratio"].mean()),
    }
    summary_path = out_root / "aggregate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    console.print(f"[green]OK[/green] Aggregated {len(slice_files)} slice files")
    console.print(f"  denominators   : ep={ep_denom:.0f} z={z_denom:.0f} ez={ez_denom:.0f}")
    console.print(f"  -> {out_path}")
    console.print(f"  -> {summary_path}")


if __name__ == "__main__":
    main()
