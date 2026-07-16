#!/usr/bin/env python3
"""
Aggregate reference-free repeat outputs into per-sample abnormality signal ratios.

Scans ``<output-base>/ref_free/abnormality_counts_*.tsv``, sums abnormal counts
across all repeat slices, and divides by combo_number * total_repeats:

    episcore_signal_ratio = episcore_abnormal_count / (n_ep_combos * total_repeats)
    zscore_signal_ratio   = zscore_abnormal_count   / (n_z_combos   * total_repeats)
    either_signal_ratio   = either_abnormal_count   / (n_union_combos * total_repeats)

Writes:
    abnormality_signal_ratio.tsv
    aggregate_summary.json
"""

from __future__ import annotations

import json
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--output-base",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Base dir containing ref_free/ outputs",
)
@click.option(
    "--total-repeats",
    default=None,
    type=int,
    help="Override repeat count (default: read from run_config.json)",
)
def main(output_base: str, total_repeats: int | None) -> None:
    """Aggregate ref_free slice abnormality counts."""
    out_root = Path(output_base) / "ref_free"
    eval_path = out_root / "eval_samples.tsv"
    config_path = out_root / "run_config.json"
    if not eval_path.is_file():
        raise click.ClickException(f"Missing {eval_path}; run ref_free.py first")
    if not config_path.is_file():
        raise click.ClickException(f"Missing {config_path}; run ref_free.py first")

    eval_info = pd.read_csv(eval_path, sep="\t")
    config = json.loads(config_path.read_text())
    n_eval = len(eval_info)

    repeats = total_repeats if total_repeats is not None else int(config["total_repeats"])
    n_ep = int(config["n_ep_combos"])
    n_z = int(config["n_z_combos"])
    n_union = int(config["n_union_combos"])
    ep_denom = float(n_ep * repeats)
    z_denom = float(n_z * repeats)
    either_denom = float(n_union * repeats)

    slice_files = sorted(out_root.glob("abnormality_counts_*.tsv"))
    if not slice_files:
        raise click.ClickException(f"No abnormality_counts_*.tsv under {out_root}")

    ep_counts = np.zeros(n_eval, dtype=np.int64)
    z_counts = np.zeros(n_eval, dtype=np.int64)
    either_counts = np.zeros(n_eval, dtype=np.int64)

    for path in slice_files:
        df = pd.read_csv(path, sep="\t")
        pos = df["eval_pos"].to_numpy(dtype=np.int64)
        ep_counts[pos] += df["episcore_abnormal_count"].to_numpy(dtype=np.int64)
        z_counts[pos] += df["zscore_abnormal_count"].to_numpy(dtype=np.int64)
        either_counts[pos] += df["either_abnormal_count"].to_numpy(dtype=np.int64)

    result = eval_info.copy()
    result["episcore_abnormal_count"] = ep_counts
    result["episcore_signal_ratio"] = ep_counts / ep_denom
    result["zscore_abnormal_count"] = z_counts
    result["zscore_signal_ratio"] = z_counts / z_denom
    result["either_abnormal_count"] = either_counts
    result["either_signal_ratio"] = either_counts / either_denom

    out_path = out_root / "abnormality_signal_ratio.tsv"
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")

    summary = {
        "total_repeats": repeats,
        "n_ep_combos": n_ep,
        "n_z_combos": n_z,
        "n_union_combos": n_union,
        "episcore_denominator": ep_denom,
        "zscore_denominator": z_denom,
        "either_denominator": either_denom,
        "n_slice_files": len(slice_files),
        "n_eval_samples": n_eval,
        "mean_episcore_signal_ratio": float(result["episcore_signal_ratio"].mean()),
        "mean_zscore_signal_ratio": float(result["zscore_signal_ratio"].mean()),
        "mean_either_signal_ratio": float(result["either_signal_ratio"].mean()),
    }
    summary_path = out_root / "aggregate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    console.print(f"[green]OK[/green] Aggregated {len(slice_files)} slice files")
    console.print(f"  denominators   : ep={ep_denom:.0f} z={z_denom:.0f} either={either_denom:.0f}")
    console.print(f"  -> {out_path}")
    console.print(f"  -> {summary_path}")


if __name__ == "__main__":
    main()
