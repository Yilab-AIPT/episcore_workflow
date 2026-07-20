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


def _ez_count_col(cutoff: float) -> str:
    return f"ezscore_abnormal_count_{cutoff:g}"


def _ez_ratio_col(cutoff: float) -> str:
    return f"ezscore_signal_ratio_{cutoff:g}"


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-base", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--total-repeats", default=None, type=int)
def main(output_base: str, total_repeats: int | None) -> None:
    out_root = Path(output_base) / "ref_free_ezscore"
    eval_path = out_root / "eval_samples.tsv"
    config_path = out_root / "run_config.json"
    if not eval_path.is_file() or not config_path.is_file():
        raise click.ClickException(f"Missing outputs under {out_root}")

    eval_info = pd.read_csv(eval_path, sep="\t")
    config = json.loads(config_path.read_text())
    n_eval = len(eval_info)

    repeats = total_repeats if total_repeats is not None else int(config["total_repeats"])
    ep_denom = float(int(config["n_ep_combos"]) * repeats)
    z_denom = float(int(config["n_z_combos"]) * repeats)
    ez_denom = float(int(config["n_ez_combos"]) * repeats)

    ez_cutoffs = [float(x) for x in config.get("ez_cutoffs", [config.get("ez_cutoff", 3.0)])]

    slice_files = sorted(out_root.glob("abnormality_counts_*.tsv"))
    if not slice_files:
        raise click.ClickException(f"No abnormality_counts_*.tsv under {out_root}")

    ep_counts = np.zeros(n_eval, dtype=np.int64)
    z_counts = np.zeros(n_eval, dtype=np.int64)
    ez_counts = {c: np.zeros(n_eval, dtype=np.int64) for c in ez_cutoffs}

    for path in slice_files:
        df = pd.read_csv(path, sep="\t")
        pos = df["eval_pos"].to_numpy(dtype=np.int64)
        ep_counts[pos] += df["episcore_abnormal_count"].to_numpy(dtype=np.int64)
        z_counts[pos] += df["zscore_abnormal_count"].to_numpy(dtype=np.int64)
        for c in ez_cutoffs:
            col = _ez_count_col(c)
            if col not in df.columns:
                # backward compat: single ezscore_abnormal_count column
                if "ezscore_abnormal_count" in df.columns and len(ez_cutoffs) == 1:
                    ez_counts[c][pos] += df["ezscore_abnormal_count"].to_numpy(dtype=np.int64)
                else:
                    raise click.ClickException(f"{path} missing column {col}")
            else:
                ez_counts[c][pos] += df[col].to_numpy(dtype=np.int64)

    result = eval_info.copy()
    result["episcore_abnormal_count"] = ep_counts
    result["episcore_signal_ratio"] = ep_counts / ep_denom
    result["zscore_abnormal_count"] = z_counts
    result["zscore_signal_ratio"] = z_counts / z_denom
    for c in ez_cutoffs:
        result[_ez_count_col(c)] = ez_counts[c]
        result[_ez_ratio_col(c)] = ez_counts[c] / ez_denom
    # Convenience alias at default cutoff 3
    if 3.0 in ez_cutoffs:
        result["ezscore_abnormal_count"] = result[_ez_count_col(3.0)]
        result["ezscore_signal_ratio"] = result[_ez_ratio_col(3.0)]

    out_path = out_root / "abnormality_signal_ratio.tsv"
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")

    summary = {
        "total_repeats": repeats,
        "n_ep_combos": int(config["n_ep_combos"]),
        "n_z_combos": int(config["n_z_combos"]),
        "n_ez_combos": int(config["n_ez_combos"]),
        "ez_cutoffs": ez_cutoffs,
        "episcore_denominator": ep_denom,
        "zscore_denominator": z_denom,
        "ezscore_denominator": ez_denom,
        "n_slice_files": len(slice_files),
        "n_eval_samples": n_eval,
    }
    summary_path = out_root / "aggregate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    console.print(f"[green]OK[/green] Aggregated {len(slice_files)} slice files")
    console.print(f"  -> {out_path}")
    console.print(f"  -> {summary_path}")


if __name__ == "__main__":
    main()
