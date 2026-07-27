#!/usr/bin/env python3
"""Scatter plots: ff_before_mq vs ezscore signal ratio for each scheme."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Set

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rich.console import Console

console = Console()

SCHEME_ORDER = [
    "baseline",
    "15-site1-J",
    "15-site1-M",
    "20-site1-M",
    "15-site2-J",
    "15-site2-M",
    "20-site2-M",
]
# Draw Normal first so trisomy (red) sits on top
STATUS_ORDER = ["Normal", "trisomy"]
PALETTE = {"Normal": "#9e9e9e", "trisomy": "#d62728"}


def _ratio_col(df: pd.DataFrame, cutoff: float) -> str:
    col = f"ezscore_signal_ratio_{cutoff:g}"
    if col in df.columns:
        return col
    if "ezscore_signal_ratio" in df.columns:
        return "ezscore_signal_ratio"
    raise click.ClickException(f"Missing ezscore ratio column for cutoff={cutoff}")


def _prepare_scheme(df: pd.DataFrame, scheme: str, ycol: str) -> pd.DataFrame:
    out = df.copy()
    out["scheme"] = scheme
    out["sample"] = out["sample"].astype(str)
    out["ff_before_mq"] = pd.to_numeric(out["ff_before_mq"], errors="coerce")
    out["signal_ratio"] = pd.to_numeric(out[ycol], errors="coerce")
    out["is_trisomy"] = out["label"].astype(str).str.startswith("T")
    out["status"] = out["is_trisomy"].map({True: "trisomy", False: "Normal"})
    return out


def _low_signal_trisomy_baseline(
    root: Path,
    ez_cutoff: float,
    ratio_max: float,
) -> Set[str]:
    """Samples that are trisomy with signal_ratio < ratio_max in baseline."""
    path = root / "baseline" / "ref_free_ezscore" / "abnormality_signal_ratio.tsv"
    if not path.is_file():
        raise click.ClickException(f"Need baseline aggregate to filter: {path}")
    df = pd.read_csv(path, sep="\t")
    ycol = _ratio_col(df, ez_cutoff)
    df = _prepare_scheme(df, "baseline", ycol)
    hide = df.loc[df["is_trisomy"] & (df["signal_ratio"] < ratio_max), "sample"]
    return set(hide.astype(str))


def _scatter_status(ax, df: pd.DataFrame) -> None:
    """Plot Normal (gray) then trisomy (red) so red is on top."""
    for status in STATUS_ORDER:
        sub = df[df["status"] == status]
        ax.scatter(
            sub["ff_before_mq"],
            sub["signal_ratio"],
            c=PALETTE[status],
            s=36,
            alpha=0.85,
            label=status,
            edgecolors="none",
            zorder=2 if status == "Normal" else 3,
        )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--ref-free-root",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Directory containing <scheme>/ref_free_ezscore/",
)
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--ez-cutoff", default=4.5, show_default=True, type=float)
@click.option(
    "--hide-baseline-trisomy-ratio-below",
    default=0.2,
    show_default=True,
    type=float,
    help="Hide trisomy samples with baseline signal_ratio below this in all plots.",
)
def main(
    ref_free_root: str,
    output_dir: str,
    ez_cutoff: float,
    hide_baseline_trisomy_ratio_below: float,
) -> None:
    """One scatter per scheme + a combined faceted figure."""
    root = Path(ref_free_root)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    hide_samples = _low_signal_trisomy_baseline(
        root, ez_cutoff, hide_baseline_trisomy_ratio_below
    )
    console.print(
        f"[cyan]Hide[/cyan] {len(hide_samples)} baseline trisomy samples with "
        f"signal_ratio < {hide_baseline_trisomy_ratio_below:g}: "
        f"{sorted(hide_samples)}"
    )

    frames: List[pd.DataFrame] = []
    for scheme in SCHEME_ORDER:
        path = root / scheme / "ref_free_ezscore" / "abnormality_signal_ratio.tsv"
        if not path.is_file():
            console.print(f"[yellow]SKIP[/yellow] missing {path}")
            continue
        df = pd.read_csv(path, sep="\t")
        ycol = _ratio_col(df, ez_cutoff)
        df = _prepare_scheme(df, scheme, ycol)
        n_before = len(df)
        df = df.loc[~df["sample"].isin(hide_samples)].copy()
        frames.append(df)

        fig, ax = plt.subplots(figsize=(6, 4.5))
        _scatter_status(ax, df)
        ax.set_xlabel("ff_before_mq")
        ax.set_ylabel(f"ezscore signal ratio (cutoff={ez_cutoff:g})")
        ax.set_title(
            f"{scheme}  (n={len(df)}; hid {n_before - len(df)} low-signal T)"
        )
        ax.set_ylim(-0.05, 1.05)
        ax.legend(frameon=True)
        fig.tight_layout()
        fig.savefig(out / f"{scheme}_ff_vs_signal_ratio.png", dpi=150)
        plt.close(fig)
        console.print(f"[green]OK[/green] {scheme} → {out / f'{scheme}_ff_vs_signal_ratio.png'}")

    if not frames:
        raise click.ClickException("No scheme aggregates found")

    all_df = pd.concat(frames, ignore_index=True)
    all_df["scheme"] = pd.Categorical(all_df["scheme"], categories=SCHEME_ORDER, ordered=True)
    # Stable draw order: Normal first, trisomy last (on top) within each facet
    all_df["status"] = pd.Categorical(all_df["status"], categories=STATUS_ORDER, ordered=True)
    all_df = all_df.sort_values(["scheme", "status"])

    g = sns.relplot(
        data=all_df,
        x="ff_before_mq",
        y="signal_ratio",
        hue="status",
        hue_order=STATUS_ORDER,
        col="scheme",
        col_wrap=4,
        palette=PALETTE,
        height=3.2,
        aspect=1.1,
        s=40,
        alpha=0.85,
        facet_kws={"sharex": True, "sharey": True},
    )
    # Ensure trisomy markers are drawn above Normal in each facet
    for ax in g.axes.flat:
        for coll in ax.collections:
            label = coll.get_label()
            if label == "trisomy":
                coll.set_zorder(3)
            elif label == "Normal":
                coll.set_zorder(2)
    g.set_axis_labels("ff_before_mq", f"ezscore signal ratio ({ez_cutoff:g})")
    g.set(ylim=(-0.05, 1.05))
    g.fig.suptitle(
        f"Ref-free ezscore signal ratio vs FF (cutoff={ez_cutoff:g}; "
        f"hid baseline T ratio<{hide_baseline_trisomy_ratio_below:g})",
        y=1.02,
    )
    g.savefig(out / "all_schemes_ff_vs_signal_ratio.png", dpi=150)
    plt.close(g.fig)

    all_df.to_csv(out / "all_schemes_signal_ratio.tsv", sep="\t", index=False)
    (out / "hidden_baseline_low_signal_trisomy.txt").write_text(
        "\n".join(sorted(hide_samples)) + ("\n" if hide_samples else "")
    )
    (out / "plot_config.json").write_text(
        json.dumps(
            {
                "ez_cutoff": ez_cutoff,
                "hide_baseline_trisomy_ratio_below": hide_baseline_trisomy_ratio_below,
                "hidden_samples": sorted(hide_samples),
                "n_hidden": len(hide_samples),
                "schemes": [s for s in SCHEME_ORDER if (root / s).is_dir()],
                "n_rows": int(len(all_df)),
            },
            indent=2,
        )
        + "\n"
    )
    console.print(f"[green]OK[/green] combined → {out / 'all_schemes_ff_vs_signal_ratio.png'}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
