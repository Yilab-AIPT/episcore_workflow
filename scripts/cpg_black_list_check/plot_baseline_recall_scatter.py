#!/usr/bin/env python3
"""Scatter plots for baseline-0.6 / baseline-0.65 / baseline-deeper-0.65 ref_free."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Set

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rich.console import Console

console = Console()

SCHEMES = ["baseline-0.6", "baseline-0.65", "baseline-deeper-0.65"]
STATUS_ORDER = ["Normal", "trisomy"]
PALETTE = {"Normal": "#9e9e9e", "trisomy": "#d62728"}
MID_RATIO_LO = 0.2
MID_RATIO_HI = 0.95


def _ratio_col(df: pd.DataFrame, cutoff: float) -> str:
    col = f"ezscore_signal_ratio_{cutoff:g}"
    if col in df.columns:
        return col
    if "ezscore_signal_ratio" in df.columns:
        return "ezscore_signal_ratio"
    raise click.ClickException(f"Missing ezscore ratio for cutoff={cutoff}")


def _prepare(df: pd.DataFrame, scheme: str, ycol: str) -> pd.DataFrame:
    out = df.copy()
    out["scheme"] = scheme
    out["sample"] = out["sample"].astype(str)
    out["ff_before_mq"] = pd.to_numeric(out["ff_before_mq"], errors="coerce")
    out["signal_ratio"] = pd.to_numeric(out[ycol], errors="coerce")
    out["is_trisomy"] = out["label"].astype(str).str.startswith("T")
    out["status"] = out["is_trisomy"].map({True: "trisomy", False: "Normal"})
    return out


def _hide_from_baseline065(
    root: Path, ez_cutoff: float, ratio_max: float
) -> Set[str]:
    for scheme in ("baseline-0.65", "baseline"):
        path = root / scheme / "ref_free_ezscore" / "abnormality_signal_ratio.tsv"
        if path.is_file():
            break
    else:
        console.print("[yellow]WARN[/yellow] no baseline for hide filter; showing all")
        return set()
    df = pd.read_csv(path, sep="\t")
    ycol = _ratio_col(df, ez_cutoff)
    df = _prepare(df, "hide_src", ycol)
    return set(
        df.loc[df["is_trisomy"] & (df["signal_ratio"] < ratio_max), "sample"].astype(str)
    )


def _scatter_status(ax, df: pd.DataFrame) -> None:
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


def _n_trisomy_mid_band(df: pd.DataFrame, lo: float = MID_RATIO_LO, hi: float = MID_RATIO_HI) -> int:
    mask = (
        df["is_trisomy"]
        & (df["signal_ratio"] > lo)
        & (df["signal_ratio"] < hi)
    )
    return int(mask.sum())


def _annotate_mid_band(ax, df: pd.DataFrame) -> None:
    n_mid = _n_trisomy_mid_band(df)
    ax.axhspan(MID_RATIO_LO, MID_RATIO_HI, color="#d62728", alpha=0.04, zorder=0)
    ax.axhline(MID_RATIO_LO, color="#d62728", ls="--", lw=0.8, alpha=0.45, zorder=1)
    ax.axhline(MID_RATIO_HI, color="#d62728", ls="--", lw=0.8, alpha=0.45, zorder=1)
    ax.text(
        0.98,
        0.02,
        f"trisomy {MID_RATIO_LO:g}–{MID_RATIO_HI:g}: n={n_mid}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#d62728",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9},
    )


def _legend_outside(ax) -> None:
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        frameon=True,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--ref-free-root", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--ez-cutoff", default=4.5, show_default=True, type=float)
@click.option("--hide-baseline-trisomy-ratio-below", default=0.2, show_default=True, type=float)
def main(
    ref_free_root: str,
    output_dir: str,
    ez_cutoff: float,
    hide_baseline_trisomy_ratio_below: float,
) -> None:
    root = Path(ref_free_root)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    hide = _hide_from_baseline065(root, ez_cutoff, hide_baseline_trisomy_ratio_below)
    console.print(f"[cyan]Hide[/cyan] {len(hide)} low-signal trisomies: {sorted(hide)}")

    frames: List[pd.DataFrame] = []
    mid_counts = {}
    for scheme in SCHEMES:
        path = root / scheme / "ref_free_ezscore" / "abnormality_signal_ratio.tsv"
        if not path.is_file():
            raise click.ClickException(f"Missing {path}")
        raw = pd.read_csv(path, sep="\t")
        df = _prepare(raw, scheme, _ratio_col(raw, ez_cutoff))
        df = df.loc[~df["sample"].isin(hide)].copy()
        frames.append(df)
        mid_counts[scheme] = _n_trisomy_mid_band(df)

        fig, ax = plt.subplots(figsize=(6.8, 4.5))
        _scatter_status(ax, df)
        _annotate_mid_band(ax, df)
        ax.set_xlabel("ff_before_mq")
        ax.set_ylabel(f"ezscore signal ratio (cutoff={ez_cutoff:g})")
        ax.set_title(f"{scheme}  (n={len(df)})")
        ax.set_ylim(-0.05, 1.05)
        _legend_outside(ax)
        fig.tight_layout()
        fig.savefig(out / f"{scheme}_ff_vs_signal_ratio.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        console.print(
            f"[green]OK[/green] {out / f'{scheme}_ff_vs_signal_ratio.png'}  "
            f"trisomy({MID_RATIO_LO:g},{MID_RATIO_HI:g})={mid_counts[scheme]}"
        )

    all_df = pd.concat(frames, ignore_index=True)
    all_df["scheme"] = pd.Categorical(all_df["scheme"], categories=SCHEMES, ordered=True)
    all_df["status"] = pd.Categorical(all_df["status"], categories=STATUS_ORDER, ordered=True)
    all_df = all_df.sort_values(["scheme", "status"])

    n_panel = len(SCHEMES)
    fig, axes = plt.subplots(1, n_panel, figsize=(5.2 * n_panel + 1.5, 4.5), sharex=True, sharey=True)
    if n_panel == 1:
        axes = [axes]
    for ax, scheme in zip(axes, SCHEMES):
        sub = all_df[all_df["scheme"] == scheme]
        _scatter_status(ax, sub)
        _annotate_mid_band(ax, sub)
        ax.set_title(f"{scheme}  (n={len(sub)})")
        ax.set_xlabel("ff_before_mq")
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel(f"ezscore signal ratio (cutoff={ez_cutoff:g})")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=True,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
    )
    fig.suptitle(
        f"Baseline recall 0.6 / 0.65 / deeper-0.65 (ez≥{ez_cutoff:g})",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 0.92, 1])
    combo_name = "baseline_0.6_vs_0.65_vs_deeper_ff_vs_signal_ratio.png"
    fig.savefig(out / combo_name, dpi=150, bbox_inches="tight")
    # keep previous filename as alias
    fig.savefig(out / "baseline_0.6_vs_0.65_ff_vs_signal_ratio.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    all_df.to_csv(out / "baseline_recall_signal_ratio.tsv", sep="\t", index=False)
    (out / "plot_config.json").write_text(
        json.dumps(
            {
                "ez_cutoff": ez_cutoff,
                "schemes": SCHEMES,
                "hidden_samples": sorted(hide),
                "n_rows": int(len(all_df)),
                "trisomy_mid_band": {
                    "lo": MID_RATIO_LO,
                    "hi": MID_RATIO_HI,
                    "counts": mid_counts,
                },
            },
            indent=2,
        )
        + "\n"
    )
    console.print(f"[green]OK[/green] {out / combo_name}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
