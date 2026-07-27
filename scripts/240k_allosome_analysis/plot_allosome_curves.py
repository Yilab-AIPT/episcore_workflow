#!/usr/bin/env python3
"""Plot chrY–FF scatter and chrX episcore/zscore–recall curves.

chrY–FF: old=circle, new=star; color by label.
Recall plots: one curve per sample; solid=old, dotted=new; color by label.
"""

from __future__ import annotations

from pathlib import Path

import click
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from rich.console import Console

import config as cfg

console = Console()

LABEL_COLOR_MAP = {
    "female": "#9E9E9E",
    "male": "#1F77B4",
    "XO": "#E74C3C",
    "XXX": "#8E44AD",
    "XXY": "#E67E22",
    "XX/XY嵌合体": "#16A085",
    "69, XYY": "#2C3E50",
    "T5,XXY": "#D35400",
}


def _color_for(label: str) -> str:
    if label in LABEL_COLOR_MAP:
        return LABEL_COLOR_MAP[label]
    palette = px.colors.qualitative.Dark24
    return palette[hash(str(label)) % len(palette)]


def _scatter_old_new(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    ylabel: str,
    x_as_pct: bool = False,
) -> go.Figure:
    """Circles for cohort=old, stars for cohort=new; color by label."""
    fig = go.Figure()
    labels = sorted(df["label"].dropna().astype(str).unique(), key=str)
    for label in labels:
        color = _color_for(label)
        for cohort, symbol, size, name_suffix in (
            ("old", "circle", 11, "old"),
            ("new", "star", 16, "new"),
        ):
            sub = df[(df["label"].astype(str) == label) & (df["cohort"] == cohort)]
            if sub.empty:
                continue
            fig.add_trace(
                go.Scatter(
                    x=sub[x],
                    y=sub[y],
                    mode="markers",
                    name=f"{label} ({name_suffix})",
                    legendgroup=str(label),
                    marker=dict(
                        symbol=symbol,
                        size=size,
                        color=color,
                        line=dict(width=1.0, color="white"),
                        opacity=0.9,
                    ),
                    customdata=np.stack([sub["sample"].astype(str)], axis=-1),
                    hovertemplate=(
                        f"label={label}<br>cohort={cohort}<br>"
                        f"{x}=%{{x}}<br>{y}=%{{y}}<br>"
                        "sample=%{customdata[0]}<extra></extra>"
                    ),
                )
            )

    if x == "ff_before_mq" and y == "chrY_ratio":
        # Fit lines span full plot x-range (left pad → right edge of data).
        all_x = df[x].to_numpy(dtype=float)
        all_x = all_x[np.isfinite(all_x)]
        x_right = float(all_x.max()) if all_x.size else 0.0
        x_left = float(all_x.min()) if all_x.size else 0.0
        for gender, dash in (("male", "dot"), ("female", "dash")):
            sub = df[df["label"].astype(str) == gender]
            if len(sub) < 2:
                continue
            xv = sub[x].to_numpy(dtype=float)
            yv = sub[y].to_numpy(dtype=float)
            mask = np.isfinite(xv) & np.isfinite(yv)
            if mask.sum() < 2:
                continue
            slope, intercept = np.polyfit(xv[mask], yv[mask], 1)
            xr = np.linspace(x_left, x_right, 100)
            fig.add_trace(
                go.Scatter(
                    x=xr,
                    y=slope * xr + intercept,
                    mode="lines",
                    line=dict(color=_color_for(gender), dash=dash, width=2),
                    name=f"{gender} fit (slope={slope:.3g})",
                    showlegend=True,
                )
            )

    xaxis: dict = dict(title=xlabel, showgrid=True, gridcolor="#ECECEC")
    if x_as_pct:
        xaxis["tickformat"] = ".1%"
    # Pin x-range to data extents so sex-fit lines reach the plot edges.
    if x == "ff_before_mq" and y == "chrY_ratio":
        all_x = df[x].to_numpy(dtype=float)
        all_x = all_x[np.isfinite(all_x)]
        if all_x.size:
            xaxis["range"] = [float(all_x.min()), float(all_x.max())]
    fig.update_layout(
        title=dict(text=title, x=0.02),
        xaxis=xaxis,
        yaxis=dict(title=ylabel, showgrid=True, gridcolor="#ECECEC"),
        template="plotly_white",
        width=1000,
        height=520,
        legend=dict(bgcolor="rgba(255,255,255,0.85)", bordercolor="#CCC", borderwidth=1),
        plot_bgcolor="#FAFAFA",
    )
    return fig


def _curves_per_sample(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    ylabel: str,
) -> go.Figure:
    """One line per sample: solid=old, dotted=new; color by label."""
    fig = go.Figure()
    legend_seen: set[str] = set()
    # Prefer stable sample order within (label, cohort)
    ordered = df.sort_values(["label", "cohort", "sample", x], kind="mergesort")
    for sample, sub in ordered.groupby("sample", sort=False):
        sub = sub.sort_values(x)
        if sub.empty or sub[y].isna().all():
            continue
        label = str(sub["label"].iloc[0])
        cohort = str(sub["cohort"].iloc[0])
        color = _color_for(label)
        dash = "solid" if cohort == "old" else "dot"
        legend_key = f"{label} ({cohort})"
        show = legend_key not in legend_seen
        if show:
            legend_seen.add(legend_key)
        ff = sub["ff_before_mq"].iloc[0] if "ff_before_mq" in sub.columns else None
        if ff is None or (isinstance(ff, float) and np.isnan(ff)):
            ff_txt = "NA"
        else:
            ff_txt = f"{float(ff):.4%}"
        fig.add_trace(
            go.Scatter(
                x=sub[x],
                y=sub[y],
                mode="lines",
                name=legend_key,
                legendgroup=legend_key,
                showlegend=show,
                line=dict(color=color, dash=dash, width=1.8),
                opacity=0.85,
                hovertemplate=(
                    f"sample={sample}<br>label={label}<br>cohort={cohort}<br>"
                    f"ff_before_mq={ff_txt}<br>"
                    f"{x}=%{{x}}<br>{y}=%{{y}}<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title=dict(text=title, x=0.02),
        xaxis=dict(title=xlabel, showgrid=True, gridcolor="#ECECEC"),
        yaxis=dict(title=ylabel, showgrid=True, gridcolor="#ECECEC"),
        template="plotly_white",
        width=1000,
        height=520,
        legend=dict(bgcolor="rgba(255,255,255,0.85)", bordercolor="#CCC", borderwidth=1),
        plot_bgcolor="#FAFAFA",
    )
    return fig


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """Write interactive HTML plots under OUTPUT_DIR/plots/."""
    cfg.PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    if cfg.CHRY_FF_TSV.is_file():
        chry = pd.read_csv(cfg.CHRY_FF_TSV, sep="\t")
        chry = chry.dropna(subset=["ff_before_mq", "chrY_ratio"])
        fig = _scatter_old_new(
            chry,
            x="ff_before_mq",
            y="chrY_ratio",
            title="chrY ratio vs fetal fraction (old=○, new=★)",
            xlabel="ff_before_mq",
            ylabel="chrY_ratio",
            x_as_pct=True,
        )
        out = cfg.PLOTS_DIR / "chry_ff.html"
        fig.write_html(str(out), include_plotlyjs="cdn")
        console.print(f"[green]Wrote[/green] {out}")
    else:
        console.print(f"[yellow]Skip chrY-FF plot; missing {cfg.CHRY_FF_TSV}[/yellow]")

    ff_map = None
    if cfg.CHRY_FF_TSV.is_file():
        ff_map = (
            pd.read_csv(cfg.CHRY_FF_TSV, sep="\t", usecols=["sample", "ff_before_mq"])
            .drop_duplicates("sample")
            .set_index("sample")["ff_before_mq"]
        )

    if cfg.EPISCORE_COLLECTED.is_file() and cfg.EPISCORE_COLLECTED.stat().st_size > 0:
        epi = pd.read_csv(cfg.EPISCORE_COLLECTED, sep="\t")
        if epi.empty or "chrX_episcore" not in epi.columns:
            console.print(
                f"[yellow]Skip episcore plot; empty/incomplete {cfg.EPISCORE_COLLECTED}[/yellow]"
            )
        else:
            if ff_map is not None:
                epi = epi.copy()
                epi["ff_before_mq"] = epi["sample"].map(ff_map)
            fig = _curves_per_sample(
                epi,
                x="recall",
                y="chrX_episcore",
                title="chrX episcore vs recall (threshold=0.5; solid=old, dotted=new)",
                xlabel="Recall",
                ylabel="chrX episcore",
            )
            fig.update_xaxes(range=[0.005, 0.99])
            out = cfg.PLOTS_DIR / "chrX_episcore_vs_recall.html"
            fig.write_html(str(out), include_plotlyjs="cdn")
            console.print(f"[green]Wrote[/green] {out}")
    else:
        console.print(
            f"[yellow]Skip episcore plot; missing {cfg.EPISCORE_COLLECTED}[/yellow]"
        )

    if cfg.ZSCORE_COLLECTED.is_file() and cfg.ZSCORE_COLLECTED.stat().st_size > 0:
        zdf = pd.read_csv(cfg.ZSCORE_COLLECTED, sep="\t")
        if zdf.empty or "chrX_zscore" not in zdf.columns:
            console.print(
                f"[yellow]Skip zscore plot; empty/incomplete {cfg.ZSCORE_COLLECTED}[/yellow]"
            )
        else:
            if ff_map is not None:
                zdf = zdf.copy()
                zdf["ff_before_mq"] = zdf["sample"].map(ff_map)
            fig = _curves_per_sample(
                zdf,
                x="recall",
                y="chrX_zscore",
                title="chrX zscore vs recall (cutoff=0.85; solid=old, dotted=new)",
                xlabel="Recall",
                ylabel="chrX zscore",
            )
            fig.update_xaxes(range=[0.005, 0.99])
            out = cfg.PLOTS_DIR / "chrX_zscore_vs_recall.html"
            fig.write_html(str(out), include_plotlyjs="cdn")
            console.print(f"[green]Wrote[/green] {out}")
    else:
        console.print(
            f"[yellow]Skip zscore plot; missing {cfg.ZSCORE_COLLECTED}[/yellow]"
        )


if __name__ == "__main__":
    main()
