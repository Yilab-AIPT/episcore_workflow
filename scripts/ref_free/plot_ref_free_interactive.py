#!/usr/bin/env python3
"""Build interactive Plotly HTML for 48+48 ref_free signal-ratio sweeps."""

from __future__ import annotations

import json
from pathlib import Path

import click
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from rich.console import Console

console = Console()

PALETTE = {"Normal": "#9e9e9e", "trisomy": "#d62728"}
FF_MIN = 0.0092


def _prepare(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["ff_before_mq"] = pd.to_numeric(out["ff_before_mq"], errors="coerce")
    out["is_trisomy"] = out["label"].astype(str).str.startswith("T")
    out["status"] = out["is_trisomy"].map({True: "trisomy", False: "Normal"})
    return out


def _ez_ratio_col(cutoff: float) -> str:
    return f"ezscore_signal_ratio_{cutoff:g}"


def _scatter_traces(
    df: pd.DataFrame,
    y_col: str,
    *,
    visible: bool = True,
    showlegend: bool = False,
) -> list:
    traces = []
    for status, color in PALETTE.items():
        sub = df[df["status"] == status]
        traces.append(
            go.Scatter(
                x=sub["ff_before_mq"],
                y=sub[y_col],
                mode="markers",
                name=status,
                legendgroup=status,
                showlegend=showlegend,
                marker=dict(color=color, size=8, opacity=0.85),
                text=sub["sample"],
                hovertemplate=(
                    "%{text}<br>ff=%{x:.4f}<br>ratio=%{y:.4f}<extra>" + status + "</extra>"
                ),
                visible=visible,
            )
        )
    return traces


def build_figure(
    df: pd.DataFrame,
    ez_cutoffs: list[float],
    title: str,
    subtitle: str,
) -> go.Figure:
    df = _prepare(df)
    default_ez = 3.0 if 3.0 in ez_cutoffs else ez_cutoffs[0]
    ez_col0 = _ez_ratio_col(default_ez)

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=("Episcore", "Zscore", "Ezscore"),
        shared_yaxes=True,
        horizontal_spacing=0.06,
    )

    # Base traces for all samples (indices 0-5): ep N/T, z N/T, ez N/T
    for tr in _scatter_traces(df, "episcore_signal_ratio", showlegend=True):
        fig.add_trace(tr, row=1, col=1)
    for tr in _scatter_traces(df, "zscore_signal_ratio"):
        fig.add_trace(tr, row=1, col=2)
    for tr in _scatter_traces(df, ez_col0):
        fig.add_trace(tr, row=1, col=3)

    # FF-filtered traces (indices 6-11), hidden by default
    df_ff = df[df["ff_before_mq"] > FF_MIN]
    for tr in _scatter_traces(df_ff, "episcore_signal_ratio", visible=False):
        fig.add_trace(tr, row=1, col=1)
    for tr in _scatter_traces(df_ff, "zscore_signal_ratio", visible=False):
        fig.add_trace(tr, row=1, col=2)
    for tr in _scatter_traces(df_ff, ez_col0, visible=False):
        fig.add_trace(tr, row=1, col=3)

    n_base = 6  # all-sample ep/z/ez
    n_ff = 6

    # Precompute ez y-data for all cutoffs (all-sample and ff-filtered)
    ez_y_all = {}
    ez_y_ff = {}
    for c in ez_cutoffs:
        col = _ez_ratio_col(c)
        if col not in df.columns:
            raise click.ClickException(f"Missing column {col}")
        ez_y_all[c] = {
            "Normal": df.loc[~df["is_trisomy"], col].tolist(),
            "trisomy": df.loc[df["is_trisomy"], col].tolist(),
        }
        ez_y_ff[c] = {
            "Normal": df_ff.loc[~df_ff["is_trisomy"], col].tolist(),
            "trisomy": df_ff.loc[df_ff["is_trisomy"], col].tolist(),
        }

    # Slider steps update only ezscore traces (indices 4,5 for all; 10,11 for ff)
    steps = []
    for c in ez_cutoffs:
        # When all-sample visible: update traces 4,5; when ff: 10,11 — update both
        steps.append(
            dict(
                method="restyle",
                args=[
                    {
                        "y": [
                            ez_y_all[c]["Normal"],
                            ez_y_all[c]["trisomy"],
                            ez_y_ff[c]["Normal"],
                            ez_y_ff[c]["trisomy"],
                        ]
                    },
                    [4, 5, 10, 11],
                ],
                label=f"{c:g}",
            )
        )

    # Filter buttons: toggle visibility of all-sample vs ff-filtered groups
    vis_all = [True] * n_base + [False] * n_ff
    vis_ff = [False] * n_base + [True] * n_ff
    updatemenus = [
        dict(
            type="buttons",
            direction="right",
            x=0.0,
            y=1.18,
            xanchor="left",
            yanchor="top",
            buttons=[
                dict(
                    label="All samples",
                    method="update",
                    args=[{"visible": vis_all}],
                ),
                dict(
                    label=f"ff > {FF_MIN * 100:.2f}%",
                    method="update",
                    args=[{"visible": vis_ff}],
                ),
            ],
        )
    ]

    default_step = ez_cutoffs.index(default_ez)
    fig.update_layout(
        title=dict(text=f"{title}<br><sup>{subtitle}</sup>", x=0.5),
        height=520,
        width=1200,
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=-0.22, x=0.5, xanchor="center"),
        margin=dict(t=100, b=80),
        updatemenus=updatemenus,
        sliders=[
            dict(
                active=default_step,
                currentvalue=dict(prefix="ezscore cutoff: "),
                pad=dict(t=30, b=10),
                steps=steps,
                x=0.15,
                len=0.7,
            )
        ],
    )
    fig.update_xaxes(title_text="ff_before_mq", tickformat=".1%")
    fig.update_yaxes(title_text="signal ratio", range=[-0.02, 1.05], row=1, col=1)
    for col in (2, 3):
        fig.update_yaxes(range=[-0.02, 1.05], row=1, col=col)
    return fig


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--result-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Dir containing ref_free_ezscore/ with aggregated TSV + run_config.json",
)
@click.option(
    "--output-html",
    default=None,
    type=click.Path(path_type=Path),
    help="Default: <result-dir>/plots/signal_ratio.html",
)
@click.option("--title", default="48+48 reference-free signal ratio", show_default=True)
def main(result_dir: Path, output_html: Path | None, title: str) -> None:
    ref_dir = result_dir / "ref_free_ezscore"
    scores = ref_dir / "abnormality_signal_ratio.tsv"
    config_path = ref_dir / "run_config.json"
    if not scores.is_file():
        raise click.ClickException(f"Missing {scores}")
    config = json.loads(config_path.read_text()) if config_path.is_file() else {}
    ez_cutoffs = [float(x) for x in config.get("ez_cutoffs", [3.0])]

    df = pd.read_csv(scores, sep="\t")
    mode = config.get("combo_mode", "?")
    if mode == "fixed":
        subtitle = (
            f"fixed combo ep {config.get('ep_threshold')}/{config.get('ep_recall')} | "
            f"z {config.get('z_threshold')}/{config.get('z_recall')}"
        )
    else:
        subtitle = (
            f"filtered combos | ep thr[{config.get('ep_threshold_min')},"
            f"{config.get('ep_threshold_max')}] "
            f"rec[{config.get('ep_recall_min')},{config.get('ep_recall_max')}] | "
            f"z thr[{config.get('z_threshold_min')},{config.get('z_threshold_max')}] "
            f"rec[{config.get('z_recall_min')},{config.get('z_recall_max')}] | "
            f"ez pairs={config.get('n_ez_combos')} ({config.get('ez_pair_mode')})"
        )

    fig = build_figure(df, ez_cutoffs, title=title, subtitle=subtitle)
    out = output_html or (result_dir / "plots" / "signal_ratio.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out), include_plotlyjs="cdn", full_html=True)
    console.print(f"[green]OK[/green] Wrote {out}")


if __name__ == "__main__":
    main()
