#!/usr/bin/env python3
"""Bootstrap AUC learning curve from fixed-combo per-repeat flag shards.

For each ``repeat_n`` on a log grid, draw ``n-boot`` random subsets of that many
repeats (with replacement if needed), compute signal ratios = mean(flags), then
ROC-AUC on ff≥ff_min Normal vs Trisomy. Plot mean AUC ± percentile band.
"""

from __future__ import annotations

import json
from pathlib import Path

import click
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from rich.console import Console

from separation import is_trisomy_label, roc_auc

console = Console()

DEFAULT_GRID = "10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000"


def _load_flags(flag_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    shards = sorted(flag_dir.glob("flags_*.npz"))
    if not shards:
        raise click.ClickException(f"No flags_*.npz under {flag_dir}")
    # sort by repeat_start
    def _start(p: Path) -> int:
        return int(p.stem.split("_")[1])

    shards = sorted(shards, key=_start)
    ep_parts, z_parts, ez_parts = [], [], []
    for p in shards:
        d = np.load(p)
        ep_parts.append(d["flags_ep"])
        z_parts.append(d["flags_z"])
        ez_parts.append(d["flags_ez"])
    flags_ep = np.concatenate(ep_parts, axis=0)
    flags_z = np.concatenate(z_parts, axis=0)
    flags_ez = np.concatenate(ez_parts, axis=0)
    cfg = {}
    cfg_path = flag_dir / "run_config.json"
    if cfg_path.is_file():
        cfg = json.loads(cfg_path.read_text())
    return flags_ep, flags_z, flags_ez, cfg


def _bootstrap_auc(
    flags: np.ndarray,
    y: np.ndarray,
    keep: np.ndarray,
    repeat_ns: list[int],
    n_boot: int,
    rng: np.random.Generator,
    lo_q: float,
    hi_q: float,
) -> pd.DataFrame:
    """flags: (n_rep, n_eval), y/keep length n_eval."""
    n_rep, n_eval = flags.shape
    y_k = y[keep]
    rows = []
    for n in repeat_ns:
        if n < 1:
            continue
        aucs = np.empty(n_boot, dtype=float)
        for b in range(n_boot):
            if n <= n_rep:
                idx = rng.choice(n_rep, size=n, replace=False)
            else:
                idx = rng.choice(n_rep, size=n, replace=True)
            ratios = flags[idx].mean(axis=0)[keep]
            aucs[b] = roc_auc(ratios, y_k)
        rows.append(
            {
                "repeat_n": n,
                "auc_mean": float(np.nanmean(aucs)),
                "auc_median": float(np.nanmedian(aucs)),
                "auc_lo": float(np.nanquantile(aucs, lo_q)),
                "auc_hi": float(np.nanquantile(aucs, hi_q)),
                "auc_std": float(np.nanstd(aucs)),
                "n_boot": n_boot,
                "n_flags_available": n_rep,
            }
        )
        console.print(
            f"  N={n:<8d} AUC={rows[-1]['auc_mean']:.4f} "
            f"[{rows[-1]['auc_lo']:.4f}, {rows[-1]['auc_hi']:.4f}]"
        )
    return pd.DataFrame(rows)


def _add_curve(
    fig: go.Figure, df: pd.DataFrame, name: str, color: str, fill: str
) -> None:
    fig.add_trace(
        go.Scatter(
            x=df["repeat_n"],
            y=df["auc_hi"],
            mode="lines",
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
            legendgroup=name,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=df["repeat_n"],
            y=df["auc_lo"],
            mode="lines",
            line=dict(width=0),
            fill="tonexty",
            fillcolor=fill,
            name=f"{name} band",
            legendgroup=name,
            showlegend=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=df["repeat_n"],
            y=df["auc_mean"],
            mode="lines+markers",
            name=name,
            legendgroup=name,
            line=dict(color=color, width=2),
            marker=dict(size=6),
            hovertemplate="N=%{x}<br>AUC=%{y:.4f}<extra>" + name + "</extra>",
        )
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--flag-dir", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Dir with flags_*.npz + eval_samples.tsv (…/fixed_flags)")
@click.option("--output-dir", default=None, type=click.Path(path_type=Path))
@click.option("--ff-min", default=0.01, show_default=True, type=float)
@click.option("--repeat-grid", default=DEFAULT_GRID, show_default=True)
@click.option("--n-boot", default=200, show_default=True, type=int)
@click.option("--seed", default=0, show_default=True, type=int)
@click.option("--lo-q", default=0.05, show_default=True, type=float)
@click.option("--hi-q", default=0.95, show_default=True, type=float)
def main(
    flag_dir: Path,
    output_dir: Path | None,
    ff_min: float,
    repeat_grid: str,
    n_boot: int,
    seed: int,
    lo_q: float,
    hi_q: float,
) -> None:
    out = output_dir or (flag_dir.parent / "auc_learning_curve")
    out.mkdir(parents=True, exist_ok=True)

    eval_info = pd.read_csv(flag_dir / "eval_samples.tsv", sep="\t")
    flags_ep, flags_z, flags_ez, cfg = _load_flags(flag_dir)
    n_rep = flags_ep.shape[0]
    if flags_ep.shape[1] != len(eval_info):
        raise click.ClickException(
            f"flag n_eval={flags_ep.shape[1]} != eval_samples {len(eval_info)}"
        )
    console.print(f"Loaded flags: n_rep={n_rep}, n_eval={len(eval_info)}")

    ff = pd.to_numeric(eval_info["ff_before_mq"], errors="coerce").to_numpy()
    labels = eval_info["label"].astype(str)
    y = labels.map(is_trisomy_label).to_numpy()
    keep = (ff >= ff_min) & (labels.eq("Normal").to_numpy() | y)
    console.print(
        f"AUC cohort ff≥{ff_min}: N={int((keep & ~y).sum())} T={int((keep & y).sum())}"
    )

    grid = [int(x) for x in repeat_grid.split(",") if x.strip()]
    grid = [n for n in grid if n <= n_rep]
    if not grid or grid[-1] != n_rep:
        if n_rep not in grid:
            grid.append(n_rep)
    grid = sorted(set(grid))
    rng = np.random.default_rng(seed)

    curves = {}
    for name, flags in [("ezscore", flags_ez), ("episcore", flags_ep), ("zscore", flags_z)]:
        console.rule(f"[cyan]{name}[/cyan]")
        curves[name] = _bootstrap_auc(
            flags, y, keep, grid, n_boot, rng, lo_q, hi_q
        )
        curves[name].to_csv(out / f"auc_vs_repeats_{name}.tsv", sep="\t", index=False, float_format="%.6f")

    style = {
        "episcore": ("rgb(31,119,180)", "rgba(31,119,180,0.20)"),
        "zscore": ("rgb(44,160,44)", "rgba(44,160,44,0.20)"),
        "ezscore": ("rgb(214,39,40)", "rgba(214,39,40,0.25)"),
    }
    fig = go.Figure()
    for name in ("episcore", "zscore", "ezscore"):
        c, f = style[name]
        _add_curve(fig, curves[name], name, c, f)

    fig.update_layout(
        title=(
            f"Fixed-combo AUC vs repeats (bootstrap n={n_boot}, "
            f"{int(lo_q*100)}–{int(hi_q*100)}% band)<br>"
            f"<sup>ff≥{ff_min*100:.0f}% Normal vs Trisomy · "
            f"ep {cfg.get('ep_threshold')}/{cfg.get('ep_recall')} · "
            f"z {cfg.get('z_threshold')}/{cfg.get('z_recall')} · "
            f"ez cutoff={cfg.get('ez_cutoff', 3)}</sup>"
        ),
        xaxis_title="repeat_n",
        yaxis_title="AUC",
        xaxis_type="log",
        template="plotly_white",
        height=520,
        width=900,
        legend=dict(orientation="h", y=-0.2),
        margin=dict(t=100, b=80),
    )
    fig.update_yaxes(range=[0.5, 1.02])
    html = out / "auc_vs_repeats.html"
    fig.write_html(str(html), include_plotlyjs="cdn", full_html=True)

    summary = {
        "n_flags": n_rep,
        "n_eval": len(eval_info),
        "n_boot": n_boot,
        "ff_min": ff_min,
        "repeat_grid": grid,
        "lo_q": lo_q,
        "hi_q": hi_q,
        "final_auc": {k: float(v.iloc[-1]["auc_mean"]) for k, v in curves.items()},
    }
    (out / "auc_curve_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    console.print(f"[green]OK[/green] Wrote {html}")
    console.print(f"  final ezscore AUC={summary['final_auc']['ezscore']:.4f}")


if __name__ == "__main__":
    main()
