#!/usr/bin/env python3
"""
Aggregate random ref-40 grid-search repeats and pick the best reference.

Scans ``<output-base>/randomly_select_ref_40/repeat_*/metrics.tsv`` and the
matching reference / combo / score files, then:

    1. Builds ``summary_all_repeats.tsv`` with one row per repeat containing the
       reference sample list and episcore/zscore/ezscore MCC/TP/TN/FP/FN.
       Supports ``dev_test_split`` (dev + test sets) and ``all`` (single all set).
    2. Selects the best repeat by a configurable score+metric criterion.
    3. Writes the best reference list, best combos, best metrics, and the
       per-sample scores under that best reference + combo to the output base.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import click
import pandas as pd
from rich.console import Console

console = Console()

SCORES = ["episcore", "zscore", "ezscore"]
METRICS = ["mcc", "tp", "tn", "fp", "fn"]


def _detect_sets(metrics: pd.DataFrame) -> List[str]:
    sets = metrics["set"].astype(str).unique().tolist()
    if "all" in sets:
        return ["all"]
    if "dev" in sets and "test" in sets:
        return ["dev", "test"]
    raise ValueError(f"Unexpected set values in metrics.tsv: {sets}")


def _load_repeat(repeat_dir: Path, eval_sets: Optional[Sequence[str]] = None) -> Optional[dict]:
    metrics_path = repeat_dir / "metrics.tsv"
    ref_path = repeat_dir / "reference_samples.txt"
    if not metrics_path.is_file() or not ref_path.is_file():
        return None
    metrics = pd.read_csv(metrics_path, sep="\t")
    if eval_sets is None:
        eval_sets = _detect_sets(metrics)

    try:
        repeat_index = int(repeat_dir.name.split("_")[-1])
    except ValueError:
        repeat_index = int(metrics["repeat_index"].iloc[0])

    row: dict = {"repeat_index": repeat_index}
    ref_samples = [s.strip() for s in ref_path.read_text().splitlines() if s.strip()]
    row["reference_list"] = ",".join(ref_samples)

    lut = {(r["score"], r["set"]): r for _, r in metrics.iterrows()}
    for score in SCORES:
        for st in eval_sets:
            r = lut.get((score, st))
            if r is None:
                return None
            for m in METRICS:
                row[f"{score}_{st}_{m}"] = r[m]
    return row


def _selection_key(df: pd.DataFrame, score: str, metric: str, eval_sets: Sequence[str]) -> pd.Series:
    """Return the ranking key (higher == better) for each repeat row."""
    if eval_sets == ["all"]:
        if metric in ("all", "mean_dev_test", "dev", "test", "min_dev_test"):
            return df[f"{score}_all_mcc"]
        raise click.ClickException(f"Unknown --select-metric for mode=all: {metric}")

    dev = df[f"{score}_dev_mcc"]
    test = df[f"{score}_test_mcc"]
    if metric == "mean_dev_test":
        return (dev + test) / 2.0
    if metric == "dev":
        return dev
    if metric == "test":
        return test
    if metric == "min_dev_test":
        return pd.concat([dev, test], axis=1).min(axis=1)
    raise click.ClickException(f"Unknown --select-metric: {metric}")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-base", required=True, type=click.Path(exists=True, file_okay=False),
              help="Output base that contains randomly_select_ref_40/repeat_*/")
@click.option("--select-score", default="ezscore", show_default=True,
              type=click.Choice(SCORES), help="Score whose MCC drives reference selection")
@click.option("--select-metric", default=None,
              help="MCC aggregation: dev_test_split -> mean_dev_test|dev|test|min_dev_test; all -> all")
def main(output_base: str, select_score: str, select_metric: Optional[str]) -> None:
    """Aggregate repeats and emit the best reference + metrics + scores."""
    out_base = Path(output_base)
    repeats_root = out_base / "randomly_select_ref_40"
    if not repeats_root.is_dir():
        raise click.ClickException(f"Missing {repeats_root}")

    repeat_dirs = sorted(
        (d for d in repeats_root.glob("repeat_*") if d.is_dir()),
        key=lambda d: int(d.name.split("_")[-1]),
    )
    if not repeat_dirs:
        raise click.ClickException(f"No repeat_* directories under {repeats_root}")

    # Detect eval sets from the first complete repeat.
    eval_sets: Optional[List[str]] = None
    for d in repeat_dirs:
        metrics_path = d / "metrics.tsv"
        if metrics_path.is_file():
            eval_sets = _detect_sets(pd.read_csv(metrics_path, sep="\t"))
            break
    if eval_sets is None:
        raise click.ClickException("Could not detect eval sets from repeat metrics")

    if select_metric is None:
        select_metric = "all" if eval_sets == ["all"] else "mean_dev_test"

    console.rule("[bold blue]Aggregate ref-40 grid search")
    console.print(f"  Output base    : {out_base}")
    console.print(f"  Repeat dirs    : {len(repeat_dirs)}")
    console.print(f"  Eval sets      : {', '.join(eval_sets)}")
    console.print(f"  Selection      : {select_score} / {select_metric}")

    rows: List[dict] = []
    skipped: List[str] = []
    for d in repeat_dirs:
        row = _load_repeat(d, eval_sets)
        if row is None:
            skipped.append(d.name)
            continue
        rows.append(row)

    if not rows:
        raise click.ClickException("No complete repeats found to aggregate")
    if skipped:
        console.print(f"[yellow]Warning[/yellow] skipped {len(skipped)} incomplete repeats: {skipped[:5]}")

    summary = pd.DataFrame(rows).sort_values("repeat_index").reset_index(drop=True)
    summary["select_key"] = _selection_key(summary, select_score, select_metric, eval_sets)

    metric_cols = [f"{s}_{st}_{m}" for s in SCORES for st in eval_sets for m in METRICS]
    summary = summary[["repeat_index", "select_key"] + metric_cols + ["reference_list"]]

    summary_path = out_base / "summary_all_repeats.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.6f")
    console.print(f"[green]OK[/green] Wrote {summary_path} ({len(summary)} repeats)")

    ranked = summary.sort_values(
        ["select_key", "repeat_index"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    best = ranked.iloc[0]
    best_index = int(best["repeat_index"])
    best_dir = repeats_root / f"repeat_{best_index}"

    console.print("\n[bold]Top 5 repeats by selection key[/bold]")
    show_cols = ["repeat_index", "select_key"]
    for s in SCORES:
        for st in eval_sets:
            show_cols.append(f"{s}_{st}_mcc")
    console.print(ranked[show_cols].head(5).to_string(index=False))

    shutil.copyfile(best_dir / "reference_samples.txt", out_base / "best_reference_samples.txt")
    shutil.copyfile(best_dir / "best_combo_episcore.csv", out_base / "best_combo_episcore.csv")
    shutil.copyfile(best_dir / "best_combo_zscore.csv", out_base / "best_combo_zscore.csv")
    shutil.copyfile(best_dir / "scores.tsv", out_base / "best_sample_scores.tsv")

    best_metrics = pd.read_csv(best_dir / "metrics.tsv", sep="\t")
    best_metrics.to_csv(out_base / "best_metrics.tsv", sep="\t", index=False, float_format="%.6f")

    best_summary = {
        "best_repeat_index": best_index,
        "select_score": select_score,
        "select_metric": select_metric,
        "select_key": float(best["select_key"]),
        "n_repeats": int(len(summary)),
        "eval_sets": list(eval_sets),
        "metrics": {
            f"{s}_{st}": {
                m: (float(best[f"{s}_{st}_{m}"]) if m == "mcc" else int(best[f"{s}_{st}_{m}"]))
                for m in METRICS
            }
            for s in SCORES for st in eval_sets
        },
    }
    (out_base / "best_summary.json").write_text(json.dumps(best_summary, indent=2))

    console.print(f"\n[bold green]Best reference[/bold green] repeat_{best_index} "
                  f"(select_key={best['select_key']:.4f})")
    for s in SCORES:
        for st in eval_sets:
            console.print(
                f"  {s:9s} {st:5s} MCC={best[f'{s}_{st}_mcc']:.4f} "
                f"(TP={int(best[f'{s}_{st}_tp'])} TN={int(best[f'{s}_{st}_tn'])} "
                f"FP={int(best[f'{s}_{st}_fp'])} FN={int(best[f'{s}_{st}_fn'])})"
            )
    console.print(f"\n[green]OK[/green] Wrote best_* outputs under {out_base}")
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
