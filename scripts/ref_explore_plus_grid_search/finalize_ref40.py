#!/usr/bin/env python3
"""
Post-aggregate ref-40 finalization: annotate scores, plot confusion matrices,
and search for the best ezscore reference subset.

Expects ``aggregate_results.py`` outputs under ``--output-base``:
    best_sample_scores.tsv
    best_reference_samples.txt
    best_combo_*.csv

Writes:
    best_sample_scores_annotated.tsv
    confusion_{episcore,zscore,ezscore}.png
    best_ezscore_ref_{n}_samples.txt
    best_sample_scores_recalc_ezscore.tsv
    confusion_ezscore_recalc.png
    finalize_summary.json
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Dict, Optional

import click
import numpy as np
import pandas as pd
from rich.console import Console

from ref40_score_eval import (
    CHR_NUMS,
    SCORE_CUTOFF,
    add_score_labels,
    apply_ezscore_pred_override,
    plot_all_confusion_matrices,
    plot_confusion_for_score,
)

console = Console()


def _mcc_from_pred(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    tp = int((y_true & y_pred).sum())
    tn = int((~y_true & ~y_pred).sum())
    fp = int((~y_true & y_pred).sum())
    fn = int((y_true & ~y_pred).sum())
    denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0


def _recalc_ezscore(combined: np.ndarray, ref_idx: np.ndarray) -> np.ndarray:
    ref_vals = combined[ref_idx]
    with np.errstate(invalid="ignore"):
        mu = np.nanmean(ref_vals, axis=0)
        sd = np.nanstd(ref_vals, axis=0, ddof=0)
    mu = np.where(np.isfinite(mu), mu, 0.0)
    sd_safe = np.where(sd > 0, sd, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (combined - mu) / sd_safe


def search_best_ezscore_ref(
    df: pd.DataFrame,
    n_ezscore_ref: int,
    n_repeats: int,
    seed: int,
) -> tuple[np.ndarray, float, list[str]]:
    """Random search over Normal-sample ezscore reference subsets."""
    ep_cols = [f"episcore_chr{n}" for n in CHR_NUMS]
    z_cols = [f"zscore_chr{n}" for n in CHR_NUMS]
    combined = df[ep_cols].to_numpy(dtype=np.float64) + df[z_cols].to_numpy(dtype=np.float64)
    labels = df["label"].astype(str).to_numpy()
    both_normal = (
        (df["pred_label_episcore"] == "Normal") & (df["pred_label_zscore"] == "Normal")
    ).to_numpy()
    y_true = np.array([s.startswith("T") for s in labels], dtype=bool)
    mask_t15 = np.zeros(len(df), dtype=bool)
    if "ff_before_mq" in df.columns:
        mask_t15 = ((df["label"] == "T15") & (df["ff_before_mq"] < 0.01)).to_numpy()

    normal_idx = np.flatnonzero(labels == "Normal")
    if normal_idx.size < n_ezscore_ref:
        raise click.ClickException(
            f"Need at least {n_ezscore_ref} Normal samples, found {normal_idx.size}"
        )

    console.print(f"Normal pool: {normal_idx.size} samples")
    console.print(
        f"Searching {n_repeats} random {n_ezscore_ref}-sample ezscore references ..."
    )

    rng = np.random.default_rng(seed)
    best_mcc = -np.inf
    best_ref_idx: Optional[np.ndarray] = None
    best_ez: Optional[np.ndarray] = None

    for _ in range(n_repeats):
        ref_idx = rng.choice(normal_idx, size=n_ezscore_ref, replace=False)
        ez = _recalc_ezscore(combined, ref_idx)
        any_pos = (ez > SCORE_CUTOFF).any(axis=1)
        y_pred = any_pos & ~both_normal
        y_pred[mask_t15] = False
        mcc = _mcc_from_pred(y_true, y_pred)
        if mcc > best_mcc:
            best_mcc = mcc
            best_ref_idx = ref_idx.copy()
            best_ez = ez.copy()

    if best_ref_idx is None or best_ez is None:
        raise click.ClickException("Ezscore reference search found no valid draw")

    best_ref_samples = df.iloc[best_ref_idx]["sample"].astype(str).tolist()
    console.print(f"Best ezscore MCC: {best_mcc:.4f}")
    console.print(f"Best {n_ezscore_ref}-sample ezscore reference:")
    console.print(", ".join(best_ref_samples))
    return best_ez, best_mcc, best_ref_samples


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-base", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--meta-csv", default=None, type=click.Path(exists=True, dir_okay=False),
              help="Default: <input-dir>/meta.csv if --input-dir set, else required")
@click.option("--input-dir", default=None, type=click.Path(exists=True, file_okay=False),
              help="Grid-search input dir (for meta.csv when --meta-csv omitted)")
@click.option("--scores-file", default="best_sample_scores.tsv", show_default=True,
              help="Annotated scores file name under output-base")
@click.option("--n-ezscore-ref", default=20, show_default=True, type=int)
@click.option("--n-repeats", default=5000, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--skip-ezscore-search", is_flag=True, default=False,
              help="Only annotate scores and plot initial confusion matrices")
def main(
    output_base: str,
    meta_csv: Optional[str],
    input_dir: Optional[str],
    scores_file: str,
    n_ezscore_ref: int,
    n_repeats: int,
    seed: int,
    skip_ezscore_search: bool,
) -> None:
    """Finalize ref-40 grid-search outputs with labels, plots, and ezscore ref search."""
    out_base = Path(output_base)
    scores_path = out_base / scores_file
    if not scores_path.is_file():
        raise click.ClickException(f"Missing {scores_path}; run aggregate_results.py first")

    if meta_csv is None:
        if input_dir is None:
            raise click.ClickException("Provide --meta-csv or --input-dir")
        meta_csv = str(Path(input_dir) / "meta.csv")
    meta_path = Path(meta_csv)
    if not meta_path.is_file():
        raise click.ClickException(f"Missing meta file: {meta_path}")

    console.rule("[bold blue]Finalize ref-40 results")
    console.print(f"  Output base : {out_base}")
    console.print(f"  Scores      : {scores_path}")
    console.print(f"  Meta        : {meta_path}")

    scores = pd.read_csv(scores_path, sep="\t")
    meta = pd.read_csv(meta_path).drop_duplicates("sample", keep="first")
    if "ff_before_mq" not in meta.columns:
        raise click.ClickException("meta.csv missing column: ff_before_mq")
    merge_cols = ["sample", "ff_before_mq"]
    if "ff_after_mq" in meta.columns:
        merge_cols.append("ff_after_mq")
    scores = scores.merge(meta[merge_cols], on="sample", how="left", suffixes=("", "_meta"))
    if "ff_before_mq_meta" in scores.columns:
        scores["ff_before_mq"] = scores["ff_before_mq"].fillna(scores["ff_before_mq_meta"])
        scores = scores.drop(columns=["ff_before_mq_meta"])

    annotated = add_score_labels(scores)
    annotated_path = out_base / "best_sample_scores_annotated.tsv"
    annotated.to_csv(annotated_path, sep="\t", index=False, float_format="%.6f")
    console.print(f"[green]OK[/green] Wrote {annotated_path}")

    initial_mccs = plot_all_confusion_matrices(annotated, out_base)

    summary: Dict[str, object] = {
        "scores_annotated": str(annotated_path),
        "initial_mcc": initial_mccs,
    }

    if skip_ezscore_search:
        summary_path = out_base / "finalize_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        console.rule("[bold green]Done (ezscore search skipped)")
        return

    best_ez, best_ez_mcc, best_ref_samples = search_best_ezscore_ref(
        annotated, n_ezscore_ref, n_repeats, seed
    )

    recalc = annotated.copy()
    for i, n in enumerate(CHR_NUMS):
        recalc[f"ezscore_chr{n}"] = best_ez[:, i]
    recalc = apply_ezscore_pred_override(recalc)

    ref_path = out_base / f"best_ezscore_ref_{n_ezscore_ref}_samples.txt"
    recalc_path = out_base / "best_sample_scores_recalc_ezscore.tsv"
    ref_path.write_text("\n".join(best_ref_samples) + "\n")
    recalc.to_csv(recalc_path, sep="\t", index=False, float_format="%.6f")
    console.print(f"[green]OK[/green] Wrote {ref_path}")
    console.print(f"[green]OK[/green] Wrote {recalc_path}")

    console.print("\nEzscore match_status counts (recalculated):")
    console.print(recalc["match_status_ezscore"].value_counts().to_string())

    _, recalc_mcc = plot_confusion_for_score(
        recalc,
        "ezscore",
        out_base / "confusion_ezscore_recalc.png",
        pred_col="pred_label_ezscore",
    )

    summary.update(
        {
            "best_ezscore_ref_samples": best_ref_samples,
            "best_ezscore_ref_file": str(ref_path),
            "recalc_scores": str(recalc_path),
            "ezscore_search_mcc": best_ez_mcc,
            "ezscore_recalc_mcc": recalc_mcc,
            "n_ezscore_ref": n_ezscore_ref,
            "n_repeats": n_repeats,
            "seed": seed,
        }
    )
    summary_path = out_base / "finalize_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    console.print(f"[green]OK[/green] Wrote {summary_path}")
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
