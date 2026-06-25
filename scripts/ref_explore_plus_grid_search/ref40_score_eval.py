#!/usr/bin/env python3
"""Match-status assignment and confusion-matrix plots for ref-40 score tables."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from rich.console import Console  # noqa: E402
from sklearn.metrics import confusion_matrix, matthews_corrcoef  # noqa: E402

console = Console()

CHR_NUMS = list(range(1, 23))
SCORE_CUTOFF = 3.0
STRONG_CUTOFF = 3.0
SCORE_PREFIXES = ("episcore", "zscore", "ezscore")


def assign_pred_label(row: pd.Series, score_prefix: str) -> str:
    """Assign comma-separated trisomy labels from per-chr scores.

    Matches ``update_samplesheet.py`` convention:
      z > 4.5  -> T{n}
      3 < z <= 4.5 -> Gray_T{n}
      otherwise -> Normal
    """
    t_labels, gray_labels = [], []
    for n in CHR_NUMS:
        z = row[f"{score_prefix}_chr{n}"]
        if pd.isna(z):
            continue
        if z > STRONG_CUTOFF:
            t_labels.append(f"T{n}")
        elif z > SCORE_CUTOFF:
            gray_labels.append(f"Gray_T{n}")
    parts = t_labels + gray_labels
    return ",".join(parts) if parts else "Normal"


def assign_match_status(row: pd.Series, score_prefix: str) -> str:
    """Per-sample TP/TN/FP/FN from binary call (any chr score > cutoff)."""
    label = row["label"]
    ff = row.get("ff_before_mq")
    chr_cols = [f"{score_prefix}_chr{n}" for n in CHR_NUMS]
    scores = row[chr_cols]
    all_less3 = (scores < SCORE_CUTOFF).all()
    any_above3 = (scores > SCORE_CUTOFF).any()

    if label == "T15" and pd.notna(ff) and ff < 0.01:
        return "FN"

    if str(label).startswith("T"):
        if any_above3:
            return "TP"
        if all_less3:
            return "FN"
        return "UNK"
    if str(label).lower() == "normal":
        if all_less3:
            return "TN"
        if any_above3:
            return "FP"
        return "UNK"
    return "UNK"


def assign_match_status_from_pred(row: pd.Series, pred_col: str) -> str:
    """TP/TN/FP/FN from a final pred_label column (used for ezscore after override)."""
    label = row["label"]
    pred_label = row[pred_col]
    ff = row.get("ff_before_mq")

    if label == "T15" and pd.notna(ff) and ff < 0.01:
        return "FN"

    is_positive = str(pred_label) != "Normal"
    if str(label).startswith("T"):
        return "TP" if is_positive else "FN"
    if str(label).lower() == "normal":
        return "FP" if is_positive else "TN"
    return "UNK"


def add_score_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Add pred_label_* and match_status_* columns for episcore, zscore, ezscore."""
    out = df.copy()
    for prefix in ("episcore", "zscore"):
        out[f"pred_label_{prefix}"] = out.apply(
            lambda row, p=prefix: assign_pred_label(row, p), axis=1
        )
        out[f"match_status_{prefix}"] = out.apply(
            lambda row, p=prefix: assign_match_status(row, p), axis=1
        )

    out["pred_label_ezscore"] = out.apply(
        lambda row: assign_pred_label(row, "ezscore"), axis=1
    )
    both_normal = (out["pred_label_episcore"] == "Normal") & (out["pred_label_zscore"] == "Normal")
    out.loc[both_normal, "pred_label_ezscore"] = "Normal"
    out["match_status_ezscore"] = out.apply(
        lambda row: assign_match_status_from_pred(row, "pred_label_ezscore"), axis=1
    )
    return out


def apply_ezscore_pred_override(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute ezscore pred/match columns after ezscore values change."""
    out = df.copy()
    out["pred_label_ezscore"] = out.apply(
        lambda row: assign_pred_label(row, "ezscore"), axis=1
    )
    both_normal = (out["pred_label_episcore"] == "Normal") & (out["pred_label_zscore"] == "Normal")
    out.loc[both_normal, "pred_label_ezscore"] = "Normal"
    out["match_status_ezscore"] = out.apply(
        lambda row: assign_match_status_from_pred(row, "pred_label_ezscore"), axis=1
    )
    return out


def plot_confusion_for_score(
    df: pd.DataFrame,
    score_prefix: str,
    output_path: Path,
    pred_col: Optional[str] = None,
) -> Tuple[np.ndarray, float]:
    """Binary Trisomy vs Normal confusion matrix; save heatmap PNG."""
    if pred_col is not None:
        y_pred = (df[pred_col] != "Normal").astype(int).to_numpy()
    else:
        chr_cols = [f"{score_prefix}_chr{n}" for n in CHR_NUMS]
        y_pred = (df[chr_cols] > SCORE_CUTOFF).any(axis=1).astype(int).to_numpy()

    y_true = df["label"].astype(str).str.startswith("T").astype(int).to_numpy()

    if "ff_before_mq" in df.columns:
        mask_t15 = (df["label"] == "T15") & (df["ff_before_mq"] < 0.01)
        y_pred[mask_t15.to_numpy()] = 0

    cm = confusion_matrix(y_true, y_pred, labels=[1, 0])
    mcc = matthews_corrcoef(y_true, y_pred)
    tp, fn = int(cm[0, 0]), int(cm[0, 1])
    fp, tn = int(cm[1, 0]), int(cm[1, 1])

    console.print(f"\n==== {score_prefix} ====")
    console.print(f"MCC: {mcc:.4f}  TP={tp} FN={fn} FP={fp} TN={tn}")

    fig, ax = plt.subplots(figsize=(4, 4))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        ax=ax,
        xticklabels=["Trisomy", "Normal"],
        yticklabels=["Trisomy", "Normal"],
        cbar=False,
        linewidths=0.5,
        linecolor="grey",
        square=True,
    )
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Actual")
    ax.set_title(f"Confusion matrix: {score_prefix}")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    console.print(f"[green]OK[/green] Wrote {output_path}")
    return cm, mcc


def plot_all_confusion_matrices(
    df: pd.DataFrame,
    output_dir: Path,
    suffix: str = "",
) -> Dict[str, float]:
    """Plot episcore, zscore, ezscore confusion matrices under output_dir."""
    tag = f"_{suffix}" if suffix else ""
    mccs: Dict[str, float] = {}
    for prefix in ("episcore", "zscore"):
        _, mcc = plot_confusion_for_score(
            df, prefix, output_dir / f"confusion_{prefix}{tag}.png"
        )
        mccs[prefix] = mcc
    _, mcc = plot_confusion_for_score(
        df,
        "ezscore",
        output_dir / f"confusion_ezscore{tag}.png",
        pred_col="pred_label_ezscore",
    )
    mccs["ezscore"] = mcc
    return mccs


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--scores-tsv", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--meta-csv", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--output-tsv", default="scores_annotated.tsv", show_default=True)
def main(scores_tsv: str, meta_csv: str, output_dir: str, output_tsv: str) -> None:
    """Annotate scores with pred labels / match status and write confusion PNGs."""
    out_dir = Path(output_dir)
    scores = pd.read_csv(scores_tsv, sep="\t")
    meta = pd.read_csv(meta_csv).drop_duplicates("sample", keep="first")
    if "ff_before_mq" not in meta.columns:
        raise click.ClickException("meta.csv missing column: ff_before_mq")
    scores = scores.merge(meta[["sample", "ff_before_mq"]], on="sample", how="left")
    annotated = add_score_labels(scores)
    out_path = out_dir / output_tsv
    out_dir.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    console.print(f"[green]OK[/green] Wrote {out_path}")
    plot_all_confusion_matrices(annotated, out_dir)


if __name__ == "__main__":
    main()
