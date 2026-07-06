#!/usr/bin/env python3
"""Plot per-chromosome score distributions for a single AIPT ref-40 sample.

Produces ``{output_prefix}_scores.pdf`` with three pages (episcore, zscore,
ezscore). Each page is a 6x4 grid of scatter subplots, one per autosome
(chr1-chr22). In every subplot:

    * x-axis = ``ff_before_mq``, y-axis = the score for that chromosome.
    * Pre-computed reference samples (from ``best_sample_scores_recalc_ezscore``)
      are drawn as round points (alpha 0.8): red when the sample's ``label`` is
      ``T{n}`` and the panel is ``chr{n}``, otherwise gray.
    * The processed single sample is drawn as a star marker.
    * A dotted black horizontal line marks ``y = threshold`` (default 3).
"""

import sys
from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from rich.console import Console  # noqa: E402

console = Console()

CHR_NUMS = list(range(1, 23))
SCORES = ["episcore", "zscore", "ezscore"]
N_ROWS, N_COLS = 6, 4


def melt_precomputed(df: pd.DataFrame, score: str) -> pd.DataFrame:
    """Long-format (sample, label, ff_before_mq, chr_num, value) for one score."""
    records = []
    has_ff = "ff_before_mq" in df.columns
    has_label = "label" in df.columns
    for _, row in df.iterrows():
        ff = row["ff_before_mq"] if has_ff else None
        label = str(row["label"]) if has_label else ""
        for n in CHR_NUMS:
            col = f"{score}_chr{n}"
            if col not in df.columns:
                continue
            records.append(
                {
                    "ff_before_mq": ff,
                    "label": label,
                    "chr_num": n,
                    "value": row[col],
                }
            )
    return pd.DataFrame(records)


def point_colors(sub: pd.DataFrame, n: int) -> list:
    """Red when label == T{n} for the chr{n} panel, else gray."""
    target = f"T{n}"
    return ["red" if lbl == target else "gray" for lbl in sub["label"]]


def plot_score_page(
    pdf: PdfPages,
    score: str,
    precomputed_long: pd.DataFrame,
    sample_df: pd.DataFrame,
    sample_id: str,
    threshold: float,
) -> None:
    """Render one 6x4 page for a single score."""
    fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(20, 24))
    axes = axes.flatten()

    for idx, n in enumerate(CHR_NUMS):
        ax = axes[idx]
        chrom = f"chr{n}"

        sub = precomputed_long[precomputed_long["chr_num"] == n]
        if not sub.empty:
            ax.scatter(
                sub["ff_before_mq"],
                sub["value"],
                c=point_colors(sub, n),
                alpha=0.8,
                s=35,
                edgecolors="none",
                zorder=2,
            )

        srow = sample_df[sample_df["chr"] == chrom]
        if not srow.empty:
            ax.scatter(
                srow["ff_before_mq"],
                srow[score],
                marker="*",
                c="blue",
                s=320,
                edgecolors="black",
                linewidths=0.6,
                zorder=4,
                label=sample_id,
            )

        ax.axhline(threshold, color="black", linestyle=":", linewidth=1.2, zorder=1)
        ax.set_title(chrom, fontsize=12)
        ax.set_xlabel("ff_before_mq", fontsize=9)
        ax.set_ylabel(score, fontsize=9)
        ax.tick_params(labelsize=8)

    # blank the unused panels
    for idx in range(len(CHR_NUMS), len(axes)):
        axes[idx].axis("off")

    legend_handles = [
        Line2D([0], [0], marker="*", color="w", markerfacecolor="blue",
               markeredgecolor="black", markersize=18, label=f"{sample_id} (this sample)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="red",
               markersize=10, label="reference T{n} on chr{n}"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markersize=10, label="other reference samples"),
        Line2D([0], [0], color="black", linestyle=":", label=f"y = {threshold:g}"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, fontsize=12,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(f"{sample_id} — {score}", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))
    pdf.savefig(fig)
    plt.close(fig)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--scores-tsv", required=True, type=click.Path(exists=True),
              help="{sample}_scores.tsv (long: sample,ff_before_mq,chr,episcore,zscore,ezscore).")
@click.option("--precomputed-tsv", required=True, type=click.Path(exists=True),
              help="best_sample_scores_recalc_ezscore.tsv (wide reference scores).")
@click.option("--output-prefix", required=True, type=str,
              help="Output prefix; writes {prefix}_scores.pdf.")
@click.option("--threshold", type=float, default=3.0, show_default=True,
              help="Horizontal reference line (positive-call cutoff).")
def main(
    scores_tsv: str,
    precomputed_tsv: str,
    output_prefix: str,
    threshold: float,
) -> None:
    """Render the 3-page per-chromosome score distribution PDF."""
    try:
        sample_df = pd.read_csv(scores_tsv, sep="\t")
        sample_df["chr"] = sample_df["chr"].astype(str)
        sample_id = str(sample_df["sample"].iloc[0]) if "sample" in sample_df.columns else "sample"
        precomputed = pd.read_csv(precomputed_tsv, sep="\t")

        out_path = f"{output_prefix}_scores.pdf"
        with PdfPages(out_path) as pdf:
            for score in SCORES:
                precomputed_long = melt_precomputed(precomputed, score)
                plot_score_page(pdf, score, precomputed_long, sample_df, sample_id, threshold)

        console.print(f"[green]OK[/green] Wrote {out_path}")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
