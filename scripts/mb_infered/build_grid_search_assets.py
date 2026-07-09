#!/usr/bin/env python3
"""Build aipt_ref_40 grid-search assets for MB inferred scoring.

Uses the episcore reference matrix sample list to pull matching per-sample
zscore CSVs from a zscore grid-search result directory, then writes the
``best_*`` files expected by ``params.grid_search_result`` (without ezscore
assets).
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]


def _sample_mean_std(values: np.ndarray) -> tuple[float, float]:
    with np.errstate(invalid="ignore"):
        mean = float(np.nanmean(values))
        std = float(np.nanstd(values, ddof=0))
    if not np.isfinite(mean):
        mean = 0.0
    if not np.isfinite(std):
        std = 0.0
    return mean, std


def _write_fixed_combo(
    path: Path,
    *,
    threshold: float,
    recall: float,
    min_recall: float,
    has_target: bool,
) -> None:
    rows = [
        {
            "chr": chrom,
            "threshold": threshold,
            "recall": recall,
            "has_target": has_target,
            "min_recall": min_recall,
        }
        for chrom in CHR_LIST
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _load_episcore_ref_stats(ep_matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for chrom in CHR_LIST:
        num = chrom.replace("chr", "")
        hypo_col = f"chr{num}_hypo_z_intra"
        hyper_col = f"chr{num}_hyper_z_intra"
        for col in (hypo_col, hyper_col):
            if col not in ep_matrix.columns:
                raise ValueError(f"Episcore reference matrix missing column: {col}")
        hypo_mean, hypo_std = _sample_mean_std(ep_matrix[hypo_col].to_numpy(dtype=float))
        hyper_mean, hyper_std = _sample_mean_std(ep_matrix[hyper_col].to_numpy(dtype=float))
        rows.append(
            {
                "chr": chrom,
                "hypo_z_intra_mean": hypo_mean,
                "hypo_z_intra_std": hypo_std,
                "hyper_z_intra_mean": hyper_mean,
                "hyper_z_intra_std": hyper_std,
            }
        )
    return pd.DataFrame(rows)


def _find_zscore_file(zscore_dir: Path, sample: str, zscore_threshold: float) -> Path | None:
    thr_token = f"{zscore_threshold:g}"
    pattern = str(zscore_dir / f"{sample}.{thr_token}.*.zscore.csv")
    matches = sorted(glob.glob(pattern))
    if matches:
        return Path(matches[0])
    for path in sorted(zscore_dir.glob(f"{sample}*.zscore.csv")):
        return path
    return None


def _load_zscore_percentages(
    zscore_dir: Path,
    samples: list[str],
    zscore_threshold: float,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    missing: list[str] = []
    for sample in samples:
        zpath = _find_zscore_file(zscore_dir, sample, zscore_threshold)
        if zpath is None:
            missing.append(sample)
            continue
        df = pd.read_csv(zpath)
        needed = {"chr", "percentage"}
        if not needed.issubset(df.columns):
            raise ValueError(f"{zpath} missing columns {sorted(needed - set(df.columns))}")
        for _, row in df.iterrows():
            records.append(
                {
                    "sample": sample,
                    "chr": str(row["chr"]),
                    "percentage": float(row["percentage"]),
                }
            )
    if missing:
        raise click.ClickException(
            f"No zscore CSV found for {len(missing)} episcore reference sample(s): {missing[:5]}"
        )
    long_df = pd.DataFrame(records)
    long_df = long_df[long_df["chr"].isin(CHR_LIST)]
    return long_df


def _load_zscore_ref_stats(
    zscore_dir: Path,
    samples: list[str],
    zscore_threshold: float,
) -> pd.DataFrame:
    long_df = _load_zscore_percentages(zscore_dir, samples, zscore_threshold)
    rows = []
    for chrom in CHR_LIST:
        sub = long_df.loc[long_df["chr"] == chrom, "percentage"].to_numpy(dtype=float)
        pct_mean, pct_std = _sample_mean_std(sub)
        rows.append({"chr": chrom, "percentage_mean": pct_mean, "percentage_std": pct_std})
    return pd.DataFrame(rows)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--episcore-reference-matrix",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Wide episcore reference matrix (one row per reference sample).",
)
@click.option(
    "--zscore-reference-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Directory with per-sample *.zscore.csv files from MB grid search.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory to write best_combo_* and best_reference_matrix.tsv.",
)
@click.option("--episcore-threshold", default=0.67, show_default=True, type=float)
@click.option("--episcore-recall", default=0.77, show_default=True, type=float)
@click.option("--zscore-threshold", default=0.55, show_default=True, type=float)
@click.option("--zscore-recall", default=0.92, show_default=True, type=float)
@click.option("--has-target", is_flag=True, default=True, show_default=True,
              help="Value for has_target in combo CSVs (unused by Nextflow scoring).")
def main(
    episcore_reference_matrix: str,
    zscore_reference_dir: str,
    output_dir: str,
    episcore_threshold: float,
    episcore_recall: float,
    zscore_threshold: float,
    zscore_recall: float,
    has_target: bool,
) -> None:
    """Build MB fixed-combo grid-search assets without ezscore references."""
    ep_path = Path(episcore_reference_matrix)
    zdir = Path(zscore_reference_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    ep_matrix = pd.read_csv(ep_path, sep="\t")
    if "sample" not in ep_matrix.columns:
        raise click.ClickException(f"Missing 'sample' column in {ep_path}")
    ref_samples = ep_matrix["sample"].astype(str).tolist()
    if not ref_samples:
        raise click.ClickException(f"No reference samples found in {ep_path}")

    console.print(f"Episcore reference samples: {len(ref_samples)}")

    ep_stats = _load_episcore_ref_stats(ep_matrix)
    z_stats = _load_zscore_ref_stats(zdir, ref_samples, zscore_threshold)
    ref_matrix = ep_stats.merge(z_stats, on="chr", how="inner")

    _write_fixed_combo(
        out / "best_combo_episcore.csv",
        threshold=episcore_threshold,
        recall=episcore_recall,
        min_recall=episcore_recall,
        has_target=has_target,
    )
    _write_fixed_combo(
        out / "best_combo_zscore.csv",
        threshold=zscore_threshold,
        recall=zscore_recall,
        min_recall=zscore_recall,
        has_target=has_target,
    )
    ref_matrix.to_csv(out / "best_reference_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    (out / "best_reference_samples.txt").write_text("\n".join(ref_samples) + "\n")

    console.print(f"[green]OK[/green] Wrote {out / 'best_combo_episcore.csv'}")
    console.print(f"[green]OK[/green] Wrote {out / 'best_combo_zscore.csv'}")
    console.print(f"[green]OK[/green] Wrote {out / 'best_reference_matrix.tsv'}")
    console.print(f"[green]OK[/green] Wrote {out / 'best_reference_samples.txt'} ({len(ref_samples)} samples)")
    console.print("[yellow]Note[/yellow] No ezscore assets written; workflow will skip ezscore.")


if __name__ == "__main__":
    try:
        main()
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
