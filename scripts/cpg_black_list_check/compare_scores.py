#!/usr/bin/env python3
"""Compare per-scheme aipt_ref_40 merge_scores against meta ground truth.

For each scheme under ``--pipeline-root/<scheme>/merge_scores/``, compare
episcore/zscore/ezscore to meta ``beta_zscores`` / ``rc_zscores`` /
``final_zscores``. Also compute delta vs the baseline scheme when present.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
META_MAP = {
    "episcore": "beta_zscores",
    "zscore": "rc_zscores",
    "ezscore": "final_zscores",
}


def _parse_vec(raw) -> np.ndarray:
    return np.asarray([float(x) for x in str(raw).split(",")], dtype=float)


def _load_rebuilt(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    sample = path.name.replace("_scores.tsv", "")
    rows = []
    for score in META_MAP:
        if score not in df.columns:
            continue
        for chrom in CHR_LIST:
            sub = df.loc[df["chr"].astype(str) == chrom, score]
            rows.append(
                {
                    "sample": sample,
                    "score": score,
                    "chr": chrom,
                    "rebuilt": float(sub.iloc[0]) if not sub.empty else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, score), grp in detail.groupby(["sample", "score"], sort=True):
        valid = grp["rebuilt"].notna() & grp["meta"].notna()
        n = int(valid.sum())
        if n == 0:
            rows.append(
                dict(
                    sample=sample,
                    score=score,
                    n_valid_chr=0,
                    max_abs_diff=np.nan,
                    mean_abs_diff=np.nan,
                    rmse=np.nan,
                    pearson_r=np.nan,
                )
            )
            continue
        sub = grp.loc[valid]
        r = (
            float(np.corrcoef(sub["meta"], sub["rebuilt"])[0, 1])
            if n > 1 and sub["meta"].std() > 0 and sub["rebuilt"].std() > 0
            else np.nan
        )
        rows.append(
            dict(
                sample=sample,
                score=score,
                n_valid_chr=n,
                max_abs_diff=sub["abs_diff"].max(),
                mean_abs_diff=sub["abs_diff"].mean(),
                rmse=float(np.sqrt(np.mean(sub["diff"] ** 2))),
                pearson_r=r,
            )
        )
    return pd.DataFrame(rows)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--prepared-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option(
    "--pipeline-root",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Directory containing <scheme>/merge_scores/ outputs.",
)
@click.option("--schemes-file", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(
    input_dir: str,
    prepared_dir: str,
    pipeline_root: str,
    schemes_file: str,
    output_dir: str,
) -> None:
    """Write per-scheme detail/summary TSVs plus an all-scheme aggregate."""
    meta = pd.read_csv(Path(input_dir) / "meta.csv")
    meta["sample"] = meta["sample"].astype(str)
    meta_idx = meta.set_index("sample")

    samples = [
        s.strip()
        for s in (Path(prepared_dir) / "sample_list.txt").read_text().splitlines()
        if s.strip()
    ]
    schemes = [
        s.strip()
        for s in Path(schemes_file).read_text().splitlines()
        if s.strip()
    ]
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    root = Path(pipeline_root)

    all_detail: List[pd.DataFrame] = []
    all_summary: List[pd.DataFrame] = []
    missing_report: Dict[str, List[str]] = {}

    # Prefetch baseline rebuilt scores so scheme order does not matter.
    baseline_rebuild: Dict[str, pd.Series] = {}
    baseline_dir = root / "baseline" / "merge_scores"
    if baseline_dir.is_dir():
        for sample in samples:
            path = baseline_dir / f"{sample}_scores.tsv"
            if path.is_file():
                baseline_rebuild[sample] = _load_rebuilt(path).set_index(
                    ["score", "chr"]
                )["rebuilt"]

    for scheme in schemes:
        score_dir = root / scheme / "merge_scores"
        if not score_dir.is_dir():
            console.print(f"[yellow]SKIP[/yellow] missing {score_dir}")
            missing_report[scheme] = ["__missing_dir__"]
            continue

        detail_rows: List[dict] = []
        missing = []
        for sample in samples:
            path = score_dir / f"{sample}_scores.tsv"
            if not path.is_file():
                missing.append(sample)
                continue
            rebuilt = _load_rebuilt(path)

            mrow = meta_idx.loc[sample]
            for score, mcol in META_MAP.items():
                meta_vec = _parse_vec(mrow[mcol])
                sub = rebuilt[rebuilt["score"] == score].set_index("chr")
                for i, chrom in enumerate(CHR_LIST):
                    r = float(sub.loc[chrom, "rebuilt"]) if chrom in sub.index else np.nan
                    m = meta_vec[i]
                    diff = r - m
                    row = {
                        "scheme": scheme,
                        "sample": sample,
                        "score": score,
                        "chr": chrom,
                        "meta": m,
                        "rebuilt": r,
                        "diff_vs_meta": diff,
                        "abs_diff_vs_meta": abs(diff) if np.isfinite(diff) else np.nan,
                        "baseline": np.nan,
                        "diff_vs_baseline": np.nan,
                        "abs_diff_vs_baseline": np.nan,
                    }
                    if sample in baseline_rebuild:
                        try:
                            b = float(baseline_rebuild[sample].loc[(score, chrom)])
                            row["baseline"] = b
                            row["diff_vs_baseline"] = r - b
                            row["abs_diff_vs_baseline"] = (
                                abs(r - b) if np.isfinite(r - b) else np.nan
                            )
                        except KeyError:
                            pass
                    detail_rows.append(row)

        if missing:
            missing_report[scheme] = missing
            console.print(
                f"[yellow]WARN[/yellow] {scheme}: missing {len(missing)}/{len(samples)} scores"
            )

        if not detail_rows:
            continue

        detail = pd.DataFrame(detail_rows)
        # rename for _summarize compatibility
        detail_for_sum = detail.rename(
            columns={"diff_vs_meta": "diff", "abs_diff_vs_meta": "abs_diff"}
        )
        summary = _summarize(detail_for_sum)
        summary.insert(0, "scheme", scheme)

        scheme_dir = out / scheme
        scheme_dir.mkdir(parents=True, exist_ok=True)
        detail.to_csv(
            scheme_dir / "comparison_detail.tsv",
            sep="\t",
            index=False,
            float_format="%.6f",
        )
        summary.to_csv(
            scheme_dir / "comparison_summary.tsv",
            sep="\t",
            index=False,
            float_format="%.6f",
        )
        all_detail.append(detail)
        all_summary.append(summary)

        med = summary.groupby("score")["rmse"].median()
        console.print(f"[green]OK[/green] {scheme}: median RMSE {med.to_dict()}")

    if all_detail:
        pd.concat(all_detail, ignore_index=True).to_csv(
            out / "all_schemes_detail.tsv", sep="\t", index=False, float_format="%.6f"
        )
    if all_summary:
        pd.concat(all_summary, ignore_index=True).to_csv(
            out / "all_schemes_summary.tsv", sep="\t", index=False, float_format="%.6f"
        )

        # scheme-level rollup
        roll = (
            pd.concat(all_summary, ignore_index=True)
            .groupby(["scheme", "score"], sort=True)
            .agg(
                n_samples=("sample", "nunique"),
                median_rmse=("rmse", "median"),
                mean_rmse=("rmse", "mean"),
                median_max_abs=("max_abs_diff", "median"),
                median_pearson_r=("pearson_r", "median"),
            )
            .reset_index()
        )
        roll.to_csv(out / "scheme_rollup.tsv", sep="\t", index=False, float_format="%.6f")

    (out / "compare_config.json").write_text(
        json.dumps(
            {
                "n_samples": len(samples),
                "schemes": schemes,
                "meta_map": META_MAP,
                "missing": missing_report,
            },
            indent=2,
        )
        + "\n"
    )
    console.print(f"[green]OK[/green] wrote comparisons under {out}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
