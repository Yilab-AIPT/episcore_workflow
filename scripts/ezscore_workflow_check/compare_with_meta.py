#!/usr/bin/env python3
"""Compare aipt_ref_40 merge_scores outputs against meta ground-truth scores."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List

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


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--prepared-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--pipeline-outdir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(
    input_dir: str,
    prepared_dir: str,
    pipeline_outdir: str,
    output_dir: str,
) -> None:
    """Write check_comparison_detail/summary TSVs."""
    meta = pd.read_csv(Path(input_dir) / "meta.csv")
    meta["sample"] = meta["sample"].astype(str)
    meta_idx = meta.set_index("sample")

    check = [
        s.strip()
        for s in (Path(prepared_dir) / "check_samples.txt").read_text().splitlines()
        if s.strip()
    ]
    score_dir = Path(pipeline_outdir) / "merge_scores"
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    detail_rows: List[dict] = []
    missing = []
    for sample in check:
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
                detail_rows.append(
                    {
                        "sample": sample,
                        "score": score,
                        "chr": chrom,
                        "meta": m,
                        "rebuilt": r,
                        "diff": diff,
                        "abs_diff": abs(diff) if np.isfinite(diff) else np.nan,
                    }
                )

    if missing:
        raise click.ClickException(f"Missing merge_scores for: {missing}")

    detail = pd.DataFrame(detail_rows)
    summary_rows = []
    for (sample, score), grp in detail.groupby(["sample", "score"], sort=True):
        valid = grp["rebuilt"].notna() & grp["meta"].notna()
        n = int(valid.sum())
        if n == 0:
            summary_rows.append(
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
        summary_rows.append(
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
    summary = pd.DataFrame(summary_rows)

    detail.to_csv(out / "check_comparison_detail.tsv", sep="\t", index=False, float_format="%.6f")
    summary.to_csv(out / "check_comparison_summary.tsv", sep="\t", index=False, float_format="%.6f")
    (out / "compare_config.json").write_text(
        json.dumps({"check_samples": check, "meta_map": META_MAP}, indent=2) + "\n"
    )

    console.print(f"[green]OK[/green] Compared {len(check)} samples")
    for _, row in summary.sort_values(["sample", "score"]).iterrows():
        console.print(
            f"  {row['sample']:16s} {row['score']:8s} "
            f"max={row['max_abs_diff']:.4f} rmse={row['rmse']:.4f} r={row['pearson_r']:.4f}"
        )


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
