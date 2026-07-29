#!/usr/bin/env python3
"""Build samplesheet + task TSV for filtered-range zscore cal on val samples."""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
from rich.console import Console

from val_blacklist import VAL_BLACKLIST, drop_blacklisted

console = Console()

DEFAULT_VAL_META = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260702-ref_40_20260625_samples"
)
Z_THRESHOLDS = [0.8, 0.85, 0.9, 0.95]
Z_RECALLS = [round(0.9 + i * 0.01, 2) for i in range(10)]  # 0.9 .. 0.99


@click.command()
@click.option("--val-meta-dir", default=DEFAULT_VAL_META, type=click.Path(exists=True))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(val_meta_dir: str, output_dir: str) -> None:
    meta_dir = Path(val_meta_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(meta_dir / "meta.csv").drop_duplicates("sample")
    meta["sample"] = meta["sample"].astype(str)
    keep = meta["label"].astype(str).eq("Normal") | meta["label"].astype(str).str.match(
        r"^T\d"
    )
    meta = drop_blacklisted(meta.loc[keep])
    samples = set(meta["sample"])

    mq = pd.read_csv(meta_dir / "mqres.csv")
    mq["sample"] = mq["sample"].astype(str)
    mq = mq[mq["sample"].isin(samples)].copy()
    # fill_missing_zscore.slurm awk: col1=sample, col3=deconv_res path
    ss = mq[["sample", "deconv_res"]].drop_duplicates().copy()
    ss.insert(1, "tissue", "plasma")
    ss_path = out / "samplesheet.csv"
    ss.to_csv(ss_path, index=False)

    rows = []
    for sample in sorted(samples):
        for thr in Z_THRESHOLDS:
            for rec in Z_RECALLS:
                rows.append(
                    {
                        "score_type": "zscore",
                        "sample": sample,
                        "threshold": thr,
                        "recall": rec,
                    }
                )
    tasks = pd.DataFrame(rows)
    tasks_path = out / "zscore_tasks.tsv"
    tasks.to_csv(tasks_path, sep="\t", index=False)
    console.print(f"[green]OK[/green] samples={len(samples)} tasks={len(tasks)}")
    console.print(f"  blacklist: {sorted(VAL_BLACKLIST)}")
    console.print(f"  -> {ss_path}")
    console.print(f"  -> {tasks_path}")


if __name__ == "__main__":
    main()
