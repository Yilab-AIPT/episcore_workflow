#!/usr/bin/env python3
"""Collect ff_before_mq + chrY_ratio for old 240k samples and new allosome cohort."""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
import pysam
from rich.console import Console

import config as cfg

console = Console()


def _chr_ratios(bam_path: str) -> tuple[float | None, float | None]:
    try:
        idxstats = pysam.idxstats(bam_path)
    except Exception as exc:  # noqa: BLE001
        console.print(f"[yellow]idxstats failed[/yellow] {bam_path}: {exc}")
        return None, None
    chrY_count = chrX_count = total = 0
    for line in idxstats.strip().split("\n"):
        if not line.strip():
            continue
        chrom, _length, count, _unmapped = line.split("\t")[:4]
        count = int(count)
        total += count
        if chrom in ("chrY", "Y"):
            chrY_count = count
        if chrom in ("chrX", "X"):
            chrX_count = count
    if total <= 0:
        return None, None
    return chrY_count / total, chrX_count / total


def _load_ff() -> pd.DataFrame:
    frames = []
    for path in (cfg.OLD_FF_20260416, cfg.OLD_FF_20260507):
        if path.is_file():
            df = pd.read_csv(path, sep="\t", usecols=["sample", "ff_before_mq"])
            frames.append(df)
    new_summary = cfg.OUTPUT_DIR / "collect_reports" / "summary_report.tsv"
    if new_summary.is_file():
        frames.append(
            pd.read_csv(new_summary, sep="\t", usecols=["sample", "ff_before_mq"])
        )
    if not frames:
        raise click.ClickException("No summary_report.tsv with ff_before_mq found")
    ff = pd.concat(frames, ignore_index=True).drop_duplicates("sample", keep="last")
    return ff


def _bam_table() -> pd.DataFrame:
    frames = [
        pd.read_csv(cfg.OLD_SAMPLESHEET_20260416),
        pd.read_csv(cfg.OLD_SAMPLESHEET_20260507),
        pd.read_csv(cfg.MQRES),
    ]
    ss = pd.concat(frames, ignore_index=True)
    # one clean_bam per sample
    return (
        ss.groupby("sample", as_index=False)
        .agg(clean_bam=("clean_bam", "first"))
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """Write tables/chry_ff.tsv with sample,label,cohort,ff,chrY/X ratios."""
    labels = pd.read_csv(cfg.COHORT_LABELS)
    ff = _load_ff()
    bams = _bam_table()
    df = labels.merge(ff, on="sample", how="left").merge(bams, on="sample", how="left")

    ratios = []
    for sample, bam in zip(df["sample"], df["clean_bam"]):
        if pd.isna(bam):
            ratios.append((None, None))
            continue
        ratios.append(_chr_ratios(str(bam)))
    df["chrY_ratio"] = [r[0] for r in ratios]
    df["chrX_ratio"] = [r[1] for r in ratios]

    out = cfg.CHRY_FF_TSV
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    n_ok = int(df["chrY_ratio"].notna().sum())
    n_ff = int(df["ff_before_mq"].notna().sum())
    console.print(f"[green]Wrote[/green] {out}")
    console.print(f"  samples={len(df)}  with_ff={n_ff}  with_chrY={n_ok}")
    missing_ff = df.loc[df["ff_before_mq"].isna(), "sample"].tolist()
    if missing_ff:
        console.print(
            f"[yellow]Missing FF (run Nextflow for new samples):[/yellow] "
            f"{', '.join(map(str, missing_ff))}"
        )


if __name__ == "__main__":
    main()
