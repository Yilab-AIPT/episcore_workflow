#!/usr/bin/env python3
"""Merge val filtered-grid episcore/zscore rows into a combined input dir."""

from __future__ import annotations

import re
from pathlib import Path

import click
import pandas as pd
from rich.console import Console

from val_blacklist import drop_blacklisted

console = Console()
CHR_LIST = [f"chr{i}" for i in range(1, 23)]
_TRISOMY_RE = re.compile(r"^T\d")


def _keep_label(label: object) -> bool:
    s = str(label)
    return s == "Normal" or bool(_TRISOMY_RE.match(s))


def _melt_wide(path: Path, threshold: float, recall: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    # wide file: one row per sample, or sample as index
    if "sample" not in df.columns:
        # first column may be sample id
        df = df.reset_index()
        if "index" in df.columns:
            df = df.rename(columns={"index": "sample"})
    records = []
    for _, row in df.iterrows():
        sample = str(row.get("sample", row.name))
        for chr_name in CHR_LIST:
            records.append(
                {
                    "sample": sample,
                    "chr": chr_name,
                    "threshold": threshold,
                    "recall": recall,
                    "hypo_z_intra": float(row[f"{chr_name}_hypo_z_intra"]),
                    "hyper_z_intra": float(row[f"{chr_name}_hyper_z_intra"]),
                    "hypo_cpgs_count": float(row[f"{chr_name}_hypo_cpgs_count"]),
                    "hyper_cpgs_count": float(row[f"{chr_name}_hyper_cpgs_count"]),
                }
            )
    return pd.DataFrame.from_records(records)


@click.command()
@click.option("--main-input", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--val-meta-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--val-zscore-root", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--val-episcore-root", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
def main(
    main_input: str,
    val_meta_dir: str,
    val_zscore_root: str,
    val_episcore_root: str,
    output_dir: str,
) -> None:
    main_path = Path(main_input)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    val_meta = pd.read_csv(Path(val_meta_dir) / "meta.csv").drop_duplicates("sample")
    val_meta["sample"] = val_meta["sample"].astype(str)
    val_meta = val_meta[val_meta["label"].map(_keep_label)]
    val_meta = drop_blacklisted(val_meta)
    val_meta["set"] = "val"
    val_samples = set(val_meta["sample"])

    # --- zscore ---
    z_rows = []
    z_root = Path(val_zscore_root) / "zscore_results"
    for combo_dir in sorted(z_root.glob("recall.*_cutoff.*")):
        m = re.match(r"recall\.(.+)_cutoff\.(.+)$", combo_dir.name)
        if not m:
            continue
        rec, thr = float(m.group(1)), float(m.group(2))
        for sample in val_samples:
            csvs = list(combo_dir.glob(f"{sample}.{thr:g}.*.zscore.csv"))
            if not csvs:
                csvs = list(combo_dir.glob(f"{sample}*.zscore.csv"))
            if not csvs:
                continue
            df = pd.read_csv(csvs[0])
            sub = df[["chr", "percentage"]].copy()
            sub["sample"] = sample
            sub["threshold"] = thr
            sub["recall"] = rec
            z_rows.append(sub[["sample", "chr", "threshold", "recall", "percentage"]])
    val_z = pd.concat(z_rows, ignore_index=True) if z_rows else pd.DataFrame()
    console.print(f"Val zscore rows: {len(val_z)} samples={val_z['sample'].nunique() if len(val_z) else 0}")

    # --- episcore from beta_to_episcore outputs ---
    ep_rows = []
    ep_root = Path(val_episcore_root)
    for combo_dir in sorted(ep_root.glob("threshold_*_recall_*")):
        m = re.match(r"threshold_(.+)_recall_(.+)$", combo_dir.name)
        if not m:
            continue
        thr, rec = float(m.group(1)), float(m.group(2))
        wide = combo_dir / "_analyze_zscore.tsv.gz"
        if not wide.is_file():
            console.print(f"[yellow]missing[/yellow] {wide}")
            continue
        ep_rows.append(_melt_wide(wide, thr, rec))
    val_ep = pd.concat(ep_rows, ignore_index=True) if ep_rows else pd.DataFrame()
    if len(val_ep):
        val_ep = val_ep[val_ep["sample"].isin(val_samples)]
    console.print(
        f"Val episcore rows: {len(val_ep)} samples={val_ep['sample'].nunique() if len(val_ep) else 0}"
    )

    have = set()
    if len(val_ep) and len(val_z):
        # require full filtered combo coverage later; for merge keep intersection
        have = set(val_ep["sample"].unique()) & set(val_z["sample"].unique())
    val_meta = val_meta[val_meta["sample"].isin(have)].copy()
    val_ep = val_ep[val_ep["sample"].isin(have)] if len(val_ep) else val_ep
    val_z = val_z[val_z["sample"].isin(have)] if len(val_z) else val_z
    console.print(f"Val samples with both: {len(val_meta)}")

    main_meta = pd.read_csv(main_path / "meta.csv")
    main_ep = pd.read_parquet(main_path / "episcore_grid_search.parquet")
    main_z = pd.read_parquet(main_path / "zscore_grid_search.parquet")
    overlap = set(main_meta["sample"].astype(str)) & have
    if overlap:
        main_meta = main_meta[~main_meta["sample"].astype(str).isin(overlap)]
        main_ep = main_ep[~main_ep["sample"].astype(str).isin(overlap)]
        main_z = main_z[~main_z["sample"].astype(str).isin(overlap)]

    merged_meta = pd.concat([main_meta, val_meta], ignore_index=True, sort=False)
    merged_ep = pd.concat([main_ep, val_ep], ignore_index=True)
    merged_z = pd.concat([main_z, val_z], ignore_index=True)

    merged_meta.to_csv(out / "meta.csv", index=False)
    merged_ep.to_parquet(out / "episcore_grid_search.parquet", index=False, compression="snappy")
    merged_z.to_parquet(out / "zscore_grid_search.parquet", index=False, compression="snappy")
    val_meta.to_csv(out / "val_samples.tsv", sep="\t", index=False)
    console.print(f"[green]OK[/green] Wrote {out}")


if __name__ == "__main__":
    main()
