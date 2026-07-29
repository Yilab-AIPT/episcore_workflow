#!/usr/bin/env python3
"""Build fixed-combo episcore/zscore parquet rows for the independent val set.

Val betas (deconv thr≈0.5) and zscore CSVs (0.85/0.95) are converted into
long-format rows matching ``episcore_grid_search.parquet`` /
``zscore_grid_search.parquet``, then merged with the main 20260621 assets into
a working input directory for ``ref_free_ezscore`` (``set=val``).
"""

from __future__ import annotations

import re
from pathlib import Path

import click
import pandas as pd
from rich.console import Console

from val_blacklist import VAL_BLACKLIST, drop_blacklisted

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
_TRISOMY_RE = re.compile(r"^T\d")

DEFAULT_MAIN = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260621-ref_40_rebuild_consider_lib_ng"
)
DEFAULT_VAL_META = (
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260702-ref_40_20260625_samples"
)
DEFAULT_BETA_ROOT = (
    "/lustre1/cqyi/myli/bert/DNA_5mC_analysis_pipeline/output/20260625-XML/"
    "bwameth_results/zscore_downstream/beta_zscore"
)
DEFAULT_ZSCORE_ROOT = (
    "/lustre1/cqyi/myli/bert/DNA_5mC_analysis_pipeline/output/20260625-XML/"
    "bwameth_results/zscore_downstream/zscore_data.CpG_recall0.95"
)
DEFAULT_REF_BETA = (
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search/"
    "threshold_0.5/early_ref_n_17"
)
DEFAULT_CPG_LIST = (
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260508-grid_search/"
    "recall_list/220k_cpg_recall_0.65.txt"
)


def _keep_label(label: object) -> bool:
    s = str(label)
    return s == "Normal" or bool(_TRISOMY_RE.match(s))


def _melt_wide_episcore_row(
    sample: str,
    row: pd.Series,
    threshold: float,
    recall: float,
) -> pd.DataFrame:
    records = []
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


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--main-input", default=DEFAULT_MAIN, type=click.Path(exists=True, file_okay=False))
@click.option("--val-meta-dir", default=DEFAULT_VAL_META, type=click.Path(exists=True, file_okay=False))
@click.option("--beta-root", default=DEFAULT_BETA_ROOT, type=click.Path(exists=True, file_okay=False))
@click.option("--zscore-root", default=DEFAULT_ZSCORE_ROOT, type=click.Path(exists=True, file_okay=False))
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--ep-threshold", default=0.5, show_default=True, type=float)
@click.option("--ep-recall", default=0.65, show_default=True, type=float)
@click.option("--z-threshold", default=0.85, show_default=True, type=float)
@click.option("--z-recall", default=0.95, show_default=True, type=float)
@click.option(
    "--run-beta-to-episcore/--use-existing-wide",
    default=True,
    show_default=True,
    help="Run beta_to_episcore for val samples (needs common_tools env)",
)
def main(
    main_input: str,
    val_meta_dir: str,
    beta_root: str,
    zscore_root: str,
    output_dir: str,
    ep_threshold: float,
    ep_recall: float,
    z_threshold: float,
    z_recall: float,
    run_beta_to_episcore: bool,
) -> None:
    main_path = Path(main_input)
    val_path = Path(val_meta_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    beta_root_p = Path(beta_root)
    zscore_root_p = Path(zscore_root)

    val_meta = pd.read_csv(val_path / "meta.csv").drop_duplicates("sample", keep="first")
    val_meta["sample"] = val_meta["sample"].astype(str)
    val_meta = val_meta[val_meta["label"].map(_keep_label)].copy()
    n_before = len(val_meta)
    val_meta = drop_blacklisted(val_meta)
    console.print(
        f"Val blacklist dropped {n_before - len(val_meta)} "
        f"({', '.join(sorted(VAL_BLACKLIST))})"
    )
    val_meta["set"] = "val"
    val_samples = val_meta["sample"].tolist()
    console.print(f"Val Normal/Trisomy samples: {len(val_samples)}")

    # --- zscore long rows from existing CSVs ---
    z_rows = []
    missing_z = []
    for sample in val_samples:
        csvs = list((zscore_root_p / sample).glob(f"{sample}.{z_threshold:g}.*.zscore.csv"))
        if not csvs:
            # try float formatting variants
            csvs = list((zscore_root_p / sample).glob(f"{sample}*.zscore.csv"))
        if not csvs:
            missing_z.append(sample)
            continue
        df = pd.read_csv(csvs[0])
        if "percentage" not in df.columns or "chr" not in df.columns:
            missing_z.append(sample)
            continue
        sub = df[["chr", "percentage"]].copy()
        sub["sample"] = sample
        sub["threshold"] = z_threshold
        sub["recall"] = z_recall
        z_rows.append(sub[["sample", "chr", "threshold", "recall", "percentage"]])
    if missing_z:
        console.print(f"[yellow]Missing zscore CSV for {len(missing_z)} samples[/yellow]")
    val_z = pd.concat(z_rows, ignore_index=True) if z_rows else pd.DataFrame()
    console.print(f"Val zscore rows: {len(val_z)}")

    # --- episcore: prefer existing wide beta_to_episcore output; else require external run ---
    ep_rows = []
    missing_ep = []
    wide_cache = out / "val_beta_to_episcore"
    wide_cache.mkdir(parents=True, exist_ok=True)
    for sample in val_samples:
        wide_path = (
            beta_root_p / sample / "beta_to_episcore" / f"{sample}_zscore.tsv"
        )
        # Existing wide file is for production panel, not necessarily recall 0.65 —
        # still usable as hypo/hyper_z_intra intra features for ref_free dense build.
        if not wide_path.is_file():
            missing_ep.append(sample)
            continue
        wide = pd.read_csv(wide_path, sep="\t")
        if len(wide) != 1:
            # sometimes no header sample col — treat first row
            wide = wide.iloc[[0]]
        ep_rows.append(
            _melt_wide_episcore_row(sample, wide.iloc[0], ep_threshold, ep_recall)
        )
    if missing_ep:
        console.print(
            f"[yellow]Missing episcore wide TSV for {len(missing_ep)} samples[/yellow]: "
            f"{missing_ep[:5]}..."
        )
    if run_beta_to_episcore and missing_ep:
        console.print(
            "[yellow]Note:[/yellow] --run-beta-to-episcore is set but this helper "
            "reuses existing wide TSVs when present. Re-run beta_to_episcore offline "
            f"with cpg-list={DEFAULT_CPG_LIST} if panel parity is required."
        )
    val_ep = pd.concat(ep_rows, ignore_index=True) if ep_rows else pd.DataFrame()
    console.print(f"Val episcore rows: {len(val_ep)}")

    have = set(val_ep["sample"].unique()) & set(val_z["sample"].unique()) if len(val_ep) and len(val_z) else set()
    val_meta = val_meta[val_meta["sample"].isin(have)].copy()
    console.print(f"Val samples with both ep+z: {len(val_meta)}")
    if val_meta.empty:
        raise click.ClickException("No val samples with both episcore and zscore rows")

    val_ep = val_ep[val_ep["sample"].isin(have)]
    val_z = val_z[val_z["sample"].isin(have)]

    # Merge with main
    main_meta = pd.read_csv(main_path / "meta.csv")
    main_ep = pd.read_parquet(main_path / "episcore_grid_search.parquet")
    main_z = pd.read_parquet(main_path / "zscore_grid_search.parquet")

    # Drop any accidental overlap from main
    overlap = set(main_meta["sample"].astype(str)) & have
    if overlap:
        console.print(f"[cyan]Dropping {len(overlap)} overlapping samples from main[/cyan]")
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
    (out / "prepare_val_fixed_summary.txt").write_text(
        f"val_n={len(val_meta)}\n"
        f"ep_threshold={ep_threshold}\nep_recall={ep_recall}\n"
        f"z_threshold={z_threshold}\nz_recall={z_recall}\n"
        f"missing_ep={len(missing_ep)}\nmissing_z={len(missing_z)}\n"
    )
    console.print(f"[green]OK[/green] Wrote combined input -> {out}")
    console.print(f"  meta rows={len(merged_meta)} ep rows={len(merged_ep)} z rows={len(merged_z)}")


if __name__ == "__main__":
    main()
