#!/usr/bin/env python3
"""Prepare samplesheets / meta tables for 240k allosome analysis.

Writes under INPUT_DIR:
  samplesheet_nf.csv          – Nextflow split_bam input (new 10 samples, PE+SE)
  episcore_samples_meta.csv   – beta paths; 20260416 females = early_ref
  zscore_samples_meta.csv     – deconv paths; same female early_ref
  cohort_labels.csv           – refreshed from meta.csv + existing old labels

Episcore betas for new samples are filled after Nextflow finishes
(``extract_beta_value/{sample}_beta_value.tsv.gz``).
"""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
from rich.console import Console

import config as cfg

console = Console()


def _beta_path(sample: str) -> str | None:
    for directory in (
        cfg.OLD_BETA_20260416,
        cfg.OLD_BETA_20260507,
        cfg.OUTPUT_DIR / "extract_beta_value",
    ):
        path = directory / f"{sample}_beta_value.tsv.gz"
        if path.is_file():
            return str(path)
    return None


def _load_labels() -> pd.DataFrame:
    """Merge old cohort labels with new meta.csv labels."""
    frames = []
    if cfg.COHORT_LABELS.is_file():
        frames.append(pd.read_csv(cfg.COHORT_LABELS))
    if cfg.META.is_file():
        new = pd.read_csv(cfg.META)
        new["cohort"] = "new"
        frames.append(new[["sample", "label", "cohort"]])
    if not frames:
        raise click.ClickException(
            f"Need {cfg.COHORT_LABELS} and/or {cfg.META} for labels"
        )
    labels = pd.concat(frames, ignore_index=True)
    # Prefer new-cohort rows when sample collides
    labels = labels.drop_duplicates(subset=["sample"], keep="last")
    return labels.reset_index(drop=True)


def _group_deconv(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Collapse PE+SE rows → one row per sample with comma-joined deconv paths."""
    rows = []
    for sample, grp in samplesheet.groupby("sample", sort=False):
        deconv = sorted({str(p) for p in grp["deconv_res"].tolist() if pd.notna(p)})
        bams = sorted({str(p) for p in grp["clean_bam"].tolist() if pd.notna(p)})
        rows.append(
            {
                "sample": sample,
                "clean_bam": bams[0] if bams else None,
                "deconv_paths": ",".join(deconv),
            }
        )
    return pd.DataFrame(rows)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """Build Nextflow samplesheet and episcore/zscore samples_meta tables."""
    cfg.INPUT_DIR.mkdir(parents=True, exist_ok=True)
    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.OUTPUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    cfg.PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    labels = _load_labels()
    labels.to_csv(cfg.COHORT_LABELS, index=False)
    console.print(f"[green]Wrote[/green] {cfg.COHORT_LABELS}  ({len(labels)} samples)")

    # --- Nextflow samplesheet for new 10 (keep PE+SE rows) ---
    mqres = pd.read_csv(cfg.MQRES)
    # Prefer PE (non single_end) first but keep both rows for MERGE_DECONV_RES
    mqres.to_csv(cfg.NF_SAMPLESHEET, index=False)
    console.print(
        f"[green]Wrote[/green] {cfg.NF_SAMPLESHEET}  "
        f"({mqres['sample'].nunique()} unique samples, {len(mqres)} rows)"
    )

    # --- Deconv path tables for old + new ---
    old_ss = pd.concat(
        [pd.read_csv(cfg.OLD_SAMPLESHEET_20260416), pd.read_csv(cfg.OLD_SAMPLESHEET_20260507)],
        ignore_index=True,
    )
    new_ss = mqres.copy()
    all_ss = pd.concat([old_ss, new_ss], ignore_index=True)
    deconv_tbl = _group_deconv(all_ss)
    deconv_tbl = deconv_tbl.merge(labels, on="sample", how="inner")

    # --- Episcore meta: 20260416 females = early_ref; others with beta = analyze ---
    female_20260416 = set(
        labels.loc[
            (labels["label"] == "female") & (labels["cohort"] == "old"),
            "sample",
        ]
    )
    # Restrict early_ref to samples that actually have 20260416 betas
    female_ref = {
        s
        for s in female_20260416
        if (cfg.OLD_BETA_20260416 / f"{s}_beta_value.tsv.gz").is_file()
    }

    epi_rows = []
    for _, row in deconv_tbl.iterrows():
        sample = row["sample"]
        beta = _beta_path(sample)
        if beta is None:
            continue
        # Skip males for chrX episcore curves (keep females + abnormals)
        label = row["label"]
        if pd.isna(label) or label == "male":
            continue
        ref_type = "early_ref" if sample in female_ref else "analyze"
        # Only 20260416 females as reference (user requirement)
        if ref_type == "early_ref" and not beta.startswith(str(cfg.OLD_BETA_20260416)):
            ref_type = "analyze"
        epi_rows.append(
            {
                "sample": sample,
                "beta_path": beta,
                "ref_type": ref_type,
                "label": label,
                "cohort": row["cohort"],
            }
        )
    epi_df = pd.DataFrame(epi_rows)
    if epi_df.empty:
        raise click.ClickException("No episcore samples with beta paths found")
    # Drop analyze females that are also early_ref duplicates — already unique by sample
    n_ref = int((epi_df["ref_type"] == "early_ref").sum())
    n_an = int((epi_df["ref_type"] == "analyze").sum())
    missing_new = [
        s
        for s in labels.loc[labels["cohort"] == "new", "sample"]
        if _beta_path(s) is None
    ]
    epi_df.to_csv(cfg.EPISCORE_SAMPLES_META, index=False)
    console.print(
        f"[green]Wrote[/green] {cfg.EPISCORE_SAMPLES_META}  "
        f"(early_ref={n_ref}, analyze={n_an})"
    )
    if missing_new:
        console.print(
            f"[yellow]New samples missing beta (run Nextflow first):[/yellow] "
            f"{', '.join(missing_new)}"
        )

    # --- Zscore meta: same female early_ref, deconv paths ---
    z_rows = []
    for _, row in deconv_tbl.iterrows():
        sample = row["sample"]
        label = row["label"]
        cohort = row["cohort"]
        if pd.isna(label) or label == "male" or pd.isna(row.get("deconv_paths")):
            continue
        if sample in female_ref:
            ref_type = "early_ref"
        elif label == "female":
            # leftover non-ref females are not plotted on chrX zscore curve
            continue
        else:
            ref_type = "analyze"
        z_rows.append(
            {
                "sample": sample,
                "deconv_paths": row["deconv_paths"],
                "ref_type": ref_type,
                "label": label,
                "cohort": cohort,
            }
        )
    z_df = pd.DataFrame(z_rows)
    z_df.to_csv(cfg.ZSCORE_SAMPLES_META, index=False)
    console.print(
        f"[green]Wrote[/green] {cfg.ZSCORE_SAMPLES_META}  "
        f"(early_ref={int((z_df.ref_type == 'early_ref').sum())}, "
        f"analyze={int((z_df.ref_type == 'analyze').sum())})"
    )
    console.print("[bold green]Done[/bold green]")


if __name__ == "__main__":
    main()
