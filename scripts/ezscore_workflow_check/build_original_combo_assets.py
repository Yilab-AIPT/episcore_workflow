#!/usr/bin/env python3
"""Build aipt_ref_40 grid-search assets matching the original ezscore workflow combo.

Original Snakemake (dna_5mc_analysis_pipeline) fixed combo:
  episcore : threshold=0.5, recall=0.65  (early profile)
  zscore   : threshold=0.85, recall=0.95  (Zrc rule)

Reference sources (same early_ref set used by production):
  episcore : assets/early_reference_beta_zscore.tsv
  zscore   : workflow/zscore/.../igtc_220k.early17.*.reference.csv
  ezscore  : chr_stats.csv (mu/sigma of episcore+zscore over ezscore-ref samples)

Writes under --output-dir:
  best_combo_episcore.csv
  best_combo_zscore.csv
  best_reference_matrix.tsv
  best_reference_samples.txt
  best_ezscore_ref_20_matrix.tsv
  best_ezscore_ref_20_samples.txt
  best_sample_scores_recalc_ezscore.tsv   (stub for plot_scores)
  LOGIC_NOTES.md
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
PROJECT = Path("/lustre1/cqyi/AIPT_2.0/workflow/episcore")
ZSCORE_WF = Path("/lustre1/cqyi/AIPT_2.0/workflow/zscore")

DEFAULT_EP_REF = PROJECT / "assets/early_reference_beta_zscore.tsv"
DEFAULT_Z_REF = (
    ZSCORE_WF
    / "reference/igtc_220k.early17"
    / "igtc_220k.early17.ref_list.0.85.1.0.CpG_final_filtered_recall0.95.NoLen.reference.csv"
)
DEFAULT_EZ_STATS = Path("/lustre1/cqyi/myli/bert/analysis_nipt/multiomics/chr_stats.csv")
DEFAULT_EZ_SAMPLES = Path(
    "/lustre1/cqyi/myli/bert/analysis_nipt/multiomics/chr_stats_reference_samples.txt"
)

EP_THRESHOLD, EP_RECALL = 0.5, 0.65
Z_THRESHOLD, Z_RECALL = 0.85, 0.95


def _norm_hcpt(sample: str) -> str:
    s = str(sample)
    if s.startswith("HCPT") and len(s) > 8:
        return s[:8]
    return s


def _write_fixed_combo(path: Path, threshold: float, recall: float) -> None:
    rows = [
        {
            "chr": chrom,
            "threshold": threshold,
            "recall": recall,
            "has_target": False,
            "min_recall": recall,
        }
        for chrom in CHR_LIST
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _episcore_stats(ep_ref: Path) -> tuple[pd.DataFrame, List[str]]:
    df = pd.read_csv(ep_ref, sep="\t")
    if "sample" not in df.columns:
        raise click.ClickException(f"Missing sample column in {ep_ref}")
    samples = [_norm_hcpt(s) for s in df["sample"].astype(str)]
    rows = []
    for chrom in CHR_LIST:
        num = chrom.removeprefix("chr")
        hypo = df[f"chr{num}_hypo_z_intra"].to_numpy(dtype=float)
        hyper = df[f"chr{num}_hyper_z_intra"].to_numpy(dtype=float)
        # Match beta_to_episcore.py / NIPT: ddof=0
        rows.append(
            {
                "chr": chrom,
                "hypo_z_intra_mean": float(np.nanmean(hypo)),
                "hypo_z_intra_std": float(np.nanstd(hypo, ddof=0)),
                "hyper_z_intra_mean": float(np.nanmean(hyper)),
                "hyper_z_intra_std": float(np.nanstd(hyper, ddof=0)),
            }
        )
    return pd.DataFrame(rows), samples


def _zscore_stats(z_ref: Path) -> pd.DataFrame:
    df = pd.read_csv(z_ref)
    needed = {"chr", "percentage"}
    if not needed.issubset(df.columns):
        raise click.ClickException(f"{z_ref} missing {sorted(needed - set(df.columns))}")
    rows = []
    for chrom in CHR_LIST:
        sub = df.loc[df["chr"].astype(str) == chrom, "percentage"].to_numpy(dtype=float)
        # Match mq_zscore_analyzer.py / polars .std(): ddof=1
        mean = float(np.nanmean(sub))
        std = float(np.nanstd(sub, ddof=1)) if len(sub) > 1 else 0.0
        rows.append({"chr": chrom, "percentage_mean": mean, "percentage_std": std})
    return pd.DataFrame(rows)


def _ezscore_matrix(ez_stats: Path) -> pd.DataFrame:
    df = pd.read_csv(ez_stats)
    # Expected: chr,mu,sigma[,count]
    colmap = {}
    if {"chr", "mu", "sigma"}.issubset(df.columns):
        colmap = {"mu": "mean", "sigma": "std"}
    elif {"chr", "mean", "std"}.issubset(df.columns):
        colmap = {}
    else:
        raise click.ClickException(
            f"{ez_stats} needs columns chr,mu,sigma or chr,mean,std; got {list(df.columns)}"
        )
    out = df.rename(columns=colmap)[["chr", "mean", "std"]].copy()
    out["chr"] = out["chr"].astype(str)
    out = out[out["chr"].isin(CHR_LIST)].set_index("chr").reindex(CHR_LIST).reset_index()
    return out


def _read_sample_list(path: Path) -> List[str]:
    samples = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            samples.append(_norm_hcpt(s))
    return samples


LOGIC_NOTES = """# Logic comparison notes (original vs aipt_ref_40)

## Episcore (fixed combo 0.5 / 0.65)

| Step | NIPT (`nipt.nf` / `CALC_EPISCORE`) | aipt_ref_40 (`AIPT_EPISCORE`) |
|------|-----------------------------------|-------------------------------|
| BAM split | `SPLIT_BAM` → target + background | `SPLIT_BAM_BY_THRESHOLDS` → target per combo thr + **background at params.threshold (0.5)** |
| Methylation | MethylDackel on target **and** background | MethylDackel on all target thr BAMs **and** thr=0.5 background |
| Beta / depth | `extract_beta_value`: beta from target counts; depth filter on **raw_total = target+background** | `FILTER_CPG_BY_DEPTH` on thr=0.5 target+bg bedGraphs → filtered CpG list; `calc_episcore_flexible` beta from target counts intersected with that list |
| CpG panel | `CpG_recall0.65.txt` (same sites as `220k_cpg_recall_0.65.txt`) | `220k_cpg_recall_{recall}.txt` |
| Formula | hypo/hyper z_intra across chr → s_inter vs ref mean/std (ddof=0) | **Same formula** |
| Reference | `early_reference_beta_zscore.tsv` (wide) | `best_reference_matrix.tsv` mean/std (built from the same early matrix here) |

**Verdict:** Math and raw_total depth filtering now match NIPT for the fixed combo. Remaining minor differences: background-only depth-passing sites inflate NIPT CpG counts without target beta contribution; flexible path only keeps sites present in the target bedGraph ∩ filtered list ∩ recall.

## Zscore (fixed combo 0.85 / 0.95)

| Step | `workflow/zscore/mq_zscore_analyzer.py` | aipt_ref_40 (`CALC_ZSCORE_FLEXIBLE`) |
|------|----------------------------------------|--------------------------------------|
| Input | Merged MQ deconv files (dual+single) | Merged `deconv_res` from PREPARE_INPUTS |
| Filters | `prob_class_1 >= cutoff`, `mTcount >= mtcount` | Same |
| CpG filter | pybedtools intersect with recall BED | searchsorted overlap on CpG **start** positions |
| Percentage | `readscount / sum(autosome readscount)` | Same |
| Z-score | `(pct - mean) / std` with polars `.std()` → **ddof=1** | `(pct - mean) / std` from matrix; this build uses **ddof=1** to match |
| Reference | `igtc_220k.early17...reference.csv` | `best_reference_matrix.tsv` percentage_* (from that CSV here) |

**Verdict:** Core filter + percentage + z-score formula match. Overlap implementation differs (BED intersect vs start-containment) but for single-base CpG sites should be equivalent. With this asset build, reference stats match the original early17 file (ddof=1).

## Ezscore

Original: z-normalize `(episcore + zscore)` using `chr_stats.csv` (mu/sigma).
aipt_ref_40: same, via `best_ezscore_ref_20_matrix.tsv` (populated from `chr_stats.csv` here).
"""


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--episcore-reference", default=str(DEFAULT_EP_REF), show_default=True)
@click.option("--zscore-reference", default=str(DEFAULT_Z_REF), show_default=True)
@click.option("--ezscore-stats", default=str(DEFAULT_EZ_STATS), show_default=True)
@click.option("--ezscore-samples", default=str(DEFAULT_EZ_SAMPLES), show_default=True)
def main(
    output_dir: str,
    episcore_reference: str,
    zscore_reference: str,
    ezscore_stats: str,
    ezscore_samples: str,
) -> None:
    """Build original-combo grid-search assets for aipt_ref_40."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    ep_stats, ep_samples = _episcore_stats(Path(episcore_reference))
    z_stats = _zscore_stats(Path(zscore_reference))
    ref_matrix = ep_stats.merge(z_stats, on="chr", how="inner")
    ez_matrix = _ezscore_matrix(Path(ezscore_stats))
    ez_samples = _read_sample_list(Path(ezscore_samples))

    _write_fixed_combo(out / "best_combo_episcore.csv", EP_THRESHOLD, EP_RECALL)
    _write_fixed_combo(out / "best_combo_zscore.csv", Z_THRESHOLD, Z_RECALL)
    ref_matrix.to_csv(out / "best_reference_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    (out / "best_reference_samples.txt").write_text("\n".join(ep_samples) + "\n")
    ez_matrix.to_csv(out / "best_ezscore_ref_20_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    (out / "best_ezscore_ref_20_samples.txt").write_text("\n".join(ez_samples) + "\n")

    # Stub precomputed table so hasEzscoreAssets() is true (plot background optional).
    stub = pd.DataFrame(
        {
            "sample": pd.Series(dtype=str),
            "label": pd.Series(dtype=str),
            "ff_before_mq": pd.Series(dtype=float),
            **{f"ezscore_chr{i}": pd.Series(dtype=float) for i in range(1, 23)},
        }
    )
    stub.to_csv(out / "best_sample_scores_recalc_ezscore.tsv", sep="\t", index=False)
    (out / "LOGIC_NOTES.md").write_text(LOGIC_NOTES)

    console.print(f"[green]OK[/green] Wrote assets under {out}")
    console.print(f"  episcore combo : {EP_THRESHOLD}/{EP_RECALL}  (n_ref={len(ep_samples)})")
    console.print(f"  zscore combo   : {Z_THRESHOLD}/{Z_RECALL}")
    console.print(f"  ezscore ref    : {len(ez_samples)} samples from chr_stats")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
