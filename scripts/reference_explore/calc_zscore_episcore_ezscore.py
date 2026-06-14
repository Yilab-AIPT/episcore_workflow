#!/usr/bin/env python3
"""
Compute per-chromosome zscore, episcore, and ezscore for reference-explore runs.

Inputs (under --input-dir):
    meta.csv        : sample metadata with ``ref_type`` (early_ref / analyze) and ``set``
    beta.csv        : per-sample hypo/hyper beta and z_intra values (episcore inputs)
    percentage.csv  : per-sample chromosome percentage values (traditional zscore)

Reference definitions:
    zscore / episcore : reference samples (``ref_type == early_ref`` by default)
    ezscore           : fixed sample list (chr_stats_reference_samples.txt),
                        using (zscore + episcore)

Outputs (under --output-dir):
    score.tsv                    : analyze samples only
    zscore_reference_matrix.tsv    : chr, mean, sd
    episcore_reference_matrix.tsv  : chr, hypo_mean, hypo_sd, hyper_mean, hyper_sd
    ezscore_reference_matrix.tsv   : chr, mean, sd
    reference_samples.tsv          : optional, when reference list is written
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

DEFAULT_EZSCORE_REF_SAMPLES = Path(
    "/lustre1/cqyi/myli/bert/analysis_nipt/multiomics/chr_stats_reference_samples.txt"
)


def load_ezscore_ref_samples(samples_file: Path) -> List[str]:
    """Load one sample ID per line from the ezscore reference list."""
    if not samples_file.is_file():
        raise FileNotFoundError(f"Missing ezscore reference sample list: {samples_file}")
    samples: List[str] = []
    with samples_file.open() as handle:
        for line in handle:
            sample = line.strip()
            if sample and not sample.startswith("#"):
                if 'HCPT' in sample:
                    sample = sample[0:8]
                samples.append(sample)
    if not samples:
        raise ValueError(f"No samples found in ezscore reference list: {samples_file}")
    return samples


def _parse_chr_spec(spec: str) -> List[str]:
    """Parse '1-22' or 'chr1,chr2' into ['chr1', 'chr2', ...]."""
    spec = spec.strip()
    tokens = [s.strip() for s in spec.split(",") if s.strip()]
    out: List[str] = []
    for token in tokens:
        if "-" in token and not token.startswith("chr"):
            try:
                start, end = token.split("-")
                out.extend([f"chr{i}" for i in range(int(start), int(end) + 1)])
                continue
            except Exception:
                pass
        out.append(token if token.startswith("chr") else f"chr{token}")
    return out


def _stack_per_chr(
    df: pd.DataFrame,
    chr_list: List[str],
    suffix: str,
    *,
    dtype: type = np.float64,
) -> np.ndarray:
    cols = [f"{c}_{suffix}" for c in chr_list]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Columns missing from beta table: {missing[:5]}")
    arr = df.loc[:, cols].to_numpy(copy=False)
    if np.issubdtype(dtype, np.integer):
        arr = np.nan_to_num(arr, nan=0).astype(dtype, copy=False)
    else:
        arr = arr.astype(dtype, copy=False)
    return arr


def build_ezscore_ref_mask(
    merged: pd.DataFrame,
    *,
    samples_file: Path = DEFAULT_EZSCORE_REF_SAMPLES,
) -> np.ndarray:
    """Return mask for ezscore normalization from a fixed sample list."""
    ref_samples = load_ezscore_ref_samples(samples_file)
    sample_col = merged["sample"].astype(str)
    mask = sample_col.isin(ref_samples).to_numpy()
    if mask.sum() == 0:
        raise ValueError(
            f"No samples from {samples_file} found in merged meta/beta table"
        )
    missing = sorted(set(ref_samples) - set(sample_col[mask]))
    if missing:
        console.print(
            f"[yellow]Warning[/yellow] {len(missing)} ezscore-ref samples missing "
            f"from input data (e.g. {missing[:3]})"
        )
    return mask


def _ref_mean_std(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    with np.errstate(invalid="ignore"):
        means = np.nanmean(values, axis=0)
        stds = np.nanstd(values, axis=0, ddof=0)
    means = np.where(np.isfinite(means), means, 0.0)
    return means, stds


def _zscore_from_ref(values: np.ndarray, means: np.ndarray, stds: np.ndarray) -> np.ndarray:
    std_safe = np.where(stds > 0, stds, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (values - means) / std_safe
    return z


def _compute_episcore(
    hypo_z_intra: np.ndarray,
    hyper_z_intra: np.ndarray,
    hypo_counts: np.ndarray,
    hyper_counts: np.ndarray,
    ref_hypo_mean: np.ndarray,
    ref_hypo_std: np.ndarray,
    ref_hyper_mean: np.ndarray,
    ref_hyper_std: np.ndarray,
) -> np.ndarray:
    hypo_std_safe = np.where(ref_hypo_std > 0, ref_hypo_std, np.nan)
    hyper_std_safe = np.where(ref_hyper_std > 0, ref_hyper_std, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        hypo_z_inter = (hypo_z_intra - ref_hypo_mean) / hypo_std_safe
        hyper_z_inter = (hyper_z_intra - ref_hyper_mean) / hyper_std_safe

    w_hypo = np.sqrt(hypo_counts.astype(np.float64))
    w_hyper = np.sqrt(hyper_counts.astype(np.float64))
    total_w = np.sqrt(w_hypo ** 2 + w_hyper ** 2)
    total_w = np.where(total_w > 0, total_w, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        s_inter = (hyper_z_inter * w_hyper - hypo_z_inter * w_hypo) / total_w
    return np.where(np.isnan(s_inter), 0.0, s_inter)


def _load_percentage_matrix(
    percentage_csv: Path,
    samples: List[str],
    chr_list: List[str],
) -> np.ndarray:
    pct = pd.read_csv(percentage_csv, sep="\t", usecols=["sample", "chr", "percentage"])
    pct["sample"] = pct["sample"].astype(str)
    pct["chr"] = pct["chr"].astype(str)

    dup = pct.duplicated(subset=["sample", "chr"], keep=False)
    if dup.any():
        raise ValueError(f"percentage.csv has {int(dup.sum())} duplicate (sample, chr) rows")

    pivot = pct.pivot(index="sample", columns="chr", values="percentage")
    pivot = pivot.reindex(index=samples, columns=chr_list)
    missing_samples = pivot.index[pivot.isna().all(axis=1)].tolist()
    if missing_samples:
        console.print(
            f"[yellow]Warning[/yellow] {len(missing_samples)} samples have no "
            f"percentage rows (e.g. {missing_samples[:3]})"
        )
    return pivot.to_numpy(dtype=np.float64)


def _build_score_dataframe(
    samples: Sequence[str],
    chr_list: List[str],
    zscore: np.ndarray,
    episcore: np.ndarray,
    ezscore: np.ndarray,
    row_idx: np.ndarray,
) -> pd.DataFrame:
    data = {"sample": list(samples)}
    for i, chr_name in enumerate(chr_list):
        num = chr_name.removeprefix("chr")
        data[f"zscore_chr{num}"] = zscore[row_idx, i]
        data[f"episcore_chr{num}"] = episcore[row_idx, i]
        data[f"ezscore_chr{num}"] = ezscore[row_idx, i]
    return pd.DataFrame(data)


def load_score_inputs(input_dir: Path) -> Tuple[pd.DataFrame, Path, List[str]]:
    """Load meta+beta merged table, percentage path, and chromosome list."""
    meta_path = input_dir / "meta.csv"
    beta_path = input_dir / "beta.csv"
    pct_path = input_dir / "percentage.csv"
    for path in (meta_path, beta_path, pct_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing required input: {path}")

    meta = pd.read_csv(meta_path)
    for col in ("sample", "set"):
        if col not in meta.columns:
            raise ValueError(f"meta.csv missing column: {col}")
    meta = meta.drop_duplicates(subset="sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)

    beta = pd.read_csv(beta_path)
    if "sample" not in beta.columns:
        raise ValueError("beta.csv missing column: sample")
    beta = beta.drop_duplicates(subset="sample", keep="first")
    beta["sample"] = beta["sample"].astype(str)

    merged = meta.merge(beta, on="sample", how="inner", validate="one_to_one")
    if len(merged) != len(meta):
        missing = set(meta["sample"]) - set(merged["sample"])
        console.print(
            f"[yellow]Warning[/yellow] {len(missing)} meta samples missing from "
            f"beta.csv (e.g. {list(missing)[:3]})"
        )
    return merged, pct_path, _parse_chr_spec("1-22")


def compute_scores(
    merged: pd.DataFrame,
    pct_path: Path,
    chr_list: List[str],
    ref_mask: np.ndarray,
    analyze_mask: np.ndarray,
    ezscore_ref_mask: np.ndarray,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute score and reference matrices for given ref/analyze partition."""
    if ref_mask.sum() == 0:
        raise ValueError("Reference mask selects zero samples")
    if analyze_mask.sum() == 0:
        raise ValueError("Analyze mask selects zero samples")

    sample_order = merged["sample"].tolist()
    hypo_z = _stack_per_chr(merged, chr_list, "hypo_z_intra")
    hyper_z = _stack_per_chr(merged, chr_list, "hyper_z_intra")
    hypo_counts = _stack_per_chr(merged, chr_list, "hypo_cpgs_count", dtype=np.int64)
    hyper_counts = _stack_per_chr(merged, chr_list, "hyper_cpgs_count", dtype=np.int64)

    ref_hypo_mean, ref_hypo_std = _ref_mean_std(hypo_z[ref_mask])
    ref_hyper_mean, ref_hyper_std = _ref_mean_std(hyper_z[ref_mask])

    episcore = _compute_episcore(
        hypo_z,
        hyper_z,
        hypo_counts,
        hyper_counts,
        ref_hypo_mean,
        ref_hypo_std,
        ref_hyper_mean,
        ref_hyper_std,
    )

    pct_values = _load_percentage_matrix(pct_path, sample_order, chr_list)
    pct_ref_mean, pct_ref_std = _ref_mean_std(pct_values[ref_mask])
    zscore = _zscore_from_ref(pct_values, pct_ref_mean, pct_ref_std)

    combined = zscore + episcore
    if ezscore_ref_mask.sum() == 0:
        raise ValueError("ezscore reference mask selects zero samples")
    ez_mean, ez_std = _ref_mean_std(combined[ezscore_ref_mask])
    ezscore = _zscore_from_ref(combined, ez_mean, ez_std)

    episcore_ref_df = pd.DataFrame(
        {
            "chr": chr_list,
            "hypo_mean": ref_hypo_mean,
            "hypo_sd": ref_hypo_std,
            "hyper_mean": ref_hyper_mean,
            "hyper_sd": ref_hyper_std,
        }
    )
    zscore_ref_df = pd.DataFrame({"chr": chr_list, "mean": pct_ref_mean, "sd": pct_ref_std})
    ezscore_ref_df = pd.DataFrame({"chr": chr_list, "mean": ez_mean, "sd": ez_std})

    analyze_samples = merged.loc[analyze_mask, "sample"].tolist()
    analyze_idx = np.flatnonzero(analyze_mask)
    score_df = _build_score_dataframe(
        analyze_samples,
        chr_list,
        zscore,
        episcore,
        ezscore,
        analyze_idx,
    )
    return score_df, zscore_ref_df, episcore_ref_df, ezscore_ref_df


def write_score_outputs(
    out_dir: Path,
    score_df: pd.DataFrame,
    zscore_ref_df: pd.DataFrame,
    episcore_ref_df: pd.DataFrame,
    ezscore_ref_df: pd.DataFrame,
    reference_samples: Optional[Sequence[str]] = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    score_df.to_csv(out_dir / "score.tsv", sep="\t", index=False, float_format="%.6f")
    zscore_ref_df.to_csv(out_dir / "zscore_reference_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    episcore_ref_df.to_csv(out_dir / "episcore_reference_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    ezscore_ref_df.to_csv(out_dir / "ezscore_reference_matrix.tsv", sep="\t", index=False, float_format="%.6f")
    if reference_samples is not None:
        pd.DataFrame({"sample": list(reference_samples)}).to_csv(
            out_dir / "reference_samples.tsv",
            sep="\t",
            index=False,
        )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing meta.csv, beta.csv, percentage.csv",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Directory for score.tsv and reference matrices",
)
@click.option(
    "--chr-spec",
    default="1-22",
    show_default=True,
    help='Chromosomes to analyze (e.g. "1-22").',
)
@click.option(
    "--ezscore-ref-samples",
    default=str(DEFAULT_EZSCORE_REF_SAMPLES),
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Sample list (one ID per line) for ezscore mean/sd reference",
)
def main(input_dir: str, output_dir: str, chr_spec: str, ezscore_ref_samples: str) -> None:
    """Compute zscore / episcore / ezscore for analyze samples."""
    input_path = Path(input_dir)
    out_path = Path(output_dir)
    chr_list = _parse_chr_spec(chr_spec)
    console.rule("[bold blue]zscore / episcore / ezscore calculation")

    merged, pct_path, _ = load_score_inputs(input_path)
    if "ref_type" not in merged.columns:
        raise click.ClickException("meta.csv missing column: ref_type")

    ref_mask = (merged["ref_type"].astype(str) == "early_ref").to_numpy()
    analyze_mask = (merged["ref_type"].astype(str) == "analyze").to_numpy()
    if "label" not in merged.columns:
        raise click.ClickException("meta.csv missing column: label")
    ezscore_ref_mask = build_ezscore_ref_mask(merged, samples_file=Path(ezscore_ref_samples))

    if ref_mask.sum() == 0:
        raise click.ClickException("No early_ref samples found in meta.csv")
    if analyze_mask.sum() == 0:
        raise click.ClickException("No analyze samples found in meta.csv")

    console.print(f"  Samples total   : {len(merged)}")
    console.print(f"  early_ref       : {int(ref_mask.sum())}")
    console.print(f"  analyze         : {int(analyze_mask.sum())}")
    console.print(
        f"  ezscore ref     : {int(ezscore_ref_mask.sum())} "
        f"(from {Path(ezscore_ref_samples).name})"
    )
    console.print(f"  Chromosomes     : {len(chr_list)} ({chr_list[0]}..{chr_list[-1]})")

    outputs = compute_scores(merged, pct_path, chr_list, ref_mask, analyze_mask, ezscore_ref_mask)
    write_score_outputs(
        out_path,
        *outputs,
        reference_samples=merged.loc[ref_mask, "sample"].tolist(),
    )

    console.print(f"\n[green]OK[/green] Wrote {out_path / 'score.tsv'} ({len(outputs[0])} analyze samples)")
    console.print(f"[green]OK[/green] Wrote reference matrices under {out_path}")
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
