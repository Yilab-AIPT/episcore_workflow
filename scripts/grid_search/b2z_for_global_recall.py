#!/usr/bin/env python3
"""
Global-recall per-chromosome z-score stitching.

For each chromosome, a "best threshold" (under a single shared recall) is
read from a CSV. Unlike ``b2z_for_best_combo.py`` (which recomputes
``z_intra`` / ``z_inter`` from beta values and counts), this script directly
copies the per-chromosome columns (``chr{n}_*``) out of each combo's existing
``beta_to_episcore.py`` output, because the values are already the ones we want
to keep -- ``recall`` is fixed across chromosomes so the per-sample z_intra
mean/std baseline is consistent.

Inputs:
    --best-combo-csv : CSV with columns ``chr,threshold,recall``. All rows
                        normally share the same ``recall`` (e.g. ``0.64``)
                        but may have different ``threshold`` values.
    --output-base    : Directory containing the
                        ``threshold_{t}_recall_{r}/_analyze_zscore.tsv.gz`` and
                        ``_reference_zscore.tsv.gz`` files produced by
                        ``beta_to_episcore.py``.

Outputs (under ``--output-base/global_recall/`` by default):
    _reference_zscore.tsv.gz
    _analyze_zscore.tsv.gz

The output schemas exactly match those of ``beta_to_episcore.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import click
import pandas as pd
from rich.console import Console

from beta_to_episcore import _write_tsv

console = Console()

# Per-chromosome column suffixes copied straight from each combo's output.
_REFERENCE_SUFFIXES: Tuple[str, ...] = (
    "hypo_beta",
    "hyper_beta",
    "hypo_z_intra",
    "hyper_z_intra",
    "s_intra",
    "hypo_cpgs_count",
    "hyper_cpgs_count",
)
_ANALYZE_SUFFIXES: Tuple[str, ...] = _REFERENCE_SUFFIXES + (
    "hypo_z_inter",
    "hyper_z_inter",
    "s_inter",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_combo_dirname(threshold: float, recall: float) -> str:
    """Match the ``%g``-formatted directory names produced by submit_*.sh."""
    return f"threshold_{threshold:g}_recall_{recall:g}"


def _load_best_combo(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"chr", "threshold", "recall"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"best-combo CSV is missing columns: {sorted(missing)}")
    df = df[["chr", "threshold", "recall"]].copy()
    df["chr"] = df["chr"].astype(str)
    df["threshold"] = df["threshold"].astype(float)
    df["recall"] = df["recall"].astype(float)
    if df["chr"].duplicated().any():
        dups = df["chr"][df["chr"].duplicated()].tolist()
        raise ValueError(f"best-combo CSV has duplicate chr rows: {dups}")
    return df


def _load_combo_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing combo output: {path}")
    return pd.read_csv(path, sep="\t", compression="gzip")


def _stitch_one_side(
    output_base: Path,
    combo_df: pd.DataFrame,
    *,
    filename: str,
    suffixes: Tuple[str, ...],
    side_label: str,
) -> pd.DataFrame:
    """Stitch one side (analyze or reference) by copying chr columns per combo.

    Returns a DataFrame with ``sample`` first followed by, for each chr in
    ``combo_df`` order, the ``chr{n}_<suffix>`` columns listed in ``suffixes``.
    """
    chr_list = combo_df["chr"].tolist()

    # combo -> list of chrs that pick that combo (preserve CSV order across chrs)
    combo_to_chrs: Dict[Tuple[float, float], List[str]] = {}
    for row in combo_df.itertuples(index=False):
        combo_to_chrs.setdefault((row.threshold, row.recall), []).append(row.chr)

    samples: List[str] = []
    per_chr_frames: Dict[str, pd.DataFrame] = {}

    for combo_idx, ((thres, recall), chrs) in enumerate(combo_to_chrs.items(), start=1):
        combo_dir = output_base / _format_combo_dirname(thres, recall)
        path = combo_dir / filename

        console.print(
            f"  [{side_label} {combo_idx}/{len(combo_to_chrs)}] {combo_dir.name}  "
            f"-> {len(chrs)} chr(s): {', '.join(chrs)}"
        )

        df = _load_combo_table(path).set_index("sample")

        if not samples:
            samples = df.index.astype(str).tolist()
        else:
            df = df.reindex(samples)

        for chr_name in chrs:
            cols = [f"{chr_name}_{suf}" for suf in suffixes]
            missing = [c for c in cols if c not in df.columns]
            if missing:
                raise KeyError(f"Columns missing in {path}: {missing}")
            per_chr_frames[chr_name] = df[cols].reset_index(drop=True)

    if not samples:
        raise RuntimeError("No combo files were loaded; cannot build output.")

    ordered_frames = [per_chr_frames[chr_name] for chr_name in chr_list]
    sample_col = pd.DataFrame({"sample": samples})
    out = pd.concat([sample_col] + ordered_frames, axis=1)

    out = out.sort_values("sample", kind="mergesort").reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--best-combo-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help=(
        "CSV with columns chr,threshold,recall identifying the best threshold "
        "per chr under a shared recall."
    ),
)
@click.option(
    "--output-base",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help=(
        "Directory containing the threshold_{t}_recall_{r}/ subdirectories "
        "produced by beta_to_episcore.py."
    ),
)
@click.option(
    "--output-subdir",
    default="global_recall",
    type=str,
    help="Subdirectory under --output-base where outputs are written.",
)
@click.option(
    "--no-gzip",
    is_flag=True,
    default=False,
    help="Write outputs as plain TSV instead of TSV.gz.",
)
def main(
    best_combo_csv: str,
    output_base: str,
    output_subdir: str,
    no_gzip: bool,
) -> None:
    """Stitch per-chromosome z-score columns across combos under a shared recall.

    \b
    1. Read the best-combo CSV.
    2. For each unique (threshold, recall) referenced by the CSV, open the
       analyze and reference TSVs and copy out the chr{n}_* columns for the
       chrs that pick that combo. No recomputation is performed.
    3. Write global_recall/_reference_zscore.tsv.gz and
       global_recall/_analyze_zscore.tsv.gz with the same schema as
       beta_to_episcore.py's outputs.
    """
    console.rule("[bold blue]Global-recall z-score stitching")

    output_base_path = Path(output_base)
    out_dir = output_base_path / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    console.print("\n[bold]Input parameters[/bold]")
    console.print(f"  Best-combo CSV : {best_combo_csv}")
    console.print(f"  Output base    : {output_base_path}")
    console.print(f"  Output dir     : {out_dir}")

    try:
        combo_df = _load_best_combo(best_combo_csv)
        unique_recalls = sorted(combo_df["recall"].unique().tolist())
        if len(unique_recalls) > 1:
            console.print(
                f"[yellow]Warning[/yellow] best-combo CSV uses multiple recalls: "
                f"{unique_recalls}  (expected one shared recall)"
            )
        console.print(
            f"  {len(combo_df)} chromosomes, "
            f"{combo_df[['threshold', 'recall']].drop_duplicates().shape[0]} "
            f"unique combos, recall(s)={unique_recalls}"
        )

        console.print("\n[bold cyan]Step 1: Stitching reference matrix[/bold cyan]")
        ref_df = _stitch_one_side(
            output_base_path,
            combo_df,
            filename="_reference_zscore.tsv.gz",
            suffixes=_REFERENCE_SUFFIXES,
            side_label="reference",
        )

        console.print("\n[bold cyan]Step 2: Stitching analyze matrix[/bold cyan]")
        analyze_df = _stitch_one_side(
            output_base_path,
            combo_df,
            filename="_analyze_zscore.tsv.gz",
            suffixes=_ANALYZE_SUFFIXES,
            side_label="analyze",
        )

        console.print("\n[bold cyan]Step 3: Writing outputs[/bold cyan]")
        suffix = ".tsv" if no_gzip else ".tsv.gz"
        ref_path = str(out_dir / f"_reference_zscore{suffix}")
        analyze_path = str(out_dir / f"_analyze_zscore{suffix}")

        _write_tsv(ref_df, ref_path)
        _write_tsv(analyze_df, analyze_path)

        console.print(
            f"[green]OK[/green] Reference matrix: {ref_path}  "
            f"({ref_df.shape[0]} samples x {ref_df.shape[1]} columns)"
        )
        console.print(
            f"[green]OK[/green] Analyze matrix : {analyze_path}  "
            f"({analyze_df.shape[0]} samples x {analyze_df.shape[1]} columns)"
        )

        console.rule("[bold green]Global-recall stitching complete")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
