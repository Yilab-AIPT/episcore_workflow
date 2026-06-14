#!/usr/bin/env python3
"""
Randomly draw 40 dev+Normal reference samples and compute score matrices.

For each repeat, 40 unique samples are drawn from:
    set == "dev" AND label == "Normal" AND sample != excluded sample
All remaining samples are treated as analyze. Outputs land in:
    <output-base>/randomly_select_ref_40/repeat_{index}/

Each repeat directory contains score.tsv and reference matrices, plus
reference_samples.tsv listing the draw.
"""

from __future__ import annotations

import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Sequence, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

from calc_zscore_episcore_ezscore import (
    DEFAULT_EZSCORE_REF_SAMPLES,
    build_ezscore_ref_mask,
    compute_scores,
    load_score_inputs,
    write_score_outputs,
)

console = Console()

DEFAULT_EXCLUDED = "PTAY0586P8S1"


def _generate_unique_ref_draws(
    pool_size: int,
    ref_n: int,
    n_repeats: int,
    rng: np.random.Generator,
) -> List[np.ndarray]:
    """Return ``n_repeats`` unique sorted index arrays into the candidate pool."""
    if ref_n <= 0 or ref_n > pool_size:
        raise ValueError(f"Invalid ref_n={ref_n} for pool_size={pool_size}")
    total = math.comb(int(pool_size), int(ref_n))
    if total <= n_repeats:
        raise ValueError(
            f"Pool too small for {n_repeats} unique draws: C({pool_size},{ref_n})={total}"
        )

    seen: set[bytes] = set()
    out: List[np.ndarray] = []
    max_attempts = n_repeats * 30 + 1000
    attempts = 0
    while len(out) < n_repeats and attempts < max_attempts:
        attempts += 1
        choice = rng.choice(pool_size, size=ref_n, replace=False)
        choice.sort()
        key = choice.tobytes()
        if key in seen:
            continue
        seen.add(key)
        out.append(choice.astype(np.int64, copy=False))
    if len(out) < n_repeats:
        raise RuntimeError(
            f"Could only generate {len(out)}/{n_repeats} unique reference draws "
            f"after {attempts} attempts"
        )
    return out


def _run_repeat(
    repeat_index: int,
    ref_global_idx: np.ndarray,
    merged: pd.DataFrame,
    pct_path: Path,
    chr_list: List[str],
    ezscore_ref_mask: np.ndarray,
    out_dir: Path,
) -> Tuple[int, List[str]]:
    ref_mask = np.zeros(len(merged), dtype=bool)
    ref_mask[ref_global_idx] = True
    analyze_mask = ~ref_mask
    outputs = compute_scores(merged, pct_path, chr_list, ref_mask, analyze_mask, ezscore_ref_mask)
    ref_samples = merged.loc[ref_mask, "sample"].tolist()
    write_score_outputs(out_dir, *outputs, reference_samples=ref_samples)
    return repeat_index, ref_samples


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing meta.csv, beta.csv, percentage.csv",
)
@click.option(
    "--output-base",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Base output directory (repeat dirs created under randomly_select_ref_40/)",
)
@click.option(
    "--total-repeats",
    default=1000,
    show_default=True,
    type=int,
    help="Total unique random reference draws planned for the full sweep",
)
@click.option(
    "--repeat-start",
    default=0,
    show_default=True,
    type=int,
    help="First repeat index (inclusive) for array-job slicing",
)
@click.option(
    "--repeat-end",
    default=None,
    type=int,
    help="Last repeat index (exclusive). Default: total-repeats",
)
@click.option(
    "--ref-n",
    default=40,
    show_default=True,
    type=int,
    help="Number of reference samples per repeat",
)
@click.option(
    "--seed",
    default=42,
    show_default=True,
    type=int,
    help="RNG seed for reproducible unique reference draws",
)
@click.option(
    "--exclude-sample",
    default=DEFAULT_EXCLUDED,
    show_default=True,
    help="Sample ID never selected as reference",
)
@click.option(
    "--threads",
    default=None,
    type=int,
    help="Worker threads. Defaults to SLURM_CPUS_PER_TASK or 8",
)
@click.option(
    "--ezscore-ref-samples",
    default=str(DEFAULT_EZSCORE_REF_SAMPLES),
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Sample list (one ID per line) for ezscore mean/sd reference",
)
def main(
    input_dir: str,
    output_base: str,
    total_repeats: int,
    repeat_start: int,
    repeat_end: int | None,
    ref_n: int,
    seed: int,
    exclude_sample: str,
    threads: int | None,
    ezscore_ref_samples: str,
) -> None:
    """Run random 40-sample reference score calculation repeats."""
    input_path = Path(input_dir)
    out_root = Path(output_base) / "randomly_select_ref_40"
    out_root.mkdir(parents=True, exist_ok=True)

    if repeat_end is None:
        repeat_end = total_repeats
    if repeat_start < 0 or repeat_end > total_repeats:
        raise click.ClickException(
            f"Repeat slice [{repeat_start}, {repeat_end}) must lie within [0, {total_repeats})"
        )
    if repeat_end <= repeat_start:
        raise click.ClickException("repeat-end must be greater than repeat-start")

    if threads is None or threads <= 0:
        import os

        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "0") or 0) or 8

    merged, pct_path, chr_list = load_score_inputs(input_path)
    if "label" not in merged.columns:
        raise click.ClickException("meta.csv missing column: label")

    pool_mask = (
        (merged["set"].astype(str) == "dev")
        & (merged["label"].astype(str) == "Normal")
        & (merged["sample"].astype(str) != exclude_sample)
    )
    pool_idx = np.flatnonzero(pool_mask.to_numpy())
    pool_size = int(pool_idx.size)
    if pool_size < ref_n:
        raise click.ClickException(
            f"Candidate pool size ({pool_size}) is smaller than ref-n ({ref_n})"
        )

    n_draws = repeat_end - repeat_start
    rng = np.random.default_rng(seed)
    all_draws = _generate_unique_ref_draws(pool_size, ref_n, total_repeats, rng)
    draws = all_draws[repeat_start:repeat_end]
    ezscore_ref_mask = build_ezscore_ref_mask(merged, samples_file=Path(ezscore_ref_samples))

    console.rule("[bold blue]Random ref-40 score calculation")
    console.print(f"  Input dir       : {input_path}")
    console.print(f"  Output root     : {out_root}")
    console.print(f"  Candidate pool  : {pool_size} dev+Normal (excluded: {exclude_sample})")
    console.print(
        f"  ezscore ref     : {int(ezscore_ref_mask.sum())} "
        f"(from {Path(ezscore_ref_samples).name})"
    )
    console.print(f"  Repeat range    : [{repeat_start}, {repeat_end})")
    console.print(f"  Threads         : {threads}")

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
    )
    failures: List[str] = []

    with progress:
        task = progress.add_task("Repeats", total=n_draws)
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futures = {}
            for offset, local_draw in enumerate(draws):
                repeat_index = repeat_start + offset
                ref_global_idx = pool_idx[np.asarray(local_draw, dtype=np.int64)]
                repeat_dir = out_root / f"repeat_{repeat_index}"
                fut = ex.submit(
                    _run_repeat,
                    repeat_index,
                    ref_global_idx,
                    merged,
                    pct_path,
                    chr_list,
                    ezscore_ref_mask,
                    repeat_dir,
                )
                futures[fut] = repeat_index

            for fut in as_completed(futures):
                repeat_index = futures[fut]
                try:
                    fut.result()
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"repeat_{repeat_index}: {exc}")
                progress.advance(task)

    manifest_rows = [
        {
            "repeat_index": repeat_start + offset,
            "reference_list": ",".join(merged.iloc[pool_idx[draw]]["sample"].astype(str).tolist()),
        }
        for offset, draw in enumerate(draws)
    ]
    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_root / f"reference_manifest_{repeat_start}_{repeat_end}.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    if failures:
        console.print(f"[bold red]Failures ({len(failures)}):[/bold red]")
        for msg in failures[:10]:
            console.print(f"  {msg}")
        raise click.ClickException(f"{len(failures)} repeats failed")

    console.print(f"\n[green]OK[/green] Completed {n_draws} repeats under {out_root}")
    console.print(f"[green]OK[/green] Wrote manifest {manifest_path}")
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}", style="bold red")
        sys.exit(1)
