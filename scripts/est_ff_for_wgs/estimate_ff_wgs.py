#!/usr/bin/env python3
"""
WGS fetal fraction estimation from a single SNP pileup.

Simplified variant of ``bin/estimate_ff_with_higher_precision.py`` for WGS data:
- cfDNA read counts only (no model / MQ columns)
- output column ``ff`` instead of ``ff_before_mq`` / ``ff_after_mq``
- supports ``all``, ``chr_only``, and ``chr_exclude`` analysis modes
"""

import math
import re
import sys
from multiprocessing import cpu_count
from pathlib import Path
from typing import List, Tuple

import click
import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

SCRIPT_DIR = Path(__file__).resolve().parent
BIN_DIR = SCRIPT_DIR.parent.parent / 'bin'
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from FFEstimator import FFEstimator  # noqa: E402
from estimate_ff import parse_chromosome_list, parse_mode_list  # noqa: E402


console = Console()

FF_PRECISION_PATTERN = re.compile(r'^0\.0*1$')
INITIAL_STEP = 0.01


def _validate_ff_precision(ctx, param, value: str) -> str:
    if not FF_PRECISION_PATTERN.match(value):
        raise click.BadParameter(
            f"--ff-precision must match the format '0.[0*n]1 "
            f"(e.g. 0.1, 0.01, 0.001, 0.0001). Got: {value!r}"
        )
    return value


def precision_to_decimals(precision_str: str) -> int:
    return len(precision_str) - 2


def _step_decimals(step: float) -> int:
    if step <= 0:
        return 0
    return max(0, int(round(-math.log10(step))))


def load_wgs_pileup(
    input_path: Path,
    ref_col: str = 'cfDNA_ref_reads',
    alt_col: str = 'cfDNA_alt_reads',
    min_depth: int = 0,
) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    compression = 'gzip' if str(input_path).endswith('.gz') else None
    df = pd.read_csv(input_path, sep='\t', compression=compression)

    required = {'chr', 'pos', 'af', ref_col, alt_col}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required - set(df.columns)}")
    if len(df) == 0:
        raise ValueError("Input file is empty")

    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    df['af'] = pd.to_numeric(df['af'], errors='coerce')
    df[ref_col] = pd.to_numeric(df[ref_col], errors='coerce', downcast='integer')
    df[alt_col] = pd.to_numeric(df[alt_col], errors='coerce', downcast='integer')

    invalid_mask = (
        df['af'].isna() | (df['af'] < 0) | (df['af'] > 1)
        | df[ref_col].isna() | (df[ref_col] < 0)
        | df[alt_col].isna() | (df[alt_col] < 0)
    )
    if invalid_mask.any():
        console.print(f"[yellow]Warning: Removing {invalid_mask.sum()} rows with invalid data[/yellow]")
        df = df[~invalid_mask]

    if min_depth > 0:
        initial_count = len(df)
        df = df[(df[ref_col] + df[alt_col]) >= min_depth].reset_index(drop=True)
        remaining = len(df)
        pct = (remaining / initial_count * 100) if initial_count else 0.0
        console.print(
            f"[cyan]Depth filter (≥{min_depth}): {remaining:,} SNPs remaining ({pct:.1f}%)[/cyan]"
        )

    if len(df) == 0:
        raise ValueError("No valid SNPs remaining after filtering")

    return df


def estimate_ff_iterative(
    estimator: FFEstimator,
    df: pd.DataFrame,
    *,
    ff_min: float,
    ff_max: float,
    ff_precision: float,
    ncpus: int,
    ref_col: str,
    alt_col: str,
) -> Tuple[float, List[dict]]:
    current_step = INITIAL_STEP
    current_min = ff_min
    current_max = ff_max
    history: List[dict] = []
    iteration = 0

    while True:
        iteration += 1
        if current_max <= current_min:
            current_max = min(1.0, current_min + current_step)

        step_decimals = _step_decimals(current_step)
        console.print(
            f"  [cyan]Round {iteration}:[/cyan] "
            f"range [{current_min:.{step_decimals}f}, {current_max:.{step_decimals}f}], "
            f"step {current_step:.{step_decimals}f}"
        )

        best_ff, _ = estimator.estimate(
            df,
            f_min=float(current_min),
            f_max=float(current_max),
            f_step=float(current_step),
            ncpus=ncpus,
            ref_col=ref_col,
            alt_col=alt_col,
        )
        best_ff = round(float(best_ff), step_decimals)
        history.append({
            'iteration': iteration,
            'step': current_step,
            'range': (current_min, current_max),
            'best': best_ff,
        })

        if current_step <= ff_precision + 1e-15:
            break

        next_step = current_step / 10
        current_min = max(ff_min, best_ff - current_step)
        current_max = min(ff_max, best_ff + current_step)
        current_step = next_step

    return best_ff, history


def display_results_summary(
    results_df: pd.DataFrame,
    modes: List[str],
    n_decimals: int,
) -> None:
    mode_str = ', '.join(modes) if len(modes) > 1 else modes[0]
    table = Table(title=f"Fetal Fraction Results - Modes: {mode_str}")
    table.add_column("Chromosome", justify="center", style="cyan")
    table.add_column("FF", justify="right", style="green")

    for _, row in results_df.iterrows():
        table.add_row(str(row['chr']), f"{row['ff']:.{n_decimals}f}")

    if len(results_df) > 1:
        for mode in modes:
            if mode == 'all':
                continue
            mode_suffix = f"_{mode}"
            mode_rows = results_df[results_df['chr'].str.endswith(mode_suffix)]
            if len(mode_rows) > 0:
                mean_ff = mode_rows['ff'].mean()
                table.add_row(
                    f"[bold]Mean ({mode})[/bold]",
                    f"[bold]{mean_ff:.{n_decimals}f}[/bold]",
                )

    console.print(table)


@click.command()
@click.option('--input-path', '-i', required=True, type=click.Path(exists=True, path_type=Path))
@click.option('--output-prefix', '-o', required=True, type=str)
@click.option('--ff-min', type=click.FloatRange(0.0, 1.0), default=0.0, show_default=True)
@click.option('--ff-max', type=click.FloatRange(0.0, 1.0), default=0.3, show_default=True)
@click.option(
    '--ff-precision',
    type=str,
    default='0.0001',
    callback=_validate_ff_precision,
    show_default=True,
)
@click.option('--ref-col', default='cfDNA_ref_reads', show_default=True)
@click.option('--alt-col', default='cfDNA_alt_reads', show_default=True)
@click.option('--min-depth', type=click.IntRange(0, None), default=3, show_default=True)
@click.option(
    '--chromosomes',
    default='1-22',
    show_default=True,
    help='Chromosomes for chr_only / chr_exclude (e.g. "1-22", "1,2,3", "21")',
)
@click.option(
    '--mode-list',
    type=str,
    default='all',
    show_default=True,
    help='Comma-separated modes: chr_only, chr_exclude, all',
)
@click.option('--ncpus', type=click.IntRange(1, cpu_count()), default=cpu_count(), show_default=True)
@click.option('--verbose', '-v', is_flag=True)
def main(
    input_path: Path,
    output_prefix: str,
    ff_min: float,
    ff_max: float,
    ff_precision: str,
    ref_col: str,
    alt_col: str,
    min_depth: int,
    chromosomes: str,
    mode_list: str,
    ncpus: int,
    verbose: bool,
) -> None:
    ff_precision_val = float(ff_precision)
    n_decimals = precision_to_decimals(ff_precision)
    float_fmt = f"%.{n_decimals}f"

    try:
        target_chromosomes = parse_chromosome_list(chromosomes)
        analysis_modes = parse_mode_list(mode_list)

        console.print(Panel.fit(
            f"[bold green]WGS Fetal Fraction Estimator[/bold green]\n"
            f"Input: {input_path}\n"
            f"Output: {output_prefix}_ff.tsv\n"
            f"Analysis modes: {', '.join(analysis_modes)}\n"
            f"Chromosomes: {', '.join(map(str, target_chromosomes))}\n"
            f"FF range: {ff_min:.{n_decimals}f} - {ff_max:.{n_decimals}f}\n"
            f"FF precision: {ff_precision}\n"
            f"Min depth: {min_depth}\n"
            f"CPU cores: {ncpus}",
            title="Configuration",
        ))

        console.print("[cyan]Loading pileup data...[/cyan]")
        df = load_wgs_pileup(
            input_path,
            ref_col=ref_col,
            alt_col=alt_col,
            min_depth=min_depth,
        )
        console.print(
            f"[green]✓ Loaded {len(df):,} SNPs from {df['chr'].nunique()} chromosomes[/green]"
        )

        estimator = FFEstimator(mode='cfDNA')
        all_results_list: List[dict] = []

        for analysis_mode in analysis_modes:
            console.print(f"\n[bold magenta]{'=' * 60}[/bold magenta]")
            console.print(f"[bold magenta]Running analysis mode: {analysis_mode}[/bold magenta]")
            console.print(f"[bold magenta]{'=' * 60}[/bold magenta]")

            if analysis_mode == 'all':
                console.print(f"\n[bold cyan]Processing all chromosomes together[/bold cyan]")
                console.print(f"Total SNPs: {len(df):,}")
                try:
                    ff_value, _ = estimate_ff_iterative(
                        estimator,
                        df,
                        ff_min=ff_min,
                        ff_max=ff_max,
                        ff_precision=ff_precision_val,
                        ncpus=ncpus,
                        ref_col=ref_col,
                        alt_col=alt_col,
                    )
                    all_results_list.append({'chr': 'all', 'ff': ff_value})
                    console.print(f"[green]✓ all: FF = {ff_value:.{n_decimals}f}[/green]")
                except Exception as e:
                    console.print(f"[red]✗ Error processing all mode: {e}[/red]")
                    if verbose:
                        console.print_exception()
            else:
                for target_chr in target_chromosomes:
                    chr_name = f"chr{target_chr}"
                    console.print(f"\n[bold cyan]Processing chromosome {target_chr}[/bold cyan]")

                    target_data = df[df['chr'] == chr_name]
                    background_data = df[df['chr'] != chr_name]

                    if len(target_data) == 0:
                        console.print(f"[yellow]Warning: No SNPs found for {chr_name}, skipping[/yellow]")
                        continue

                    if analysis_mode == 'chr_only':
                        estimation_data = target_data
                        chr_label = f"{chr_name}_only"
                        data_description = f"Target SNPs only: {len(estimation_data):,}"
                    else:
                        if len(background_data) == 0:
                            console.print(
                                f"[red]Error: No background SNPs available for {chr_name}[/red]"
                            )
                            continue
                        estimation_data = background_data
                        chr_label = f"{chr_name}_exclude"
                        data_description = (
                            f"Target SNPs: {len(target_data):,}, "
                            f"Background SNPs: {len(background_data):,}"
                        )

                    console.print(data_description)
                    try:
                        ff_value, _ = estimate_ff_iterative(
                            estimator,
                            estimation_data,
                            ff_min=ff_min,
                            ff_max=ff_max,
                            ff_precision=ff_precision_val,
                            ncpus=ncpus,
                            ref_col=ref_col,
                            alt_col=alt_col,
                        )
                        all_results_list.append({'chr': chr_label, 'ff': ff_value})
                        console.print(f"[green]✓ {chr_label}: FF = {ff_value:.{n_decimals}f}[/green]")
                    except Exception as e:
                        console.print(f"[red]✗ Error processing {chr_name}: {e}[/red]")
                        if verbose:
                            console.print_exception()

        if not all_results_list:
            console.print("[red]✗ No results generated - check input data and parameters[/red]")
            sys.exit(1)

        output_path = Path(f'{output_prefix}_ff.tsv')
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results_df = pd.DataFrame(all_results_list, columns=['chr', 'ff'])
        results_df['ff'] = results_df['ff'].round(n_decimals)
        results_df.to_csv(output_path, sep='\t', index=False, float_format=float_fmt)

        display_results_summary(results_df, modes=analysis_modes, n_decimals=n_decimals)
        console.print(f"[bold green]✓ Results saved to {output_path}[/bold green]")

    except Exception as e:
        console.print(f"[red]✗ Fatal error: {e}[/red]")
        if verbose:
            console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    import multiprocessing as mp

    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass

    main()
