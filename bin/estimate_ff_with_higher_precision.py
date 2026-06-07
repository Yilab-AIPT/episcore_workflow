#!/usr/bin/env python3
"""
High-Precision Fetal Fraction Estimation for a single SNP pileup.

This script is a higher-precision variant of ``bin/estimate_ff.py`` adapted to
the Nextflow ``snp_est_ff`` workflow. It uses the same likelihood model
implemented in ``bin/FFEstimator.py`` (Hardy-Weinberg mixture model over
maternal/fetal genotypes), but replaces the single fixed-grid search with an
iterative range-narrowing search:

  1. Search [ff_min, ff_max] with step = 0.01 -> best candidate b1
  2. Search [b1 - 0.01, b1 + 0.01] with step = 0.001 -> b2
  3. Search [b2 - 0.001, b2 + 0.001] with step = 0.0001 -> b3
  4. ... continue until the step reaches ``--ff-precision``.

Differences vs. the offline ``scripts/ff_decimal`` variant:

- Operates on a SINGLE input pileup (``--input-path``) instead of scanning a
  directory of pileups, so it slots directly into a per-sample Nextflow process.
- Adds ``--known-sites`` to optionally restrict the pileup to panel SNP sites
  before estimation. This shares the same file format / parameter as
  ``bam_to_pileup.py`` (``params.snp_list``).

The CLI matches ``bin/estimate_ff.py`` except that ``--ff-step`` is replaced by
``--ff-precision`` (a decimal precision specifier of the form ``0.[0*n]1``,
e.g. ``0.1``, ``0.01``, ``0.001``, ``0.0001``) and ``--known-sites`` is added.
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

# ``FFEstimator.py`` and ``estimate_ff.py`` live alongside this script in
# ``bin/``. Nextflow places the ``bin/`` directory of the pipeline on PATH and
# Python adds the script's own directory to ``sys.path``, so these direct
# imports resolve both inside the container and when run locally.
from FFEstimator import FFEstimator  # noqa: E402
from estimate_ff import (  # noqa: E402
    load_and_validate_data,
    parse_chromosome_list,
    parse_mode_list,
)


console = Console()


FF_PRECISION_PATTERN = re.compile(r'^0\.0*1$')
INITIAL_STEP = 0.01


def _validate_ff_precision(ctx, param, value: str) -> str:
    """Click callback that enforces the ``0.[0*n]1`` precision format."""
    if not FF_PRECISION_PATTERN.match(value):
        raise click.BadParameter(
            f"--ff-precision must match the format '0.[0*n]1' "
            f"(e.g. 0.1, 0.01, 0.001, 0.0001). Got: {value!r}"
        )
    return value


def precision_to_decimals(precision_str: str) -> int:
    """Return the number of decimal places implied by a precision string.

    Examples:
        "0.1"    -> 1
        "0.01"   -> 2
        "0.0001" -> 4
    """
    return len(precision_str) - 2


def _step_decimals(step: float) -> int:
    """Number of decimal places for a step that is a negative power of 10."""
    if step <= 0:
        return 0
    return max(0, int(round(-math.log10(step))))


def parse_known_sites(sites_file: Path) -> pd.DataFrame:
    """Parse a known SNP sites file (VCF-like TSV) into a DataFrame.

    The format matches the ``--known-sites`` file consumed by
    ``bam_to_pileup.py``: a tab-separated, optionally ``#``-commented table
    whose columns 0/1/3/4 are ``chr``, ``pos``, ``ref``, ``alt``. Only
    single-nucleotide variants are retained.

    Args:
        sites_file: Path to the known sites TSV file.

    Returns:
        DataFrame with columns ``chr`` (str), ``pos`` (int), ``ref``, ``alt``.
    """
    sites = pd.read_csv(
        sites_file,
        sep='\t',
        comment='#',
        usecols=[0, 1, 3, 4],
        names=['chr', 'pos', 'ref', 'alt'],
        dtype={'chr': str, 'pos': int, 'ref': str, 'alt': str},
    )
    sites['ref'] = sites['ref'].str.upper()
    sites['alt'] = sites['alt'].str.upper()

    # Keep single-nucleotide variants only (consistent with bam_to_pileup.py).
    single_nuc_mask = (sites['ref'].str.len() == 1) & (sites['alt'].str.len() == 1)
    sites = sites[single_nuc_mask]
    return sites


def filter_by_known_sites(df: pd.DataFrame, sites_file: Path) -> pd.DataFrame:
    """Restrict a pileup DataFrame to positions present in ``sites_file``.

    Filtering is done on the ``(chr, pos)`` key, matching how
    ``bam_to_pileup.py`` writes the pileup ``chr``/``pos`` columns directly from
    the known sites file (so chromosome naming is consistent between the two).

    Args:
        df: Pileup DataFrame (must contain ``chr`` and ``pos`` columns).
        sites_file: Path to the known sites TSV file.

    Returns:
        Filtered DataFrame containing only rows whose ``(chr, pos)`` appear in
        the known sites panel.
    """
    sites = parse_known_sites(sites_file)

    df = df.copy()
    df['chr'] = df['chr'].astype(str)
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce').astype('Int64')

    key_df = sites[['chr', 'pos']].drop_duplicates()
    key_df['chr'] = key_df['chr'].astype(str)
    key_df['pos'] = key_df['pos'].astype('Int64')

    n_before = len(df)
    filtered = df.merge(key_df, on=['chr', 'pos'], how='inner').reset_index(drop=True)
    n_after = len(filtered)
    pct = (n_after / n_before * 100) if n_before > 0 else 0.0
    console.print(
        f"[cyan]Known-sites filter ({len(key_df):,} panel sites): "
        f"{n_after:,}/{n_before:,} pileup SNPs retained ({pct:.1f}%)[/cyan]"
    )

    if n_after == 0:
        raise ValueError("No pileup SNPs remaining after known-sites filtering")

    return filtered


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
    label: str = "",
) -> Tuple[float, List[dict]]:
    """Iteratively narrow the FF search range until ``ff_precision`` is met.

    Returns ``(best_ff, history)``. ``history`` is a list of per-round dicts
    with keys ``step``, ``range``, and ``best`` (useful for diagnostics).
    """
    if not (0 <= ff_min < ff_max <= 1):
        raise ValueError(
            f"Invalid FF search range: ff_min={ff_min}, ff_max={ff_max}"
        )
    if ff_precision <= 0:
        raise ValueError(f"ff_precision must be positive: {ff_precision}")

    # Always start at the coarse 0.01 grid as specified in the task.
    current_step = INITIAL_STEP
    current_min = ff_min
    current_max = ff_max

    history: List[dict] = []
    iteration = 0
    while True:
        iteration += 1

        # Guard against degenerate ranges from boundary candidates.
        if current_max <= current_min:
            current_max = min(1.0, current_min + current_step)

        step_decimals = _step_decimals(current_step)
        console.print(
            f"  [cyan]Round {iteration}{(' ' + label) if label else ''}:[/cyan] "
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

        # Round to the current step's decimals to avoid floating-point drift
        # propagating into the next iteration's search range.
        best_ff = round(float(best_ff), step_decimals)
        history.append({
            'iteration': iteration,
            'step': current_step,
            'range': (current_min, current_max),
            'best': best_ff,
        })

        # Stop once the search step is at or finer than the requested precision.
        if current_step <= ff_precision + 1e-15:
            break

        # Narrow the window around the best candidate and refine the step.
        next_step = current_step / 10
        current_min = max(ff_min, best_ff - current_step)
        current_max = min(ff_max, best_ff + current_step)
        current_step = next_step

    return best_ff, history


@click.command()
@click.option(
    '--input-path', '-i',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to input pileup TSV.GZ file containing SNP data',
)
@click.option(
    '--output-prefix', '-o',
    type=str,
    required=True,
    help='Output prefix for results (will create {prefix}_ff.tsv)',
)
@click.option(
    '--known-sites',
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help='Optional known SNP sites TSV (same format as bam_to_pileup.py). '
         'When given, the pileup is filtered to these panel sites before estimation.',
)
@click.option(
    '--ff-min',
    type=click.FloatRange(0.0, 1.0),
    default=0.0,
    help='Minimum fetal fraction for estimation (default: 0.0)',
)
@click.option(
    '--ff-max',
    type=click.FloatRange(0.0, 1.0),
    default=0.3,
    help='Maximum fetal fraction for estimation (default: 0.3)',
)
@click.option(
    '--ff-precision',
    type=str,
    default='0.0001',
    callback=_validate_ff_precision,
    help='Decimal precision for FF (format "0.[0*n]1"; default: 0.0001)',
)
@click.option(
    '--chromosomes',
    default='1-22',
    help='Chromosomes to analyze (e.g., "1-22", "1,2,3", or "21"). Default: 1-22',
)
@click.option(
    '--cfdna-ref-col',
    default='cfDNA_ref_reads',
    help='Column name for cfDNA reference reads (default: cfDNA_ref_reads)',
)
@click.option(
    '--cfdna-alt-col',
    default='cfDNA_alt_reads',
    help='Column name for cfDNA alternative reads (default: cfDNA_alt_reads)',
)
@click.option(
    '--model-ref-col',
    default='fetal_ref_reads_from_model',
    help='Column name for modeled fetal reference reads (default: fetal_ref_reads_from_model)',
)
@click.option(
    '--model-alt-col',
    default='fetal_alt_reads_from_model',
    help='Column name for modeled fetal alternative reads (default: fetal_alt_reads_from_model)',
)
@click.option(
    '--min-raw-depth',
    type=click.IntRange(0, None),
    default=0,
    help='Minimum raw depth filter for cfDNA reads (default: 0)',
)
@click.option(
    '--min-model-depth',
    type=click.IntRange(0, None),
    default=0,
    help='Minimum model depth filter for model filtered reads (default: 0)',
)
@click.option(
    '--mode-list',
    type=str,
    default='chr_exclude',
    help='Comma-separated list of analysis modes: chr_only, chr_exclude, all (default: chr_exclude)',
)
@click.option(
    '--ncpus',
    type=click.IntRange(1, cpu_count()),
    default=cpu_count(),
    help=f'Number of CPU cores to use for parallel processing (default: {cpu_count()})',
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output',
)
def main(
    input_path: Path,
    output_prefix: str,
    known_sites: Path,
    ff_min: float,
    ff_max: float,
    ff_precision: str,
    chromosomes: str,
    cfdna_ref_col: str,
    cfdna_alt_col: str,
    model_ref_col: str,
    model_alt_col: str,
    min_raw_depth: int,
    min_model_depth: int,
    mode_list: str,
    ncpus: int,
    verbose: bool,
) -> None:
    """High-precision Fetal Fraction Estimation for a single pileup.

    Follows the same logic as ``bin/estimate_ff.py`` (cfDNA mode + cfDNA+model
    mode, supports chr_only / chr_exclude / all analysis modes) but performs an
    iterative range-narrowing grid search controlled by ``--ff-precision``
    instead of a fixed-step grid. Optionally restricts the pileup to a known
    SNP panel via ``--known-sites`` before estimating FF.
    """
    if verbose:
        console.print("[blue]Verbose mode enabled[/blue]")

    ff_precision_str = ff_precision
    ff_precision_val = float(ff_precision_str)
    n_decimals = precision_to_decimals(ff_precision_str)
    float_fmt = f"%.{n_decimals}f"

    try:
        target_chromosomes = parse_chromosome_list(chromosomes)
        analysis_modes = parse_mode_list(mode_list)

        column_info = (
            f"cfDNA columns: {cfdna_ref_col}, {cfdna_alt_col}\n"
            f"Model columns: {model_ref_col}, {model_alt_col}"
        )
        depth_filter_info = (
            f"Raw depth filter: ≥{min_raw_depth}\n"
            f"Model depth filter: ≥{min_model_depth}"
        )

        console.print(Panel.fit(
            f"[bold green]High-Precision Fetal Fraction Estimator[/bold green]\n"
            f"Input: {input_path}\n"
            f"Output: {output_prefix}_ff.tsv\n"
            f"Known sites: {known_sites if known_sites else 'none (no filtering)'}\n"
            f"Analysis modes: {', '.join(analysis_modes)}\n"
            f"Modes: cfDNA (before MQ) + cfDNA+model (after MQ)\n"
            f"{column_info}\n"
            f"Chromosomes: {', '.join(map(str, target_chromosomes))}\n"
            f"FF Range: {ff_min:.{n_decimals}f} - {ff_max:.{n_decimals}f}\n"
            f"FF Precision: {ff_precision_str} ({n_decimals} decimals)\n"
            f"{depth_filter_info}\n"
            f"CPU Cores: {ncpus}",
            title="Configuration",
        ))

        console.print("[cyan]Loading input data...[/cyan]")
        df = load_and_validate_data(
            input_path,
            cfdna_ref_col=cfdna_ref_col,
            cfdna_alt_col=cfdna_alt_col,
            model_ref_col=model_ref_col,
            model_alt_col=model_alt_col,
            min_raw_depth=min_raw_depth,
            min_model_depth=min_model_depth,
        )

        console.print(
            f"[green]✓ Loaded {len(df)} SNPs from "
            f"{df['chr'].nunique()} chromosomes[/green]"
        )

        # Optionally restrict the pileup to the known SNP panel before estimation.
        if known_sites is not None:
            df = filter_by_known_sites(df, known_sites)
            console.print(
                f"[green]✓ {len(df)} SNPs from "
                f"{df['chr'].nunique()} chromosomes after known-sites filter[/green]"
            )

        ff_estimator_cfdna = FFEstimator(mode='cfDNA')
        ff_estimator_model = FFEstimator(mode='cfDNA+model')

        all_results_list: List[dict] = []

        for analysis_mode in analysis_modes:
            console.print(f"\n[bold magenta]{'=' * 60}[/bold magenta]")
            console.print(
                f"[bold magenta]Running analysis mode: {analysis_mode}[/bold magenta]"
            )
            console.print(f"[bold magenta]{'=' * 60}[/bold magenta]")

            if analysis_mode == 'all':
                console.print("\n[bold cyan]Processing all chromosomes together[/bold cyan]")
                console.print(f"Total SNPs: {len(df)}")

                try:
                    console.print("[cyan]Estimating FF (before MQ) using all data...[/cyan]")
                    ff_before_mq, _ = estimate_ff_iterative(
                        ff_estimator_cfdna,
                        df,
                        ff_min=ff_min,
                        ff_max=ff_max,
                        ff_precision=ff_precision_val,
                        ncpus=ncpus,
                        ref_col=cfdna_ref_col,
                        alt_col=cfdna_alt_col,
                        label="(before MQ, all)",
                    )

                    console.print("[cyan]Estimating FF (after MQ) using all data...[/cyan]")
                    ff_after_mq, _ = estimate_ff_iterative(
                        ff_estimator_model,
                        df,
                        ff_min=ff_min,
                        ff_max=ff_max,
                        ff_precision=ff_precision_val,
                        ncpus=ncpus,
                        ref_col=model_ref_col,
                        alt_col=model_alt_col,
                        label="(after MQ, all)",
                    )

                    all_results_list.append({
                        'chr': 'all',
                        'ff_before_mq': ff_before_mq,
                        'ff_after_mq': ff_after_mq,
                    })
                    console.print(
                        f"[green]✓ All: FF_before = {ff_before_mq:.{n_decimals}f}, "
                        f"FF_after = {ff_after_mq:.{n_decimals}f}[/green]"
                    )

                except Exception as e:
                    console.print(f"[red]✗ Error processing all mode: {str(e)}[/red]")
                    if verbose:
                        console.print_exception()

            else:
                for target_chr in target_chromosomes:
                    chr_name = f"chr{target_chr}"
                    console.print(
                        f"\n[bold cyan]Processing chromosome {target_chr}[/bold cyan]"
                    )

                    target_data = df[df['chr'] == chr_name]
                    background_data = df[df['chr'] != chr_name]

                    if len(target_data) == 0:
                        console.print(
                            f"[yellow]Warning: No SNPs found for {chr_name}, skipping[/yellow]"
                        )
                        continue

                    if analysis_mode == 'chr_only':
                        estimation_data = target_data
                        data_description = (
                            f"Target SNPs only: {len(estimation_data)}"
                        )
                        chr_label = f"{chr_name}_only"
                    elif analysis_mode == 'chr_exclude':
                        if len(background_data) == 0:
                            console.print(
                                f"[red]Error: No background SNPs available for {chr_name}[/red]"
                            )
                            continue
                        estimation_data = background_data
                        data_description = (
                            f"Target SNPs: {len(target_data)}, "
                            f"Background SNPs: {len(background_data)}"
                        )
                        chr_label = f"{chr_name}_exclude"

                    console.print(data_description)

                    try:
                        console.print(
                            f"[cyan]Estimating FF (before MQ) for {chr_name}...[/cyan]"
                        )
                        ff_before_mq, _ = estimate_ff_iterative(
                            ff_estimator_cfdna,
                            estimation_data,
                            ff_min=ff_min,
                            ff_max=ff_max,
                            ff_precision=ff_precision_val,
                            ncpus=ncpus,
                            ref_col=cfdna_ref_col,
                            alt_col=cfdna_alt_col,
                            label=f"(before MQ, {chr_label})",
                        )

                        console.print(
                            f"[cyan]Estimating FF (after MQ) for {chr_name}...[/cyan]"
                        )
                        ff_after_mq, _ = estimate_ff_iterative(
                            ff_estimator_model,
                            estimation_data,
                            ff_min=ff_min,
                            ff_max=ff_max,
                            ff_precision=ff_precision_val,
                            ncpus=ncpus,
                            ref_col=model_ref_col,
                            alt_col=model_alt_col,
                            label=f"(after MQ, {chr_label})",
                        )

                        all_results_list.append({
                            'chr': chr_label,
                            'ff_before_mq': ff_before_mq,
                            'ff_after_mq': ff_after_mq,
                        })

                        console.print(
                            f"[green]✓ {chr_label}: "
                            f"FF_before = {ff_before_mq:.{n_decimals}f}, "
                            f"FF_after = {ff_after_mq:.{n_decimals}f}[/green]"
                        )

                    except Exception as e:
                        console.print(
                            f"[red]✗ Error processing {chr_name}: {str(e)}[/red]"
                        )
                        if verbose:
                            console.print_exception()
                        continue

        if all_results_list:
            output_path = Path(f'{output_prefix}_ff.tsv')
            console.print(f"\n[cyan]Saving all results to {output_path}...[/cyan]")
            results_df = pd.DataFrame(
                all_results_list,
                columns=['chr', 'ff_before_mq', 'ff_after_mq'],
            )
            # Round to ff_precision decimals so the stored values exactly match
            # the requested precision (avoids surprises from float printing).
            results_df['ff_before_mq'] = results_df['ff_before_mq'].round(n_decimals)
            results_df['ff_after_mq'] = results_df['ff_after_mq'].round(n_decimals)

            output_dir = output_path.parent
            if not output_dir.exists():
                output_dir.mkdir(parents=True, exist_ok=True)
            results_df.to_csv(
                output_path,
                sep='\t',
                index=False,
                float_format=float_fmt,
            )

            display_results_summary(results_df, modes=analysis_modes, n_decimals=n_decimals)
            console.print(
                f"[bold green]✓ Analysis complete! Results saved to {output_path}[/bold green]"
            )
        else:
            console.print(
                "[red]✗ No results generated - check input data and parameters[/red]"
            )
            sys.exit(1)

    except Exception as e:
        console.print(f"[red]✗ Fatal error: {str(e)}[/red]")
        if verbose:
            console.print_exception()
        sys.exit(1)


def display_results_summary(
    results_df: pd.DataFrame,
    modes: list = None,
    n_decimals: int = 4,
) -> None:
    """Display a formatted summary table of fetal fraction estimation results."""
    if modes is None:
        modes = ['chr_exclude']

    mode_str = ', '.join(modes) if len(modes) > 1 else modes[0]
    table = Table(title=f"Fetal Fraction Estimation Results - Modes: {mode_str}")
    table.add_column("Chromosome", justify="center", style="cyan")
    table.add_column("FF Before MQ", justify="right", style="green")
    table.add_column("FF After MQ", justify="right", style="yellow")

    for _, row in results_df.iterrows():
        table.add_row(
            str(row['chr']),
            f"{row['ff_before_mq']:.{n_decimals}f}",
            f"{row['ff_after_mq']:.{n_decimals}f}",
        )

    if len(results_df) > 1:
        for mode in modes:
            if mode == 'all':
                continue
            mode_suffix = f"_{mode}"
            mode_rows = results_df[results_df['chr'].str.endswith(mode_suffix)]
            if len(mode_rows) > 0:
                mean_before = mode_rows['ff_before_mq'].mean()
                mean_after = mode_rows['ff_after_mq'].mean()
                table.add_row(
                    f"[bold]Mean ({mode})[/bold]",
                    f"[bold]{mean_before:.{n_decimals}f}[/bold]",
                    f"[bold]{mean_after:.{n_decimals}f}[/bold]",
                )

    console.print(table)


if __name__ == "__main__":
    import multiprocessing as mp

    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass

    main()
