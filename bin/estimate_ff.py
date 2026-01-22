#!/usr/bin/env python3
"""
Fetal Fraction Estimation for Aneuploidy Detection

This script analyzes cell-free DNA (cfDNA) sequencing data to estimate fetal
fraction from maternal plasma samples using two different approaches:
1. Standard cfDNA analysis (before maternal quality filtering)
2. Model-based analysis using quality-filtered reads (after maternal quality filtering)

The analysis uses SNP read counts and population allele frequencies to perform
maximum likelihood estimation of fetal fraction across chromosomes.
"""

import numpy as np
import pandas as pd
from multiprocessing import cpu_count
import click
import sys
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

# Import the FFEstimator class
from FFEstimator import FFEstimator


# Initialize rich console for beautiful output
console = Console()


@click.command()
@click.option(
    '--input-path', '-i',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to input TSV.GZ file containing SNP data'
)
@click.option(
    '--output-prefix', '-o',
    type=str,
    required=True,
    help='Output prefix for results (will create {prefix}_ff.tsv)'
)
@click.option(
    '--ff-min',
    type=click.FloatRange(0.0, 1.0),
    default=0.001,
    help='Minimum fetal fraction for estimation (default: 0.001)'
)
@click.option(
    '--ff-max',
    type=click.FloatRange(0.0, 1.0),
    default=0.3,
    help='Maximum fetal fraction for estimation (default: 0.3)'
)
@click.option(
    '--ff-step',
    type=click.FloatRange(0.0001, 0.1),
    default=0.001,
    help='Step size for fetal fraction grid search (default: 0.001)'
)
@click.option(
    '--chromosomes',
    default='1-22',
    help='Chromosomes to analyze (e.g., "1-22", "1,2,3", or "21"). Default: 1-22'
)
@click.option(
    '--cfdna-ref-col',
    default='cfDNA_ref_reads',
    help='Column name for cfDNA reference reads (default: cfDNA_ref_reads)'
)
@click.option(
    '--cfdna-alt-col',
    default='cfDNA_alt_reads',
    help='Column name for cfDNA alternative reads (default: cfDNA_alt_reads)'
)
@click.option(
    '--model-ref-col',
    default='fetal_ref_reads_from_model',
    help='Column name for modeled fetal reference reads (default: fetal_ref_reads_from_model)'
)
@click.option(
    '--model-alt-col',
    default='fetal_alt_reads_from_model',
    help='Column name for modeled fetal alternative reads (default: fetal_alt_reads_from_model)'
)
@click.option(
    '--min-raw-depth',
    type=click.IntRange(0, None),
    default=0,
    help='Minimum raw depth filter for cfDNA reads (default: 0)'
)
@click.option(
    '--min-model-depth',
    type=click.IntRange(0, None),
    default=0,
    help='Minimum model depth filter for model filtered reads (default: 0)'
)
@click.option(
    '--mode-list',
    type=str,
    default='chr_exclude',
    help='Comma-separated list of analysis modes: chr_only, chr_exclude, all (default: chr_exclude)'
)
@click.option(
    '--ncpus',
    type=click.IntRange(1, cpu_count()),
    default=cpu_count(),
    help=f'Number of CPU cores to use for parallel processing (default: {cpu_count()})'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output'
)
def main(
    input_path: Path,
    output_prefix: str,
    ff_min: float,
    ff_max: float,
    ff_step: float,
    chromosomes: str,
    cfdna_ref_col: str,
    cfdna_alt_col: str,
    model_ref_col: str,
    model_alt_col: str,
    min_raw_depth: int,
    min_model_depth: int,
    mode_list: str,
    ncpus: int,
    verbose: bool
) -> None:
    """
    Fetal Fraction Estimation Tool.
    
    This tool analyzes cell-free DNA sequencing data to estimate fetal fraction
    using two different approaches:
    
    1. cfDNA mode (ff_before_mq): Standard analysis using raw cfDNA reads
    2. cfDNA+model mode (ff_after_mq): Analysis using quality-filtered model reads
    
    The input file should be a TSV.GZ file with columns:
    chr, pos, af, and read count columns (names configurable via CLI options)
    
    Analysis Modes (--mode-list):
    - chr_only: Use only the target chromosome's data to estimate FF
    - chr_exclude: Use background chromosomes (exclude target) to estimate FF
    - all: Use all chromosomes to estimate a single global FF
    
    Multiple modes can be specified (comma-separated) and will run sequentially.
    
    Results are saved as a TSV file with columns: chr, ff_before_mq, ff_after_mq
    
    Multi-threading is used by default to speed up fetal fraction estimation.
    Use --ncpus to control the number of CPU cores used.
    """
    # Configure console output
    if verbose:
        console.print("[blue]Verbose mode enabled[/blue]")
    
    try:
        # Validate and parse chromosome specification
        target_chromosomes = parse_chromosome_list(chromosomes)
        
        # Parse and validate mode list
        analysis_modes = parse_mode_list(mode_list)
        
        # Build column info for display
        column_info = (
            f"cfDNA columns: {cfdna_ref_col}, {cfdna_alt_col}\n"
            f"Model columns: {model_ref_col}, {model_alt_col}"
        )
        
        # Build depth filter info
        depth_filter_info = f"Raw depth filter: ≥{min_raw_depth}\nModel depth filter: ≥{min_model_depth}"
        
        # Display startup information
        console.print(Panel.fit(
            f"[bold green]Fetal Fraction Estimator[/bold green]\n"
            f"Input: {input_path}\n"
            f"Output: {output_prefix}_ff.tsv\n"
            f"Analysis modes: {', '.join(analysis_modes)}\n"
            f"Modes: cfDNA (before MQ) + cfDNA+model (after MQ)\n"
            f"{column_info}\n"
            f"Chromosomes: {', '.join(map(str, target_chromosomes))}\n"
            f"FF Range: {ff_min:.3f} - {ff_max:.3f} (step: {ff_step:.3f})\n"
            f"{depth_filter_info}\n"
            f"CPU Cores: {ncpus}",
            title="Configuration"
        ))
        
        # Load and validate input data
        console.print("[cyan]Loading input data...[/cyan]")
        df = load_and_validate_data(
            input_path,
            cfdna_ref_col=cfdna_ref_col,
            cfdna_alt_col=cfdna_alt_col,
            model_ref_col=model_ref_col,
            model_alt_col=model_alt_col,
            min_raw_depth=min_raw_depth,
            min_model_depth=min_model_depth
        )
        
        console.print(f"[green]✓ Loaded {len(df)} SNPs from {df['chr'].nunique()} chromosomes[/green]")
        
        # Initialize FFEstimators for both modes
        ff_estimator_cfdna = FFEstimator(mode='cfDNA')
        ff_estimator_model = FFEstimator(mode='cfDNA+model')
        
        # Initialize combined results storage for all modes
        all_results_list = []
        
        # Run analysis for each mode
        for analysis_mode in analysis_modes:
            console.print(f"\n[bold magenta]{'='*60}[/bold magenta]")
            console.print(f"[bold magenta]Running analysis mode: {analysis_mode}[/bold magenta]")
            console.print(f"[bold magenta]{'='*60}[/bold magenta]")
            
            if analysis_mode == 'all':
                # All mode: Use all data to calculate one single FF
                console.print(f"\n[bold cyan]Processing all chromosomes together[/bold cyan]")
                console.print(f"Total SNPs: {len(df)}")
                
                try:
                    # Estimate fetal fraction using cfDNA mode (before MQ filtering)
                    console.print(f"[cyan]Estimating FF (before MQ) using all data...[/cyan]")
                    ff_before_mq, _ = ff_estimator_cfdna.estimate(
                        df,
                        f_min=ff_min,
                        f_max=ff_max,
                        f_step=ff_step,
                        ncpus=ncpus,
                        ref_col=cfdna_ref_col,
                        alt_col=cfdna_alt_col
                    )
                    
                    # Estimate fetal fraction using cfDNA+model mode (after MQ filtering)
                    console.print(f"[cyan]Estimating FF (after MQ) using all data...[/cyan]")
                    ff_after_mq, _ = ff_estimator_model.estimate(
                        df,
                        f_min=ff_min,
                        f_max=ff_max,
                        f_step=ff_step,
                        ncpus=ncpus,
                        ref_col=model_ref_col,
                        alt_col=model_alt_col
                    )
                    
                    # Store results with mode-specific name
                    all_results_list.append({
                        'chr': 'all',
                        'ff_before_mq': ff_before_mq,
                        'ff_after_mq': ff_after_mq
                    })
                    
                    console.print(
                        f"[green]✓ All: FF_before = {ff_before_mq:.3f}, "
                        f"FF_after = {ff_after_mq:.3f}[/green]"
                    )
                
                except Exception as e:
                    console.print(f"[red]✗ Error processing all mode: {str(e)}[/red]")
                    if verbose:
                        console.print_exception()
            
            else:
                # chr_only or chr_exclude mode: Process each chromosome
                for target_chr in target_chromosomes:
                    chr_name = f"chr{target_chr}"
                    console.print(f"\n[bold cyan]Processing chromosome {target_chr}[/bold cyan]")
                    
                    # Filter data based on mode
                    target_data = df[df['chr'] == chr_name]
                    background_data = df[df['chr'] != chr_name]
                    
                    if len(target_data) == 0:
                        console.print(f"[yellow]Warning: No SNPs found for {chr_name}, skipping[/yellow]")
                        continue
                    
                    # Select data based on mode
                    if analysis_mode == 'chr_only':
                        estimation_data = target_data
                        data_description = f"Target SNPs only: {len(estimation_data)}"
                        chr_label = f"{chr_name}_only"
                    elif analysis_mode == 'chr_exclude':
                        if len(background_data) == 0:
                            console.print(f"[red]Error: No background SNPs available for {chr_name}[/red]")
                            continue
                        estimation_data = background_data
                        data_description = f"Target SNPs: {len(target_data)}, Background SNPs: {len(background_data)}"
                        chr_label = f"{chr_name}_exclude"
                    
                    console.print(data_description)
                    
                    try:
                        # Estimate fetal fraction using cfDNA mode (before MQ filtering)
                        console.print(f"[cyan]Estimating FF (before MQ) for {chr_name}...[/cyan]")
                        ff_before_mq, _ = ff_estimator_cfdna.estimate(
                            estimation_data,
                            f_min=ff_min,
                            f_max=ff_max,
                            f_step=ff_step,
                            ncpus=ncpus,
                            ref_col=cfdna_ref_col,
                            alt_col=cfdna_alt_col
                        )
                        
                        # Estimate fetal fraction using cfDNA+model mode (after MQ filtering)
                        console.print(f"[cyan]Estimating FF (after MQ) for {chr_name}...[/cyan]")
                        ff_after_mq, _ = ff_estimator_model.estimate(
                            estimation_data,
                            f_min=ff_min,
                            f_max=ff_max,
                            f_step=ff_step,
                            ncpus=ncpus,
                            ref_col=model_ref_col,
                            alt_col=model_alt_col
                        )
                        
                        # Store results with mode-specific chromosome label
                        all_results_list.append({
                            'chr': chr_label,
                            'ff_before_mq': ff_before_mq,
                            'ff_after_mq': ff_after_mq
                        })
                        
                        console.print(
                            f"[green]✓ {chr_label}: FF_before = {ff_before_mq:.3f}, "
                            f"FF_after = {ff_after_mq:.3f}[/green]"
                        )

                    except Exception as e:
                        console.print(f"[red]✗ Error processing {chr_name}: {str(e)}[/red]")
                        if verbose:
                            console.print_exception()
                        continue
        
        # Save all results to a single file
        if all_results_list:
            output_path = Path(f'{output_prefix}_ff.tsv')
            console.print(f"\n[cyan]Saving all results to {output_path}...[/cyan]")
            results_df = pd.DataFrame(all_results_list, columns=['chr', 'ff_before_mq', 'ff_after_mq'])
            results_df.to_csv(output_path, sep='\t', index=False)
            
            # Display summary table
            display_results_summary(results_df, modes=analysis_modes)
            console.print(f"[bold green]✓ Analysis complete! Results saved to {output_path}[/bold green]")
        else:
            console.print("[red]✗ No results generated - check input data and parameters[/red]")
            sys.exit(1)
            
    except Exception as e:
        console.print(f"[red]✗ Fatal error: {str(e)}[/red]")
        if verbose:
            console.print_exception()
        sys.exit(1)


def parse_chromosome_list(chr_spec: str) -> list:
    """
    Parse chromosome specification string into list of chromosome numbers.
    
    Args:
        chr_spec: String like "1-22", "1,2,3", or "21"
    
    Returns:
        List of chromosome numbers (integers)
    
    Raises:
        ValueError: If chromosome specification is invalid
    """
    try:
        if '-' in chr_spec:
            # Range specification (e.g., "1-22")
            start, end = map(int, chr_spec.split('-'))
            return list(range(start, end + 1))
        elif ',' in chr_spec:
            # Comma-separated list (e.g., "1,2,3")
            return [int(x.strip()) for x in chr_spec.split(',')]
        else:
            # Single chromosome (e.g., "21")
            return [int(chr_spec)]
    except ValueError:
        raise ValueError(f"Invalid chromosome specification: {chr_spec}")


def parse_mode_list(mode_spec: str) -> list:
    """
    Parse mode specification string into list of analysis modes.
    
    Args:
        mode_spec: Comma-separated string like "chr_only,chr_exclude,all"
    
    Returns:
        List of mode strings
    
    Raises:
        ValueError: If mode specification contains invalid modes
    """
    valid_modes = {'chr_only', 'chr_exclude', 'all'}
    
    # Parse comma-separated modes
    modes = [m.strip() for m in mode_spec.split(',')]
    
    # Validate each mode
    invalid_modes = set(modes) - valid_modes
    if invalid_modes:
        raise ValueError(
            f"Invalid mode(s): {invalid_modes}. "
            f"Valid modes are: {', '.join(valid_modes)}"
        )
    
    # Remove duplicates while preserving order
    seen = set()
    unique_modes = []
    for mode in modes:
        if mode not in seen:
            seen.add(mode)
            unique_modes.append(mode)
    
    return unique_modes


def load_and_validate_data(
    input_path: Path,
    cfdna_ref_col: str = 'cfDNA_ref_reads',
    cfdna_alt_col: str = 'cfDNA_alt_reads',
    model_ref_col: str = 'fetal_ref_reads_from_model',
    model_alt_col: str = 'fetal_alt_reads_from_model',
    min_raw_depth: int = 0,
    min_model_depth: int = 0
) -> pd.DataFrame:
    """
    Load and validate input SNP data from TSV.GZ file.
    
    This function validates the presence of required columns for both cfDNA and
    cfDNA+model analysis modes. It performs data type validation, cleaning, and
    applies depth filters.
    
    Args:
        input_path (Path): Path to input file
        cfdna_ref_col (str): Column name for cfDNA reference reads
        cfdna_alt_col (str): Column name for cfDNA alternative reads
        model_ref_col (str): Column name for modeled fetal reference reads
        model_alt_col (str): Column name for modeled fetal alternative reads
        min_raw_depth (int): Minimum depth for cfDNA reads (default: 0)
        min_model_depth (int): Minimum depth for modeled reads (default: 0)
    
    Returns:
        pd.DataFrame: Validated pandas DataFrame with proper data types
    
    Raises:
        ValueError: If data validation fails or required columns are missing
        FileNotFoundError: If input file doesn't exist
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    try:
        # Load data with appropriate compression detection
        df = pd.read_csv(input_path, sep='\t', compression='gzip' if input_path.suffix == '.gz' else None)
    except Exception as e:
        raise ValueError(f"Failed to load input file: {str(e)}")
    
    # Determine required columns for both modes
    base_columns = {'chr', 'pos', 'af'}
    required_columns = base_columns.copy()
    required_columns.update({cfdna_ref_col, cfdna_alt_col, model_ref_col, model_alt_col})
    read_columns = [cfdna_ref_col, cfdna_alt_col, model_ref_col, model_alt_col]
    
    console.print(f"[cyan]Validating columns for both cfDNA and cfDNA+model modes[/cyan]")
    console.print(f"[cyan]cfDNA columns: {cfdna_ref_col}, {cfdna_alt_col}[/cyan]")
    console.print(f"[cyan]Model columns: {model_ref_col}, {model_alt_col}[/cyan]")
    
    # Check for required columns
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")
    
    # Basic data validation
    if len(df) == 0:
        raise ValueError("Input file is empty")
    
    # Ensure proper data types for base columns
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    df['af'] = pd.to_numeric(df['af'], errors='coerce')
    
    # Ensure proper data types for read count columns
    for col in read_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce', downcast='integer')
    
    # Build invalid mask for base columns
    invalid_mask = (
        df['af'].isna() | (df['af'] < 0) | (df['af'] > 1)
    )
    
    # Add read count validation to invalid mask
    for col in read_columns:
        invalid_mask |= (df[col].isna() | (df[col] < 0))
    
    if invalid_mask.any():
        console.print(f"[yellow]Warning: Removing {invalid_mask.sum()} rows with invalid data[/yellow]")
        df = df[~invalid_mask]
    
    if len(df) == 0:
        raise ValueError("No valid data remaining after filtering")
    
    # Apply depth filters
    initial_count = len(df)
    
    # Filter by minimum raw depth (cfDNA reads)
    if min_raw_depth > 0:
        df['cfdna_total_depth'] = df[cfdna_ref_col] + df[cfdna_alt_col]
        df = df[df['cfdna_total_depth'] >= min_raw_depth]
        remaining_count = len(df)
        percentage = (remaining_count / initial_count) * 100 if initial_count > 0 else 0
        console.print(f"[cyan]Raw depth filter (≥{min_raw_depth}): {remaining_count:,} SNPs remaining ({percentage:.1f}%)[/cyan]")
        df = df.drop(columns=['cfdna_total_depth']).reset_index(drop=True)
        
        if len(df) == 0:
            raise ValueError("No SNPs remaining after raw depth filtering")
    
    # Filter by minimum model depth
    if min_model_depth > 0:
        current_count = len(df)
        df['model_total_depth'] = df[model_ref_col] + df[model_alt_col]
        df = df[df['model_total_depth'] >= min_model_depth]
        remaining_count = len(df)
        percentage = (remaining_count / current_count) * 100 if current_count > 0 else 0
        console.print(f"[cyan]Model depth filter (≥{min_model_depth}): {remaining_count:,} SNPs remaining ({percentage:.1f}%)[/cyan]")
        df = df.drop(columns=['model_total_depth']).reset_index(drop=True)
        
        if len(df) == 0:
            raise ValueError("No SNPs remaining after model depth filtering")
    
    return df


def display_results_summary(results_df: pd.DataFrame, modes: list = None) -> None:
    """
    Display a formatted summary table of fetal fraction estimation results.
    
    Args:
        results_df: DataFrame containing FF results with columns 'chr', 'ff_before_mq', 'ff_after_mq'.
        modes: List of analysis modes used (e.g., ['chr_only', 'chr_exclude', 'all'])
    """
    if modes is None:
        modes = ['chr_exclude']
    
    mode_str = ', '.join(modes) if len(modes) > 1 else modes[0]
    table = Table(title=f"Fetal Fraction Estimation Results - Modes: {mode_str}")
    table.add_column("Chromosome", justify="center", style="cyan")
    table.add_column("FF Before MQ", justify="right", style="green")
    table.add_column("FF After MQ", justify="right", style="yellow")
    
    for _, row in results_df.iterrows():
        chr_name = row['chr']
        ff_before = row['ff_before_mq']
        ff_after = row['ff_after_mq']
        
        # Use the chromosome name as-is (already includes mode suffix)
        display_name = str(chr_name)
        
        table.add_row(
            display_name,
            f"{ff_before:.3f}",
            f"{ff_after:.3f}"
        )
    
    # Add summary statistics for each mode if multiple rows exist
    if len(results_df) > 1:
        # Calculate mean for each mode
        for mode in modes:
            if mode == 'all':
                # For 'all' mode, just show that single value (already displayed)
                continue
            else:
                # Filter rows for this mode
                mode_suffix = f"_{mode}"
                mode_rows = results_df[results_df['chr'].str.endswith(mode_suffix)]
                
                if len(mode_rows) > 0:
                    mean_before = mode_rows['ff_before_mq'].mean()
                    mean_after = mode_rows['ff_after_mq'].mean()
                    
                    table.add_row(
                        f"[bold]Mean ({mode})[/bold]",
                        f"[bold]{mean_before:.3f}[/bold]",
                        f"[bold]{mean_after:.3f}[/bold]"
                    )
    
    console.print(table)


if __name__ == "__main__":
    # Set multiprocessing start method for compatibility
    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass  # Start method already set
    
    main()
