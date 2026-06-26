#!/usr/bin/env python3
"""
Generate Final Report for NIPT Analysis

This script consolidates results from episcore analysis and SNP-based fetal fraction
estimation to create a comprehensive report with the following columns:
- sample: Sample identifier
- ff_before_mq: Average fetal fraction before MethylQueen filtering (across chromosomes)
- ff_after_mq: Average fetal fraction after MethylQueen filtering (across chromosomes)
- cpg_mean_coverage: Mean CpG coverage (average of raw_total_count)
- snp_mean_coverage: Mean SNP coverage (average of current_depth)
- chr{#}_s_inter: Inter-sample Z-score for each chromosome (from episcore analysis)
"""

import pandas as pd
import gzip
import click
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
from typing import Optional
import sys

# Initialize rich console for beautiful output
console = Console()


def read_episcore_file(episcore_file: Path) -> pd.DataFrame:
    """Read episcore TSV file and extract s_inter columns.
    
    Args:
        episcore_file: Path to episcore TSV file
        
    Returns:
        DataFrame with s_inter columns for each chromosome
    """
    console.print(f"[cyan]Reading episcore file: {episcore_file}[/cyan]")
    df = pd.read_csv(episcore_file, sep='\t')
    
    # Extract only the s_inter columns (chr{#}_s_inter)
    s_inter_cols = [col for col in df.columns if col.endswith('_s_inter')]
    
    if not s_inter_cols:
        console.print("[yellow]Warning: No s_inter columns found in episcore file[/yellow]")
        return pd.DataFrame()
    
    console.print(f"[green]✓ Found {len(s_inter_cols)} s_inter columns[/green]")
    return df[s_inter_cols]


def read_beta_value_file(beta_file: Path) -> float:
    """Read beta value gzipped TSV and calculate mean raw_total_count.
    
    Args:
        beta_file: Path to beta value gzipped TSV file
        
    Returns:
        Mean raw_total_count value
    """
    console.print(f"[cyan]Reading beta value file: {beta_file}[/cyan]")
    
    with gzip.open(beta_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    
    if 'raw_total_count' not in df.columns:
        console.print("[red]Error: 'raw_total_count' column not found in beta value file[/red]")
        return 0.0
    
    mean_coverage = df['raw_total_count'].mean()
    console.print(f"[green]✓ Calculated CpG mean coverage: {mean_coverage:.2f}[/green]")
    return mean_coverage


def read_snp_pileup_file(pileup_file: Path) -> float:
    """Read SNP pileup gzipped TSV and calculate mean current_depth.
    
    Args:
        pileup_file: Path to SNP pileup gzipped TSV file
        
    Returns:
        Mean current_depth value
    """
    console.print(f"[cyan]Reading SNP pileup file: {pileup_file}[/cyan]")
    
    with gzip.open(pileup_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    
    if 'current_depth' not in df.columns:
        console.print("[red]Error: 'current_depth' column not found in pileup file[/red]")
        return 0.0
    
    mean_coverage = df['current_depth'].mean()
    console.print(f"[green]✓ Calculated SNP mean coverage: {mean_coverage:.2f}[/green]")
    return mean_coverage


def read_snp_ff_file(ff_file: Path) -> tuple:
    """Read SNP fetal fraction TSV and calculate average ff values.
    
    Args:
        ff_file: Path to SNP fetal fraction TSV file
        
    Returns:
        Tuple of (avg_ff_before_mq, avg_ff_after_mq)
    """
    console.print(f"[cyan]Reading SNP FF file: {ff_file}[/cyan]")
    df = pd.read_csv(ff_file, sep='\t')
    
    if 'ff_before_mq' not in df.columns or 'ff_after_mq' not in df.columns:
        console.print("[red]Error: Required FF columns not found in file[/red]")
        return 0.0, 0.0
    
    avg_ff_before = df['ff_before_mq'].mean()
    avg_ff_after = df['ff_after_mq'].mean()
    
    console.print(f"[green]✓ Calculated FF before MQ: {avg_ff_before:.4f}[/green]")
    console.print(f"[green]✓ Calculated FF after MQ: {avg_ff_after:.4f}[/green]")
    
    return avg_ff_before, avg_ff_after


def read_meta_file(meta_file: Path) -> pd.DataFrame:
    """Read meta information TSV file.
    
    Args:
        meta_file: Path to meta TSV file
        
    Returns:
        DataFrame with meta information
    """
    console.print(f"[cyan]Reading meta file: {meta_file}[/cyan]")
    df = pd.read_csv(meta_file, sep='\t')
    
    if 'sample' not in df.columns:
        console.print("[red]Error: 'sample' column not found in meta file[/red]")
        sys.exit(1)
    
    console.print(f"[green]✓ Loaded meta data with {len(df)} rows and {len(df.columns)} columns[/green]")
    return df


@click.command()
@click.option(
    '--sample-id',
    required=True,
    type=str,
    help='Sample identifier'
)
@click.option(
    '--episcore',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to episcore TSV file'
)
@click.option(
    '--beta-value',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to beta value gzipped TSV file'
)
@click.option(
    '--snp-pileup',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to SNP pileup gzipped TSV file'
)
@click.option(
    '--snp-ff',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to SNP fetal fraction TSV file'
)
@click.option(
    '--meta',
    required=False,
    type=str,
    help='Meta information TSV file'
)
@click.option(
    '--output-prefix',
    required=True,
    type=click.Path(path_type=Path),
    help='Output prefix for final report TSV'
)
def main(
    sample_id: str,
    episcore: Path,
    beta_value: Path,
    snp_pileup: Path,
    snp_ff: Path,
    output_prefix: Path,
    meta: Optional[str]
):
    """Generate final NIPT analysis report."""
    try:
        # Print header
        console.print(Panel.fit(
            f"[bold cyan]NIPT Report Generation[/bold cyan]\n"
            f"Sample: {sample_id}",
            border_style="cyan"
        ))
        
        # Read all input files and extract metrics
        console.print("\n[bold]Step 1: Processing input files[/bold]")
        
        # 1. Read episcore file
        episcore_df = read_episcore_file(episcore)
        
        # 2. Read beta value file
        cpg_mean_coverage = read_beta_value_file(beta_value)
        
        # 3. Read SNP pileup file
        snp_mean_coverage = read_snp_pileup_file(snp_pileup)
        
        # 4. Read SNP FF file
        avg_ff_before, avg_ff_after = read_snp_ff_file(snp_ff)
        
        # Build output DataFrame
        console.print("\n[bold]Step 2: Building output report[/bold]")
        
        output_data = {
            'sample': sample_id,
            'ff_before_mq': avg_ff_before,
            'ff_after_mq': avg_ff_after,
            'cpg_mean_coverage': cpg_mean_coverage,
            'snp_mean_coverage': snp_mean_coverage
        }
        
        # Add s_inter columns from episcore
        if not episcore_df.empty:
            for col in episcore_df.columns:
                output_data[col] = episcore_df[col].iloc[0]
        
        # Create single-row DataFrame
        report_df = pd.DataFrame([output_data])
        
        # Merge with meta data if provided
        if meta:
            console.print("\n[bold]Step 2.5: Merging with meta data[/bold]")
            meta_df = read_meta_file(Path(meta))
            # Left join: keep all rows from report_df, add columns from meta_df
            report_df = report_df.merge(meta_df, on='sample', how='left')
            console.print(f"[green]✓ Merged report with meta data[/green]")
            console.print(f"[green]✓ Total columns after merge: {len(report_df.columns)}[/green]")
        
        # Write output
        console.print("\n[bold]Step 3: Writing output file[/bold]")
        output = Path(f"{output_prefix}_report.tsv")
        output.parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(output, sep='\t', index=False, float_format='%.6f')
        
        console.print(f"[green]✓ Report saved to: {output}[/green]")
        console.print(f"[green]✓ Total columns: {len(report_df.columns)}[/green]")
        
        # Display summary
        console.print("\n[bold cyan]Summary Statistics:[/bold cyan]")
        console.print(f"  FF before MQ:      {avg_ff_before:.3f}")
        console.print(f"  FF after MQ:       {avg_ff_after:.3f}")
        console.print(f"  CpG coverage:      {cpg_mean_coverage:.2f}x")
        console.print(f"  SNP coverage:      {snp_mean_coverage:.2f}x")
        console.print(f"  Chromosomes:       {len([c for c in report_df.columns if c.endswith('_s_inter')])}")
        
        console.print("\n[bold green]✓ Report generation complete![/bold green]")
        
    except Exception as e:
        console.print(f"\n[bold red]Error:[/bold red] {str(e)}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == '__main__':
    main()
