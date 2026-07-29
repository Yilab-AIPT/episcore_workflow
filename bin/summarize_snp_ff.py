#!/usr/bin/env python3
"""Summarize per-sample SNP fetal-fraction TSV files into one table."""

from pathlib import Path

import click
import pandas as pd
from rich.console import Console

console = Console()


def _sample_id_from_ff_path(ff_path: Path) -> str:
    name = ff_path.name
    if name.endswith('_ff.tsv'):
        return name[:-len('_ff.tsv')]
    return ff_path.stem


def _summarize_ff_file(ff_path: Path) -> dict:
    df = pd.read_csv(ff_path, sep='\t')
    required = {'ff_before_mq', 'ff_after_mq'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{ff_path}: missing columns {sorted(missing)}")

    return {
        'sample': _sample_id_from_ff_path(ff_path),
        'ff_before_mq': round(float(df['ff_before_mq'].mean()), 6),
        'ff_after_mq': round(float(df['ff_after_mq'].mean()), 6),
    }


@click.command()
@click.option(
    '--input-dir',
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help='Directory containing per-sample *_ff.tsv files',
)
@click.option(
    '--output',
    '-o',
    type=click.Path(path_type=Path),
    required=True,
    help='Output summary TSV path',
)
def main(input_dir: Path, output: Path) -> None:
    ff_files = sorted(input_dir.glob('*_ff.tsv'))
    if not ff_files:
        console.print(f"[red]No *_ff.tsv files found in {input_dir}[/red]")
        raise SystemExit(1)

    rows = [_summarize_ff_file(ff_path) for ff_path in ff_files]
    summary = pd.DataFrame(rows, columns=['sample', 'ff_before_mq', 'ff_after_mq'])
    summary['ff'] = summary['ff_after_mq']
    summary = summary[['sample', 'ff', 'ff_before_mq', 'ff_after_mq']].sort_values('sample')

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, sep='\t', index=False)
    console.print(f"[green]Wrote FF summary for {len(summary)} samples to {output}[/green]")


if __name__ == '__main__':
    main()
