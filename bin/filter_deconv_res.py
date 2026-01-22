#!/usr/bin/env python3
"""Filter deconvolution results and split reads by probability threshold.

This script processes deconvolution result files (txt/parquet format) and splits
read names into target and background sets based on a probability threshold.
The output files can be used with samtools to split BAM files accordingly.
"""

import sys
from pathlib import Path
from typing import Set, Tuple

import click
import pandas as pd
import pyarrow.parquet as pq
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

console = Console()


def detect_file_format(input_path: Path) -> str:
    """Detect whether the input file is txt (TSV) or parquet format.
    
    Args:
        input_path: Path to the input file.
        
    Returns:
        File format as string: 'txt' or 'parquet'.
        
    Raises:
        ValueError: If file format cannot be determined or is unsupported.
    """
    suffix = input_path.suffix.lower()
    
    if suffix in ['.txt', '.tsv', '.csv']:
        return 'txt'
    elif suffix in ['.parquet', '.pq']:
        return 'parquet'
    else:
        # Try to detect by reading first few bytes
        try:
            with open(input_path, 'rb') as f:
                magic = f.read(4)
                if magic == b'PAR1':
                    return 'parquet'
        except Exception:
            pass
        
        # Default to txt if unsure
        console.print(f"[yellow]Warning: Cannot determine file format from extension '{suffix}', assuming txt/tsv format[/yellow]")
        return 'txt'


def validate_columns(df: pd.DataFrame, required_cols: list) -> None:
    """Validate that required columns exist in the dataframe.
    
    Args:
        df: DataFrame to validate.
        required_cols: List of required column names.
        
    Raises:
        ValueError: If any required column is missing.
    """
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Required column(s) {missing_cols} not found in input file. "
            f"Available columns: {list(df.columns)}"
        )


def read_txt_in_chunks(input_path: Path, chunksize: int = 100000):
    """Read txt/tsv file in chunks for memory efficiency.
    
    Args:
        input_path: Path to the txt/tsv file.
        chunksize: Number of rows to read per chunk.
        
    Yields:
        DataFrame chunks containing only 'name' and 'prob_class_1' columns.
        
    Raises:
        ValueError: If required columns are missing.
        pd.errors.ParserError: If file parsing fails.
    """
    # Read first chunk to validate columns
    first_chunk = pd.read_csv(
        input_path,
        sep='\t',
        nrows=1,
        dtype=str
    )
    validate_columns(first_chunk, ['name', 'prob_class_1'])
    
    # Now read in chunks, selecting only required columns
    for chunk in pd.read_csv(
        input_path,
        sep='\t',
        usecols=['name', 'prob_class_1'],
        dtype={'name': str, 'prob_class_1': float},
        chunksize=chunksize
    ):
        yield chunk


def read_parquet_in_batches(input_path: Path, batch_size: int = 100000):
    """Read parquet file in batches for memory efficiency.
    
    Args:
        input_path: Path to the parquet file.
        batch_size: Number of rows to read per batch.
        
    Yields:
        DataFrame batches containing only 'name' and 'prob_class_1' columns.
        
    Raises:
        ValueError: If required columns are missing.
        Exception: If parquet reading fails.
    """
    parquet_file = pq.ParquetFile(input_path)
    
    # Validate columns exist
    schema_names = parquet_file.schema.names
    required_cols = ['name', 'prob_class_1']
    missing_cols = [col for col in required_cols if col not in schema_names]
    if missing_cols:
        raise ValueError(
            f"Required column(s) {missing_cols} not found in parquet file. "
            f"Available columns: {schema_names}"
        )
    
    # Read in batches, selecting only required columns
    for batch in parquet_file.iter_batches(batch_size=batch_size, columns=required_cols):
        yield batch.to_pandas()


def process_chunk(
    chunk: pd.DataFrame,
    threshold: float,
    target_set: Set[str],
    background_set: Set[str],
    classified_set: Set[str]
) -> Tuple[int, int, int]:
    """Process a chunk of data and update read name sets.
    
    Args:
        chunk: DataFrame chunk with 'name' and 'prob_class_1' columns.
        threshold: Probability threshold for filtering.
        target_set: Set to accumulate target read names (prob > threshold).
        background_set: Set to accumulate background read names (prob <= threshold).
        classified_set: Set to accumulate all classified read names.
        
    Returns:
        Tuple of (target_count, background_count, total_count) for this chunk.
    """
    # Remove any rows with missing values
    chunk = chunk.dropna()
    
    # Convert prob_class_1 to float if not already
    chunk['prob_class_1'] = pd.to_numeric(chunk['prob_class_1'], errors='coerce')
    chunk = chunk.dropna()
    
    # Filter by threshold
    target_mask = chunk['prob_class_1'] > threshold
    target_names = chunk.loc[target_mask, 'name'].tolist()
    background_names = chunk.loc[~target_mask, 'name'].tolist()
    all_names = chunk['name'].tolist()
    
    # Update sets (automatically handles deduplication)
    target_set.update(target_names)
    background_set.update(background_names)
    classified_set.update(all_names)
    
    return len(target_names), len(background_names), len(all_names)


def write_read_names(output_path: Path, read_names: Set[str]) -> None:
    """Write read names to a file, one per line, sorted.
    
    Args:
        output_path: Path to the output file.
        read_names: Set of read names to write.
    """
    with open(output_path, 'w') as f:
        for name in sorted(read_names):
            f.write(f"{name}\n")


@click.command()
@click.option(
    '--input',
    'input_file',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Input file path (txt/tsv or parquet format) with columns: name, prob_class_1'
)
@click.option(
    '--threshold',
    type=float,
    required=True,
    help='Probability threshold to filter reads (prob_class_1 > threshold -> target)'
)
@click.option(
    '--output-dir',
    type=click.Path(path_type=Path),
    default=Path.cwd(),
    help='Output directory for result files (default: current directory)'
)
@click.option(
    '--chunksize',
    type=int,
    default=100000,
    help='Number of rows to process per chunk (default: 100000)'
)
def main(input_file: Path, threshold: float, output_dir: Path, chunksize: int):
    """Filter deconvolution results and split reads by probability threshold.
    
    This script reads a deconvolution result file (txt/tsv or parquet format),
    filters reads based on a probability threshold, and outputs three files:
    
    \b
    - target_reads.txt: Read names with prob_class_1 > threshold
    - background_reads.txt: Read names with prob_class_1 <= threshold  
    - classified_reads.txt: All read names
    
    The output files contain unique, sorted read names (one per line) suitable
    for use with 'samtools view -N' to split BAM files.
    
    Examples:
        filter_deconv_res.py --input deconv.txt --threshold 0.5
        
        filter_deconv_res.py --input results.parquet --threshold 0.7 --output-dir ./output
    """
    try:
        # Display header
        console.print("\n[bold cyan]Deconvolution Results Filter[/bold cyan]")
        console.print(f"Input file: {input_file}")
        console.print(f"Threshold: {threshold}")
        console.print(f"Output directory: {output_dir}\n")
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Detect file format
        file_format = detect_file_format(input_file)
        console.print(f"[green]Detected file format: {file_format}[/green]")
        
        # Initialize sets for deduplication
        target_set: Set[str] = set()
        background_set: Set[str] = set()
        classified_set: Set[str] = set()
        
        # Counters
        total_target = 0
        total_background = 0
        total_classified = 0
        chunks_processed = 0
        
        # Process file in chunks
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console
        ) as progress:
            
            task = progress.add_task("[cyan]Processing chunks...", total=None)
            
            try:
                if file_format == 'txt':
                    chunk_iterator = read_txt_in_chunks(input_file, chunksize)
                else:  # parquet
                    chunk_iterator = read_parquet_in_batches(input_file, chunksize)
                
                for chunk in chunk_iterator:
                    target_count, background_count, classified_count = process_chunk(
                        chunk, threshold, target_set, background_set, classified_set
                    )
                    
                    total_target += target_count
                    total_background += background_count
                    total_classified += classified_count
                    chunks_processed += 1
                    
                    progress.update(
                        task,
                        description=f"[cyan]Processed {chunks_processed} chunks "
                                  f"({total_classified:,} reads total)..."
                    )
                
                progress.update(task, description="[green]Processing complete!")
                
            except ValueError as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
                sys.exit(1)
            except Exception as e:
                console.print(f"[bold red]Unexpected error during processing:[/bold red] {e}")
                sys.exit(1)
        
        # Write output files
        console.print("\n[cyan]Writing output files...[/cyan]")
        
        output_files = {
            'target_reads.txt': target_set,
            'background_reads.txt': background_set,
            'classified_reads.txt': classified_set
        }
        
        for filename, read_set in output_files.items():
            output_path = output_dir / filename
            write_read_names(output_path, read_set)
            console.print(f"  [green]✓[/green] {filename}: {len(read_set):,} unique reads")
        
        # Summary statistics
        console.print("\n[bold green]Summary:[/bold green]")
        console.print(f"  Total chunks processed: {chunks_processed:,}")
        console.print(f"  Total reads processed: {total_classified:,}")
        console.print(f"  Unique classified reads: {len(classified_set):,}")
        console.print(f"  Unique target reads (prob > {threshold}): {len(target_set):,} "
                     f"({len(target_set)/len(classified_set)*100:.2f}%)")
        console.print(f"  Unique background reads (prob ≤ {threshold}): {len(background_set):,} "
                     f"({len(background_set)/len(classified_set)*100:.2f}%)")
        
        console.print("\n[bold green]✓ Processing completed successfully![/bold green]\n")
        
    except KeyboardInterrupt:
        console.print("\n[yellow]Operation cancelled by user[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[bold red]Fatal error:[/bold red] {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
