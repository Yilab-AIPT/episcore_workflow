#!/usr/bin/env python3
"""
BAM to pileup converter for WGS (non-bisulfite) sequencing.

Generates raw reference/alternate allele counts at known SNP sites from a WGS BAM.
Unlike ``bin/bam_to_pileup.py``, this script does not apply bisulfite strand
filtering or methylation-aware base classification.
"""

import sys
import gzip
from pathlib import Path
from typing import Optional, Tuple, List
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed

import click
import pandas as pd
import pysam
from rich.console import Console
from rich.progress import Progress, TaskID
from rich.table import Table

console = Console()


@dataclass
class SNPSite:
    chr: str
    pos: int
    ref: str
    alt: str
    af: float


@dataclass
class PileupCounts:
    cfDNA_ref_reads: int = 0
    cfDNA_alt_reads: int = 0

    @property
    def current_depth(self) -> int:
        return self.cfDNA_ref_reads + self.cfDNA_alt_reads


def parse_known_sites(
    sites_file: Path, progress: Progress, task_id: TaskID, max_sites: Optional[int] = None
) -> List[SNPSite]:
    if not sites_file.exists():
        raise FileNotFoundError(f"Known sites file not found: {sites_file}")

    progress.update(task_id, description="Parsing known sites...")

    sites_data = pd.read_csv(
        sites_file,
        sep='\t',
        comment='#',
        usecols=[0, 1, 3, 4, 7],
        names=['chr', 'pos', 'ref', 'alt', 'info_field'],
        dtype={'chr': str, 'pos': int, 'ref': str, 'alt': str, 'info_field': str},
        nrows=max_sites,
    )

    progress.update(task_id, advance=25)

    sites_data['ref'] = sites_data['ref'].str.upper()
    sites_data['alt'] = sites_data['alt'].str.upper()

    single_nuc_mask = (sites_data['ref'].str.len() == 1) & (sites_data['alt'].str.len() == 1)
    sites_data = sites_data[single_nuc_mask]

    progress.update(task_id, advance=25)

    af_pattern = r'AF=([^;]+)'
    sites_data['af_match'] = sites_data['info_field'].str.extract(af_pattern, expand=False)
    sites_data = sites_data.dropna(subset=['af_match'])
    sites_data['af'] = pd.to_numeric(sites_data['af_match'], errors='coerce')
    sites_data = sites_data.dropna(subset=['af'])

    progress.update(task_id, advance=25)

    sites = [
        SNPSite(row['chr'], row['pos'], row['ref'], row['alt'], row['af'])
        for _, row in sites_data.iterrows()
    ]

    progress.update(task_id, advance=25)
    console.print(f"[green]✓[/green] Parsed {len(sites):,} single-nucleotide SNP sites")

    return sites


def classify_base(base: str, ref: str, alt: str) -> Optional[str]:
    base = base.upper()
    ref = ref.upper()
    alt = alt.upper()
    if base == ref:
        return 'REF'
    if base == alt:
        return 'ALT'
    return None


def process_pileup_site(
    bam_file: pysam.AlignmentFile, site: SNPSite, min_mapq: int, min_bq: int
) -> PileupCounts:
    counts = PileupCounts()
    processed_templates = set()

    try:
        for pileup_column in bam_file.pileup(
            site.chr,
            site.pos - 1,
            site.pos,
            stepper='nofilter',
            ignore_overlaps=True,
            ignore_orphans=True,
            min_base_quality=min_bq,
            min_mapping_quality=min_mapq,
        ):
            if pileup_column.pos != site.pos - 1:
                continue

            for pileup_read in pileup_column.pileups:
                read = pileup_read.alignment

                if read.mapping_quality < min_mapq:
                    continue
                if read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue
                if read.is_unmapped:
                    continue
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                if pileup_read.query_position is None:
                    continue

                query_base = read.query_sequence[pileup_read.query_position]
                base_quality = read.query_qualities[pileup_read.query_position]
                if base_quality < min_bq:
                    continue

                template_key = (read.query_name, site.chr, site.pos)
                if template_key in processed_templates:
                    continue
                processed_templates.add(template_key)

                classification = classify_base(query_base, site.ref, site.alt)
                if classification == 'REF':
                    counts.cfDNA_ref_reads += 1
                elif classification == 'ALT':
                    counts.cfDNA_alt_reads += 1

            break

    except Exception as e:
        console.print(f"[yellow]Warning: Error processing site {site.chr}:{site.pos}: {e}[/yellow]")

    return counts


def process_sites_chunk(
    bam_path: Path, sites_chunk: List[SNPSite], min_mapq: int, min_bq: int
) -> List[Tuple[SNPSite, PileupCounts]]:
    results = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for site in sites_chunk:
            counts = process_pileup_site(bam, site, min_mapq, min_bq)
            results.append((site, counts))
    return results


def generate_pileup_data(
    bam_file: Path,
    sites: List[SNPSite],
    min_mapq: int,
    min_bq: int,
    ncpus: int,
    progress: Progress,
    task_id: TaskID,
) -> List[Tuple[SNPSite, PileupCounts]]:
    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    if not sites:
        console.print("[yellow]Warning: No sites to process[/yellow]")
        return []

    progress.update(task_id, description=f"Processing pileup data with {ncpus} workers...")

    min_chunk_size = 50
    target_chunks = ncpus * 3
    chunk_size = max(min_chunk_size, len(sites) // target_chunks)
    site_chunks = [sites[i:i + chunk_size] for i in range(0, len(sites), chunk_size)]

    console.print(
        f"[blue]Processing {len(sites):,} sites in {len(site_chunks)} chunks "
        f"(~{chunk_size} sites/chunk) using {ncpus} workers[/blue]"
    )

    results = []
    completed_sites = 0

    with ProcessPoolExecutor(max_workers=ncpus) as executor:
        future_to_chunk = {
            executor.submit(process_sites_chunk, bam_file, chunk, min_mapq, min_bq): chunk
            for chunk in site_chunks
        }
        for future in as_completed(future_to_chunk):
            chunk = future_to_chunk[future]
            chunk_results = future.result()
            results.extend(chunk_results)
            completed_sites += len(chunk)
            progress_pct = (completed_sites / len(sites)) * 100
            progress.update(task_id, completed=progress_pct)

    progress.update(task_id, completed=100)
    console.print(f"[green]✓[/green] Processed pileup data for {len(results):,} sites")
    return results


def save_pileup_output(
    results: List[Tuple[SNPSite, PileupCounts]], output_prefix: str, progress: Progress, task_id: TaskID
) -> Path:
    output_file = Path(f"{output_prefix}_pileup.tsv.gz")
    progress.update(task_id, description="Saving pileup data...")

    def sort_key(item):
        site, _ = item
        chr_name = site.chr
        chr_part = chr_name[3:] if chr_name.startswith('chr') else chr_name
        try:
            return (0, int(chr_part), site.pos)
        except ValueError:
            return (1, chr_part, site.pos)

    sorted_results = sorted(results, key=sort_key)

    with gzip.open(output_file, 'wt') as f:
        header = [
            'chr', 'pos', 'ref', 'alt', 'af',
            'cfDNA_ref_reads', 'cfDNA_alt_reads', 'current_depth',
        ]
        f.write('\t'.join(header) + '\n')
        for site, counts in sorted_results:
            row = [
                site.chr,
                str(site.pos),
                site.ref,
                site.alt,
                f"{site.af:.6f}",
                str(counts.cfDNA_ref_reads),
                str(counts.cfDNA_alt_reads),
                str(counts.current_depth),
            ]
            f.write('\t'.join(row) + '\n')

    progress.update(task_id, advance=100)
    console.print(f"[green]✓[/green] Pileup data saved to: {output_file}")
    return output_file


@click.command()
@click.option('--input-bam', required=True, type=click.Path(exists=True, path_type=Path))
@click.option('--known-sites', required=True, type=click.Path(exists=True, path_type=Path))
@click.option('--output', required=True, type=str, help='Output prefix ({prefix}_pileup.tsv.gz)')
@click.option('--min-mapq', default=20, type=int, show_default=True)
@click.option('--min-bq', default=13, type=int, show_default=True)
@click.option('--ncpus', default=8, type=int, show_default=True)
@click.option(
    '--max-sites',
    default=None,
    type=int,
    help='Optional cap on number of known sites (for quick tests)',
)
def main(
    input_bam: Path,
    known_sites: Path,
    output: str,
    min_mapq: int,
    min_bq: int,
    ncpus: int,
    max_sites: Optional[int],
) -> None:
    console.print("\n[bold blue]WGS BAM to Pileup Converter[/bold blue]")
    console.print("=" * 70)

    params_table = Table(title="Input Parameters", show_header=True, header_style="bold magenta")
    params_table.add_column("Parameter", style="cyan", no_wrap=True)
    params_table.add_column("Value", style="white")
    params_table.add_row("Input BAM", str(input_bam))
    params_table.add_row("Known Sites", str(known_sites))
    params_table.add_row("Output Prefix", output)
    params_table.add_row("Min MAPQ", str(min_mapq))
    params_table.add_row("Min Base Quality", str(min_bq))
    params_table.add_row("Parallel Workers", str(ncpus))
    if max_sites is not None:
        params_table.add_row("Max Sites", str(max_sites))
    console.print(params_table)
    console.print()

    try:
        with Progress(console=console) as progress:
            sites_task = progress.add_task("Parsing known sites...", total=100)
            pileup_task = progress.add_task("Processing pileup...", total=100)
            save_task = progress.add_task("Saving output...", total=100)

            all_sites = parse_known_sites(known_sites, progress, sites_task, max_sites=max_sites)
            pileup_results = generate_pileup_data(
                input_bam, all_sites, min_mapq, min_bq, ncpus, progress, pileup_task
            )
            output_file = save_pileup_output(pileup_results, output, progress, save_task)

        total_sites = len(pileup_results)
        sites_with_coverage = sum(1 for _, counts in pileup_results if counts.current_depth > 0)
        mean_depth = (
            sum(counts.current_depth for _, counts in pileup_results) / total_sites
            if total_sites else 0.0
        )

        summary_table = Table(title="Processing Summary", show_header=True, header_style="bold green")
        summary_table.add_column("Metric", style="cyan", no_wrap=True)
        summary_table.add_column("Count", style="white", justify="right")
        summary_table.add_row("Total SNP sites processed", f"{total_sites:,}")
        summary_table.add_row("Sites with coverage", f"{sites_with_coverage:,}")
        summary_table.add_row("Mean depth", f"{mean_depth:.2f}")
        console.print(summary_table)
        console.print(f"\n[bold green]✓ Processing completed successfully![/bold green]")
        console.print(f"Output file: [cyan]{output_file}[/cyan]\n")

    except Exception as e:
        console.print(f"\n[bold red]✗ Error during processing:[/bold red] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
