#!/usr/bin/env python3
"""Combine per-chromosome episcore, zscore, ezscore and optional FF for one sample.

``ezscore`` is the z-normalised ``episcore + zscore`` using the per-chromosome
mean/std from ``best_ezscore_ref_20_matrix.tsv`` (the fixed ezscore reference).
When ``--ff`` is provided, ``ff_before_mq`` is read from the SNP fetal-fraction
table (the ``chr == 'all'`` row). With ``--skip-ff``, ``ff_before_mq`` is omitted.

Output ``{output_prefix}_scores.tsv`` (long, one row per autosome) has columns:
``sample, [ff_before_mq,] chr, episcore, zscore`` and, unless ``--skip-ezscore``
is set, ``ezscore``.
"""

import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]


def read_ff_before_mq(path: Path) -> float:
    """Return ``ff_before_mq`` from the ``chr == 'all'`` row (fallback: mean)."""
    df = pd.read_csv(path, sep="\t")
    if "ff_before_mq" not in df.columns:
        raise ValueError(f"FF table missing 'ff_before_mq' column: {path}")
    if "chr" in df.columns and (df["chr"].astype(str) == "all").any():
        return float(df.loc[df["chr"].astype(str) == "all", "ff_before_mq"].iloc[0])
    return float(df["ff_before_mq"].mean())


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--sample", required=True, type=str, help="Sample id.")
@click.option("--episcore", "episcore_tsv", required=True, type=click.Path(exists=True),
              help="{sample}_episcore.tsv with columns chr,episcore.")
@click.option("--zscore", "zscore_tsv", required=True, type=click.Path(exists=True),
              help="{sample}_zscore.tsv with columns chr,zscore.")
@click.option("--ff", "ff_tsv", default=None, type=click.Path(exists=True),
              help="SNP FF table ({sample}_ff.tsv) with an 'all' ff_before_mq row.")
@click.option("--skip-ff", is_flag=True, default=False,
              help="Omit ff_before_mq from the output.")
@click.option("--ezscore-matrix", default=None, type=click.Path(exists=True),
              help="best_ezscore_ref_20_matrix.tsv with columns chr,mean,std.")
@click.option("--skip-ezscore", is_flag=True, default=False,
              help="Omit ezscore; output episcore and zscore only.")
@click.option("--output-prefix", required=True, type=str,
              help="Output prefix; writes {prefix}_scores.tsv.")
def main(
    sample: str,
    episcore_tsv: str,
    zscore_tsv: str,
    ff_tsv: str | None,
    skip_ff: bool,
    ezscore_matrix: str | None,
    skip_ezscore: bool,
    output_prefix: str,
) -> None:
    """Merge per-chr episcore/zscore, optionally derive ezscore, attach ff_before_mq."""
    try:
        if skip_ezscore and ezscore_matrix:
            raise ValueError("Pass either --skip-ezscore or --ezscore-matrix, not both.")
        if not skip_ezscore and not ezscore_matrix:
            raise ValueError("Provide --ezscore-matrix or pass --skip-ezscore.")
        if skip_ff and ff_tsv:
            raise ValueError("Pass either --skip-ff or --ff, not both.")
        if not skip_ff and not ff_tsv:
            raise ValueError("Provide --ff or pass --skip-ff.")

        episcore = pd.read_csv(episcore_tsv, sep="\t")[["chr", "episcore"]]
        zscore = pd.read_csv(zscore_tsv, sep="\t")[["chr", "zscore"]]

        df = pd.DataFrame({"chr": CHR_LIST})
        df = df.merge(episcore, on="chr", how="left")
        df = df.merge(zscore, on="chr", how="left")

        if skip_ezscore:
            score_cols = ["chr", "episcore", "zscore"]
        else:
            ez = pd.read_csv(ezscore_matrix, sep="\t")[["chr", "mean", "std"]]
            ez = ez.rename(columns={"mean": "ez_mean", "std": "ez_std"})
            df = df.merge(ez, on="chr", how="left")
            combined = df["episcore"].to_numpy() + df["zscore"].to_numpy()
            std_safe = np.where(df["ez_std"].to_numpy() > 0, df["ez_std"].to_numpy(), np.nan)
            with np.errstate(divide="ignore", invalid="ignore"):
                ezscore = (combined - df["ez_mean"].to_numpy()) / std_safe
            df["ezscore"] = ezscore
            score_cols = ["chr", "episcore", "zscore", "ezscore"]

        df.insert(0, "sample", sample)
        if skip_ff:
            out_cols = ["sample"] + score_cols
            msg = "ff skipped"
        else:
            ff_before_mq = read_ff_before_mq(Path(ff_tsv))
            df.insert(1, "ff_before_mq", ff_before_mq)
            out_cols = ["sample", "ff_before_mq"] + score_cols
            msg = f"ff_before_mq={ff_before_mq:.4f}"

        out = df[out_cols]
        out_path = f"{output_prefix}_scores.tsv"
        out.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
        console.print(f"[green]OK[/green] Wrote {out_path} ({msg})")

    except Exception as exc:  # noqa: BLE001 - top-level reporting only
        console.print(f"\n[bold red]Error:[/bold red] {exc}", style="bold red")
        console.print_exception()
        sys.exit(1)


if __name__ == "__main__":
    main()
