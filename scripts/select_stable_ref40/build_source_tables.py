#!/usr/bin/env python3
"""
Build episcore (beta) and zscore (percentage) source tables.

Episcore source (early.config combo: recall=0.65, threshold=0.5):
    --episcore-samplesheet  CSV with columns sample,episcore_file
    -> beta.csv with hypo/hyper z_intra + cpg counts (and other zscore.tsv cols)

Zscore source (external combo: recall=0.95, cutoff/threshold=0.85):
    --percentage-csv        existing aggregated percentage table (may be incomplete)
    --meta-csv              sample list to cover
    Missing samples are fetched from original pipeline directories:
      1) chr_percentage.mTcount_1.tsv under zscore_data.CpG_recall0.95
      2) fallback: *.0.85.1.0.CpG_final_filtered_recall0.95.NoLen.zscore.csv
         (percentage as fraction -> *100)

Writes under --output-dir:
    beta.csv
    percentage.csv                 (cutoff=0.85 only, percent scale)
    episcore_source_summary.tsv
    zscore_source_summary.tsv
    zscore_fetch_report.tsv        (per missing sample: found/not_found + path/reason)
    zscore_fetch_summary.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
DEFAULT_CUTOFF = 0.85
DEFAULT_PIPELINE_ROOTS = (
    "/lustre1/cqyi/myli/bert/DNA_5mC_analysis_pipeline/output",
    "/appsnew/home/myli/lustre1/bert/DNA_5mC_analysis_pipeline/output",
)

NEEDED_BETA_SUFFIXES = (
    "hypo_z_intra",
    "hyper_z_intra",
    "hypo_cpgs_count",
    "hyper_cpgs_count",
    "hypo_beta",
    "hyper_beta",
    "s_intra",
    "hypo_z_inter",
    "hyper_z_inter",
    "s_inter",
)


def _batch_dirs(roots: List[Path]) -> List[Path]:
    out: List[Path] = []
    for root in roots:
        if root.is_dir():
            out.extend([p for p in root.iterdir() if p.is_dir()])
    return out


def _zscore_downstream_from_episcore(episcore_file: Path) -> Optional[Path]:
    for i, part in enumerate(episcore_file.parts):
        if part == "zscore_downstream":
            return Path(*episcore_file.parts[: i + 1])
    return None


def _candidate_names(sample: str, episcore_file: Optional[str]) -> List[str]:
    names = [sample]
    if episcore_file:
        base = Path(episcore_file).name.replace("_zscore.tsv", "")
        if base and base not in names:
            names.append(base)
        # P<digit> <-> B<digit> swap often used for blood vs plasma naming
        for n in list(names):
            swapped = n.replace("P", "B") if "P" in n else n.replace("B", "P")
            # only swap the first P/B after the numeric id block: keep conservative
            if "P" in n[5:]:
                # e.g. JPTAY1404P7H1 -> JPTAY1404B7H1
                idx = n.rfind("P")
                if idx > 0:
                    alt = n[:idx] + "B" + n[idx + 1 :]
                    if alt not in names:
                        names.append(alt)
            if "B" in n[5:]:
                idx = n.rfind("B")
                if idx > 0:
                    alt = n[:idx] + "P" + n[idx + 1 :]
                    if alt not in names:
                        names.append(alt)
            if swapped != n and swapped not in names:
                names.append(swapped)
    return names


def _find_percentage_paths(
    sample: str,
    *,
    episcore_file: Optional[str],
    batch_dirs: List[Path],
) -> Tuple[Optional[Path], Optional[Path], str]:
    """Return (chr_percentage_path, zscore_csv_path, note)."""
    names = _candidate_names(sample, episcore_file)
    search_roots: List[Path] = []
    if episcore_file:
        zd = _zscore_downstream_from_episcore(Path(episcore_file))
        if zd is not None:
            search_roots.append(zd)

    rel_pct = [
        "futheranalysis/zscore_data.CpG_recall0.95/{name}/chr_percentage.mTcount_1.tsv",
        "zscore_data.CpG_recall0.95/{name}/chr_percentage.mTcount_1.tsv",
        "futheranalysis/zscore_data/{name}/chr_percentage.mTcount_1.tsv",
    ]
    rel_zs = [
        "futheranalysis/zscore_data.CpG_recall0.95/{name}/{name}.0.85.1.0.CpG_final_filtered_recall0.95.NoLen.zscore.csv",
        "zscore_data.CpG_recall0.95/{name}/{name}.0.85.1.0.CpG_final_filtered_recall0.95.NoLen.zscore.csv",
        "futheranalysis/zscore_data.CpG_recall0.95/{name}/{name}.early.zscore.cutoff0.85.zscore.csv",
    ]

    pct_hit: Optional[Path] = None
    zs_hit: Optional[Path] = None

    def scan(root: Path) -> None:
        nonlocal pct_hit, zs_hit
        for name in names:
            if pct_hit is None:
                for pattern in rel_pct:
                    cand = root / pattern.format(name=name)
                    if cand.is_file():
                        pct_hit = cand
                        break
            if zs_hit is None:
                # exact patterns first
                for pattern in rel_zs:
                    cand = root / pattern.format(name=name)
                    if cand.is_file():
                        zs_hit = cand
                        break
                if zs_hit is None:
                    # glob fallback inside sample dir
                    for sub in (
                        f"futheranalysis/zscore_data.CpG_recall0.95/{name}",
                        f"zscore_data.CpG_recall0.95/{name}",
                    ):
                        d = root / sub
                        if d.is_dir():
                            hits = sorted(d.glob("*.0.85*.zscore.csv")) + sorted(
                                d.glob("*cutoff0.85*.zscore.csv")
                            )
                            # prefer NoLen files that contain percentage
                            for h in hits:
                                zs_hit = h
                                if "NoLen" in h.name:
                                    break
                        if zs_hit is not None:
                            break

    for root in search_roots:
        scan(root)
        if pct_hit is not None:
            return pct_hit, zs_hit, "via_episcore_zscore_downstream"

    for batch in batch_dirs:
        zd = batch / "bwameth_results" / "zscore_downstream"
        if not zd.is_dir():
            continue
        scan(zd)
        if pct_hit is not None or zs_hit is not None:
            return pct_hit, zs_hit, "via_batch_scan"
    return None, None, "not_found"


def _rows_from_chr_percentage(path: Path, sample: str, cutoff: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "cutoff" not in df.columns:
        raise ValueError(f"No cutoff column in {path}")
    sub = df[np.isclose(df["cutoff"].astype(float), cutoff)].copy()
    if sub.empty:
        raise ValueError(f"No rows at cutoff={cutoff} in {path}")
    # unify sample column
    if "sample_group" in sub.columns:
        sub = sub.rename(columns={"sample_group": "sample"})
    sub["sample"] = sample
    keep = ["sample", "chr", "cutoff", "count", "percentage", "min_cpg"]
    for c in keep:
        if c not in sub.columns:
            sub[c] = np.nan if c != "min_cpg" else 1
    sub["source_file"] = str(path)
    sub["source_kind"] = "chr_percentage"
    return sub[keep + ["source_file", "source_kind"]]


def _rows_from_zscore_csv(path: Path, sample: str, cutoff: float) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "percentage" not in df.columns:
        raise ValueError(f"No percentage column in {path}")
    if "chr" not in df.columns:
        raise ValueError(f"No chr column in {path}")
    out = df.copy()
    out["sample"] = sample
    out["cutoff"] = cutoff
    # fraction (sum~1) -> percent (sum~100)
    pct = out["percentage"].astype(float)
    if pct.sum() <= 1.5:
        out["percentage"] = pct * 100.0
    if "readscount" in out.columns:
        out["count"] = out["readscount"].astype(float)
    else:
        out["count"] = np.nan
    out["min_cpg"] = 1
    out["source_file"] = str(path)
    out["source_kind"] = "zscore_csv"
    return out[
        ["sample", "chr", "cutoff", "count", "percentage", "min_cpg", "source_file", "source_kind"]
    ]


def build_beta_table(episcore_samplesheet: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sheet = pd.read_csv(episcore_samplesheet)
    if not {"sample", "episcore_file"}.issubset(sheet.columns):
        raise click.ClickException(
            "episcore samplesheet must have columns: sample, episcore_file"
        )
    rows = []
    summary_rows = []
    for _, r in sheet.iterrows():
        sample = str(r["sample"])
        path = Path(str(r["episcore_file"]))
        status = "ok"
        n_cols = 0
        err = ""
        try:
            if not path.is_file():
                status = "missing_file"
                raise FileNotFoundError(path)
            df = pd.read_csv(path, sep="\t")
            if len(df) != 1:
                # still take first row
                status = f"warn_nrows_{len(df)}"
            row = df.iloc[0].to_dict()
            row["sample"] = sample
            row["episcore_file"] = str(path)
            rows.append(row)
            n_cols = len(df.columns)
        except Exception as exc:  # noqa: BLE001
            status = "error" if status == "ok" else status
            err = str(exc)
            summary_rows.append(
                {
                    "sample": sample,
                    "status": status,
                    "n_cols": n_cols,
                    "episcore_file": str(path),
                    "error": err,
                }
            )
            continue
        # check required suffixes present for chr1 at least
        missing_suf = [
            s for s in ("hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count")
            if f"chr1_{s}" not in df.columns
        ]
        if missing_suf:
            status = "missing_columns"
            err = ",".join(missing_suf)
        summary_rows.append(
            {
                "sample": sample,
                "status": status,
                "n_cols": n_cols,
                "episcore_file": str(path),
                "error": err,
            }
        )
    beta = pd.DataFrame(rows)
    if not beta.empty:
        # put sample first
        cols = ["sample"] + [c for c in beta.columns if c != "sample"]
        beta = beta[cols]
    summary = pd.DataFrame(summary_rows)
    return beta, summary


def build_percentage_table(
    meta: pd.DataFrame,
    existing_pct: pd.DataFrame,
    episcore_map: Dict[str, str],
    *,
    cutoff: float,
    pipeline_roots: List[Path],
) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    if "cutoff" in existing_pct.columns:
        base = existing_pct[np.isclose(existing_pct["cutoff"].astype(float), cutoff)].copy()
    else:
        base = existing_pct.copy()
        base["cutoff"] = cutoff
    base["sample"] = base["sample"].astype(str)
    base["chr"] = base["chr"].astype(str)
    have = set(base["sample"].unique())
    want = set(meta["sample"].astype(str))
    missing = sorted(want - have)

    batches = _batch_dirs(pipeline_roots)
    fetch_rows = []
    added_frames: List[pd.DataFrame] = []

    for sample in missing:
        ep = episcore_map.get(sample)
        pct_path, zs_path, note = _find_percentage_paths(
            sample, episcore_file=ep, batch_dirs=batches
        )
        status = "not_found"
        used = ""
        err = ""
        try:
            if pct_path is not None:
                frame = _rows_from_chr_percentage(pct_path, sample, cutoff)
                # keep autosomes only
                frame = frame[frame["chr"].isin(CHR_LIST)]
                added_frames.append(frame)
                status = "fetched_chr_percentage"
                used = str(pct_path)
            elif zs_path is not None:
                frame = _rows_from_zscore_csv(zs_path, sample, cutoff)
                frame = frame[frame["chr"].isin(CHR_LIST)]
                if frame.empty or "percentage" not in frame.columns:
                    raise ValueError("zscore csv lacked usable percentage")
                # early.zscore.cutoff files may lack percentage — already guarded
                if frame["percentage"].isna().all():
                    raise ValueError("percentage all-NA in zscore csv")
                added_frames.append(frame)
                status = "fetched_zscore_csv"
                used = str(zs_path)
            else:
                status = "not_found"
        except Exception as exc:  # noqa: BLE001
            status = "error"
            err = str(exc)
            used = str(pct_path or zs_path or "")
        fetch_rows.append(
            {
                "sample": sample,
                "status": status,
                "source_file": used,
                "search_note": note,
                "episcore_file": ep or "",
                "error": err,
                "set": meta.loc[meta["sample"].astype(str) == sample, "set"].iloc[0]
                if (meta["sample"].astype(str) == sample).any()
                else "",
                "label": meta.loc[meta["sample"].astype(str) == sample, "label"].iloc[0]
                if (meta["sample"].astype(str) == sample).any()
                else "",
            }
        )

    keep_cols = ["sample", "chr", "cutoff", "count", "percentage", "min_cpg"]
    out = base[keep_cols].copy() if set(keep_cols).issubset(base.columns) else base.copy()
    for c in keep_cols:
        if c not in out.columns:
            out[c] = np.nan
    out = out[keep_cols]
    out["source_kind"] = "existing_percentage_csv"

    if added_frames:
        add = pd.concat(added_frames, ignore_index=True)
        add = add[keep_cols + ["source_kind"]]
        out = pd.concat([out, add[keep_cols + ["source_kind"]]], ignore_index=True)

    # one row per sample-chr
    out = out.drop_duplicates(subset=["sample", "chr"], keep="first")
    fetch_df = pd.DataFrame(fetch_rows)
    summary = {
        "meta_samples": int(len(want)),
        "existing_at_cutoff": int(len(have)),
        "missing_requested": int(len(missing)),
        "fetched_chr_percentage": int((fetch_df["status"] == "fetched_chr_percentage").sum())
        if not fetch_df.empty
        else 0,
        "fetched_zscore_csv": int((fetch_df["status"] == "fetched_zscore_csv").sum())
        if not fetch_df.empty
        else 0,
        "not_found": int((fetch_df["status"] == "not_found").sum()) if not fetch_df.empty else 0,
        "error": int((fetch_df["status"] == "error").sum()) if not fetch_df.empty else 0,
        "final_unique_samples": int(out["sample"].nunique()),
        "cutoff": cutoff,
    }
    return out, fetch_df, summary


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--episcore-samplesheet",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="CSV with sample,episcore_file (early combo recall=0.65 / thr=0.5)",
)
@click.option(
    "--percentage-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Existing percentage.csv (may be incomplete)",
)
@click.option(
    "--meta-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Meta samplesheet defining the sample universe",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Output directory for beta/percentage/summaries",
)
@click.option("--cutoff", default=DEFAULT_CUTOFF, show_default=True, type=float)
@click.option(
    "--pipeline-root",
    multiple=True,
    default=DEFAULT_PIPELINE_ROOTS,
    show_default=True,
    help="DNA_5mC pipeline output roots to search for missing percentage files",
)
def main(
    episcore_samplesheet: str,
    percentage_csv: str,
    meta_csv: str,
    output_dir: str,
    cutoff: float,
    pipeline_root: Tuple[str, ...],
) -> None:
    """Build beta.csv + percentage.csv and fetch reports."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    console.rule("[bold blue]Build episcore / zscore source tables")

    console.print("Building beta from episcore zscore.tsv files ...")
    beta, epi_summary = build_beta_table(Path(episcore_samplesheet))
    beta_path = out_dir / "beta.csv"
    # Drop helper path col from scoring table but keep in summary
    beta_out = beta.drop(columns=["episcore_file"], errors="ignore")
    beta_out.to_csv(beta_path, index=False)
    epi_summary.to_csv(out_dir / "episcore_source_summary.tsv", sep="\t", index=False)
    console.print(
        f"[green]OK[/green] beta.csv: {len(beta_out)} samples "
        f"(ok={(epi_summary.status == 'ok').sum()}, "
        f"issues={(epi_summary.status != 'ok').sum()})"
    )

    meta = pd.read_csv(meta_csv).drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)
    existing = pd.read_csv(percentage_csv, sep="\t")
    epi_sheet = pd.read_csv(episcore_samplesheet)
    epi_map = {
        str(s): str(p)
        for s, p in zip(epi_sheet["sample"], epi_sheet["episcore_file"])
    }

    console.print(f"Completing percentage at cutoff={cutoff} ...")
    pct, fetch_df, fetch_summary = build_percentage_table(
        meta,
        existing,
        epi_map,
        cutoff=cutoff,
        pipeline_roots=[Path(p) for p in pipeline_root],
    )
    # scoring input: no source_kind required by downstream, but keep a full copy
    pct_score = pct[["sample", "chr", "cutoff", "count", "percentage", "min_cpg"]].copy()
    pct_score.to_csv(out_dir / "percentage.csv", sep="\t", index=False)
    pct.to_csv(out_dir / "percentage_with_source.csv", sep="\t", index=False)
    fetch_df.to_csv(out_dir / "zscore_fetch_report.tsv", sep="\t", index=False)
    (out_dir / "zscore_fetch_summary.json").write_text(json.dumps(fetch_summary, indent=2))

    # compact zscore coverage summary vs meta
    covered = set(pct_score["sample"])
    zsum_rows = []
    for _, r in meta.iterrows():
        s = str(r["sample"])
        zsum_rows.append(
            {
                "sample": s,
                "set": r.get("set", ""),
                "label": r.get("label", ""),
                "ref_type": r.get("ref_type", ""),
                "has_percentage": s in covered,
                "n_chr": int((pct_score["sample"] == s).sum()) if s in covered else 0,
            }
        )
    zsum = pd.DataFrame(zsum_rows)
    zsum.to_csv(out_dir / "zscore_source_summary.tsv", sep="\t", index=False)

    console.print("[green]OK[/green] percentage.csv:", fetch_summary)
    if not fetch_df.empty:
        console.print("Fetch status counts:", fetch_df["status"].value_counts().to_dict())
        nf = fetch_df[fetch_df["status"].isin(["not_found", "error"])]
        if not nf.empty:
            console.print(
                f"[yellow]Unresolved[/yellow] {len(nf)} samples "
                f"(sets={nf['set'].value_counts().to_dict()})"
            )
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
