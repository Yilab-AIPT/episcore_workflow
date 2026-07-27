#!/usr/bin/env python3
"""Reference-free abnormality signal from per-scheme merge_scores.

Adapts ``scripts/ref_free/ref_free_ezscore.py`` (fixed-combo path) to the
cpg-blacklist check outputs: each scheme already has per-chr episcore/zscore
in ``pipeline/<scheme>/merge_scores/``. Those matrices are re-referenced each
repeat against a random half of the dev-Normal pool; ezscore is then
z-normalised from ``episcore + zscore`` against the other half.

Because this cohort has 79 analyze/dev/Normal samples (not 80), the default
``--ref-n 40`` is auto-clamped to ``pool // 2`` (= 39) unless overridden.

Pool: all ``set=dev`` / ``Normal`` samples with scores (any ``ref_type``).
Eval: ``dev`` trisomies + all ``test`` samples with ``ff_before_mq > --min-ff``
(pool excluded).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

CHR_LIST = [f"chr{i}" for i in range(1, 23)]
SCHEME_ORDER = [
    "baseline",
    "15-site1-J",
    "15-site1-M",
    "20-site1-M",
    "15-site2-J",
    "15-site2-M",
    "20-site2-M",
]


def _ez_count_col(cutoff: float) -> str:
    return f"ezscore_abnormal_count_{cutoff:g}"


def _ez_cutoff_grid(lo: float, hi: float, step: float) -> List[float]:
    if step <= 0 or hi < lo:
        raise click.ClickException("Invalid ez cutoff grid")
    n = int(round((hi - lo) / step)) + 1
    return [round(lo + i * step, 10) for i in range(n)]


def _generate_half_partitions(
    pool_size: int,
    half: int,
    n_repeats: int,
    rng: np.random.Generator,
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    if pool_size < 2 * half:
        raise ValueError(f"pool_size={pool_size} < 2 * half={half}")
    ref_draws: List[np.ndarray] = []
    ez_draws: List[np.ndarray] = []
    for _ in range(n_repeats):
        # Draw a subset of size 2*half when pool is larger, then split.
        if pool_size == 2 * half:
            perm = rng.permutation(pool_size)
        else:
            perm = rng.choice(pool_size, size=2 * half, replace=False)
        ref_draws.append(perm[:half].astype(np.int64, copy=False))
        ez_draws.append(perm[half:].astype(np.int64, copy=False))
    return ref_draws, ez_draws


def _ref_mean_std(values: np.ndarray, ref_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """values: [n_chr, n_sample] → mean/std over ref_idx → [n_chr]."""
    ref = values[:, ref_idx]
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(ref, axis=1)
        std = np.nanstd(ref, axis=1, ddof=0)
    mean = np.where(np.isfinite(mean), mean, 0.0)
    std = np.where(np.isfinite(std), std, 0.0)
    return mean, std


def _rereference(values: np.ndarray, ref_idx: np.ndarray) -> np.ndarray:
    """Z-score each chromosome vs ref samples (same idea as compute_zscore)."""
    mean, std = _ref_mean_std(values, ref_idx)
    std_safe = np.where(std > 0, std, np.nan)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        return (values - mean[:, None]) / std_safe


def _compute_ezscore(
    episcore: np.ndarray,
    zscore: np.ndarray,
    ez_ref_idx: np.ndarray,
) -> np.ndarray:
    combined = episcore + zscore
    n_chr, _ = combined.shape
    ez = np.empty_like(combined)
    for hi in range(n_chr):
        ref_vals = combined[hi, ez_ref_idx]
        with np.errstate(invalid="ignore"):
            mu = np.nanmean(ref_vals)
            sd = np.nanstd(ref_vals, ddof=0)
        mu = mu if np.isfinite(mu) else 0.0
        sd_safe = sd if sd > 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            ez[hi] = (combined[hi] - mu) / sd_safe
    return ez


def _flag_abnormal(scores: np.ndarray, eval_idx: np.ndarray, cutoff: float) -> np.ndarray:
    sub = scores[:, eval_idx]
    with np.errstate(invalid="ignore"):
        return (np.nanmax(sub, axis=0) > cutoff).astype(np.int64)


def _flag_abnormal_multi(
    scores: np.ndarray,
    eval_idx: np.ndarray,
    cutoffs: Sequence[float],
) -> np.ndarray:
    sub = scores[:, eval_idx]
    with np.errstate(invalid="ignore"):
        max_chr = np.nanmax(sub, axis=0)
    return np.stack([(max_chr > c).astype(np.int64) for c in cutoffs], axis=0)


def load_scheme_matrices(
    score_dir: Path,
    samples: List[str],
) -> Tuple[np.ndarray, np.ndarray]:
    """Return episcore/zscore arrays shaped [n_chr, n_sample]."""
    sample_index = {s: i for i, s in enumerate(samples)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}
    ep = np.full((len(CHR_LIST), len(samples)), np.nan)
    z = np.full((len(CHR_LIST), len(samples)), np.nan)
    missing = []
    for sample in samples:
        path = score_dir / f"{sample}_scores.tsv"
        if not path.is_file():
            missing.append(sample)
            continue
        df = pd.read_csv(path, sep="\t")
        j = sample_index[sample]
        for _, row in df.iterrows():
            chrom = str(row["chr"])
            if chrom not in chr_index:
                continue
            i = chr_index[chrom]
            ep[i, j] = float(row["episcore"])
            z[i, j] = float(row["zscore"])
    if missing:
        raise click.ClickException(
            f"Missing merge_scores for {len(missing)} samples "
            f"(e.g. {missing[:5]}) under {score_dir}"
        )
    return ep, z


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--meta-csv", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--pipeline-root",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help=".../pipeline containing <scheme>/merge_scores/",
)
@click.option("--scheme", required=True, type=str)
@click.option("--output-base", required=True, type=click.Path(file_okay=False))
@click.option("--total-repeats", default=10000, show_default=True, type=int)
@click.option("--repeat-start", default=0, show_default=True, type=int)
@click.option("--repeat-end", default=None, type=int)
@click.option("--ref-n", default=40, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--cutoff", default=3.0, show_default=True, type=float)
@click.option("--ez-cutoff-min", default=3.0, show_default=True, type=float)
@click.option("--ez-cutoff-max", default=4.5, show_default=True, type=float)
@click.option("--ez-cutoff-step", default=0.1, show_default=True, type=float)
@click.option("--min-ff", default=0.01, show_default=True, type=float)
def main(
    meta_csv: str,
    pipeline_root: str,
    scheme: str,
    output_base: str,
    total_repeats: int,
    repeat_start: int,
    repeat_end: Optional[int],
    ref_n: int,
    seed: int,
    cutoff: float,
    ez_cutoff_min: float,
    ez_cutoff_max: float,
    ez_cutoff_step: float,
    min_ff: float,
) -> None:
    """Run ref-free abnormality sweep for one blacklist / baseline scheme."""
    ez_cutoffs = _ez_cutoff_grid(ez_cutoff_min, ez_cutoff_max, ez_cutoff_step)
    if repeat_end is None:
        repeat_end = total_repeats
    if repeat_start < 0 or repeat_end > total_repeats or repeat_end <= repeat_start:
        raise click.ClickException(
            f"Repeat slice [{repeat_start}, {repeat_end}) must lie within [0, {total_repeats})"
        )

    score_dir = Path(pipeline_root) / scheme / "merge_scores"
    if not score_dir.is_dir():
        raise click.ClickException(f"Missing {score_dir}")

    out_root = Path(output_base) / "ref_free_ezscore"
    out_root.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(meta_csv)
    meta = meta.drop_duplicates("sample", keep="first").copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")
    meta["is_trisomy"] = meta["label"].astype(str).str.startswith("T")

    scored = {p.name.replace("_scores.tsv", "") for p in score_dir.glob("*_scores.tsv")}
    meta = meta.loc[meta["sample"].isin(scored)].copy()

    # Pool: all set=dev / Normal with scores (any ref_type, e.g. analyze + early_ref)
    pool_mask = (meta["set"].astype(str) == "dev") & (meta["label"].astype(str) == "Normal")
    pool_samples = sorted(meta.loc[pool_mask, "sample"].tolist())
    if len(pool_samples) < 2:
        raise click.ClickException(f"Dev Normal pool too small: {len(pool_samples)}")

    # Clamp ref_n so 2*ref_n fits in the pool
    max_half = len(pool_samples) // 2
    if ref_n > max_half:
        console.print(
            f"[yellow]WARN[/yellow] ref_n={ref_n} but pool={len(pool_samples)}; "
            f"clamping to {max_half}"
        )
        ref_n = max_half
    if ref_n < 1:
        raise click.ClickException("ref_n became < 1 after clamping")

    # Eval: (dev trisomy | all test) & ff filter; exclude pool
    pool_set = set(pool_samples)
    is_dev_trisomy = (meta["set"].astype(str) == "dev") & meta["is_trisomy"]
    is_test = meta["set"].astype(str) == "test"
    eval_mask = (
        (is_dev_trisomy | is_test)
        & (meta["ff_before_mq"] > min_ff)
        & (~meta["sample"].isin(pool_set))
    )
    eval_samples = sorted(meta.loc[eval_mask, "sample"].tolist())
    if not eval_samples:
        raise click.ClickException("No eval samples after set/ff/pool filters")

    # Universe = pool ∪ eval (stable order: pool then eval)
    universe = pool_samples + [s for s in eval_samples if s not in pool_set]
    sample_index = {s: i for i, s in enumerate(universe)}
    pool_idx = np.asarray([sample_index[s] for s in pool_samples], dtype=np.int64)
    eval_idx = np.asarray([sample_index[s] for s in eval_samples], dtype=np.int64)

    console.rule(f"[bold blue]ref_free scheme={scheme}")
    console.print(f"  score dir   : {score_dir}")
    console.print(f"  output      : {out_root}")
    console.print(f"  repeats     : [{repeat_start}, {repeat_end}) of {total_repeats}")
    console.print(f"  pool        : {len(pool_samples)} dev Normal → {ref_n}+{ref_n}")
    console.print(
        f"  eval        : {len(eval_samples)} "
        f"(dev trisomy | test, ff>{min_ff})"
    )
    console.print(f"  ep/z cutoff : {cutoff}")
    console.print(
        f"  ez cutoffs  : {ez_cutoffs[0]:g} .. {ez_cutoffs[-1]:g} (n={len(ez_cutoffs)})"
    )

    ep_mat, z_mat = load_scheme_matrices(score_dir, universe)

    if repeat_start == 0:
        meta_idx = meta.set_index("sample")
        eval_info = pd.DataFrame(
            {
                "sample": eval_samples,
                "set": [meta_idx.loc[s, "set"] for s in eval_samples],
                "label": [meta_idx.loc[s, "label"] for s in eval_samples],
                "ff_before_mq": [meta_idx.loc[s, "ff_before_mq"] for s in eval_samples],
            }
        )
        if "ref_type" in meta_idx.columns:
            eval_info["ref_type"] = [meta_idx.loc[s, "ref_type"] for s in eval_samples]
        eval_info.to_csv(out_root / "eval_samples.tsv", sep="\t", index=False)
        run_config = {
            "scheme": scheme,
            "ref_n": ref_n,
            "ez_ref_n": ref_n,
            "normal_pool_size": len(pool_samples),
            "pool_note": (
                f"set=dev & Normal (any ref_type) n={len(pool_samples)}; "
                f"using {ref_n}+{ref_n} per repeat"
            ),
            "cutoff": cutoff,
            "ez_cutoffs": ez_cutoffs,
            "ez_cutoff_min": ez_cutoff_min,
            "ez_cutoff_max": ez_cutoff_max,
            "ez_cutoff_step": ez_cutoff_step,
            "total_repeats": total_repeats,
            "seed": seed,
            "min_ff": min_ff,
            "eval_filter": "(set==dev & trisomy) | set==test; ff_before_mq>min_ff; not in pool",
            "n_ep_combos": 1,
            "n_z_combos": 1,
            "n_ez_combos": 1,
            "n_eval_samples": len(eval_samples),
            "score_source": "merge_scores re-referenced (fixed combo already applied)",
        }
        (out_root / "run_config.json").write_text(json.dumps(run_config, indent=2) + "\n")
        (out_root / "pool_samples.txt").write_text("\n".join(pool_samples) + "\n")

    rng = np.random.default_rng(seed)
    ref_local, ez_local = _generate_half_partitions(
        pool_size=pool_idx.size, half=ref_n, n_repeats=total_repeats, rng=rng
    )

    n_eval = eval_idx.size
    ep_counts = np.zeros(n_eval, dtype=np.int64)
    z_counts = np.zeros(n_eval, dtype=np.int64)
    ez_counts = np.zeros((len(ez_cutoffs), n_eval), dtype=np.int64)

    for repeat_index in range(repeat_start, repeat_end):
        ref_idx = pool_idx[ref_local[repeat_index]]
        ez_ref_idx = pool_idx[ez_local[repeat_index]]

        episcore = _rereference(ep_mat, ref_idx)
        zscore = _rereference(z_mat, ref_idx)
        ep_counts += _flag_abnormal(episcore, eval_idx, cutoff)
        z_counts += _flag_abnormal(zscore, eval_idx, cutoff)
        ez = _compute_ezscore(episcore, zscore, ez_ref_idx)
        ez_counts += _flag_abnormal_multi(ez, eval_idx, ez_cutoffs)

        done = repeat_index - repeat_start + 1
        if done % 50 == 0 or done == repeat_end - repeat_start:
            console.print(f"  completed repeat {repeat_index + 1}/{repeat_end}")

    slice_data: Dict[str, np.ndarray] = {
        "eval_pos": np.arange(n_eval, dtype=np.int64),
        "episcore_abnormal_count": ep_counts,
        "zscore_abnormal_count": z_counts,
    }
    for i, c in enumerate(ez_cutoffs):
        slice_data[_ez_count_col(c)] = ez_counts[i]
    slice_path = out_root / f"abnormality_counts_{repeat_start}_{repeat_end}.tsv"
    pd.DataFrame(slice_data).to_csv(slice_path, sep="\t", index=False)
    console.print(f"[green]Done[/green] {repeat_end - repeat_start} repeats → {slice_path}")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
