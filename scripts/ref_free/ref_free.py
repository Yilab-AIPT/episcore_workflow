#!/usr/bin/env python3
"""
Reference-free episcore / zscore abnormality signal sweep.

For each repeat:
    1. Draw ``ref-n`` (default 50) reference samples from dev-set Normal samples.
    2. For **every** (threshold, recall) combo in the episcore and zscore grid-search
       parquets, compute reference-normalized episcore / zscore.
    3. Flag eval samples as abnormal when any chromosome score exceeds ``--cutoff``.

Evaluation set (fixed across repeats):
    dev trisomy samples + all test-set samples (from meta.csv ``set`` / ``label``).

Per-sample signal ratio denominators (``--combo-mode all``):
    episcore : n_ep_combos * total_repeats
    zscore   : n_z_combos * total_repeats
    either   : n_union_combos * total_repeats
               (union of episcore + zscore combo keys; abnormal when ep or z flags)

With ``--combo-mode fixed``, one episcore and one zscore combo are used; all
denominators reduce to ``total_repeats``.

Inputs (under ``--input-dir``):
    episcore_grid_search.parquet
    zscore_grid_search.parquet
    meta.csv

Outputs (under ``--output-base``):
    ref_free/eval_samples.tsv
    ref_free/run_config.json
    ref_free/abnormality_counts_{start}_{end}.tsv
    ref_free/manifest_{start}_{end}.tsv
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
from rich.console import Console

SCRIPT_DIR = Path(__file__).resolve().parent
REF40_DIR = SCRIPT_DIR.parent / "ref_explore_plus_grid_search"
if str(REF40_DIR) not in sys.path:
    sys.path.insert(0, str(REF40_DIR))

from grid_search_ref40 import (  # noqa: E402
    CHR_LIST,
    _build_dense,
    compute_episcore,
    compute_zscore,
    generate_unique_ref_draws,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)
console = Console()

DEFAULT_CUTOFF = 3.0
Combo = Tuple[float, float]


def _combo_list(df: pd.DataFrame) -> List[Combo]:
    thr = df["threshold"].astype(float).to_numpy()
    rec = df["recall"].astype(float).to_numpy()
    return sorted({(float(t), float(r)) for t, r in zip(thr, rec)})


def _combo_index(combos: List[Combo]) -> Dict[Combo, int]:
    return {c: i for i, c in enumerate(combos)}


def _accumulate_combo_flags(
    scores: np.ndarray,
    eval_idx: np.ndarray,
    cutoff: float,
) -> np.ndarray:
    """Sum abnormal flags over combos for each eval sample.

    ``scores`` shape: [n_combo, n_chr, n_sample]
    Returns shape: [n_eval]
    """
    sub = scores[:, :, eval_idx]
    flags = (sub > cutoff).any(axis=1)
    return flags.sum(axis=0).astype(np.int64)


def _load_fixed_combo_arrays(
    ep_df: pd.DataFrame,
    z_df: pd.DataFrame,
    ep_threshold: float,
    ep_recall: float,
    z_threshold: float,
    z_recall: float,
    sample_index: Dict[str, int],
    chr_index: Dict[str, int],
) -> Tuple[List[np.ndarray], np.ndarray]:
    """Return episcore inputs and zscore percentage for one combo each."""
    ep_sub = ep_df[
        (ep_df["threshold"].astype(float) == ep_threshold)
        & (ep_df["recall"].astype(float) == ep_recall)
    ]
    z_sub = z_df[
        (z_df["threshold"].astype(float) == z_threshold)
        & (z_df["recall"].astype(float) == z_recall)
    ]
    if ep_sub.empty:
        raise click.ClickException(
            f"No episcore rows for threshold={ep_threshold}, recall={ep_recall}"
        )
    if z_sub.empty:
        raise click.ClickException(
            f"No zscore rows for threshold={z_threshold}, recall={z_recall}"
        )
    _, ep_arrays = _build_dense(
        ep_sub,
        ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
        sample_index,
        chr_index,
    )
    _, z_arrays = _build_dense(z_sub, ["percentage"], sample_index, chr_index)
    return [arr[0] for arr in ep_arrays], z_arrays[0][0]


def _flag_abnormal(scores: np.ndarray, eval_idx: np.ndarray, cutoff: float) -> np.ndarray:
    """Return bool vector over eval samples: True when any chr score > cutoff."""
    sub = scores[:, eval_idx]
    return (sub > cutoff).any(axis=0)


def _accumulate_either_union(
    episcore_all: np.ndarray,
    zscore_all: np.ndarray,
    eval_idx: np.ndarray,
    cutoff: float,
    ep_combos: List[Combo],
    z_combos: List[Combo],
    union_combos: List[Combo],
) -> np.ndarray:
    """Sum either-abnormal flags over union combos for each eval sample."""
    ep_index = _combo_index(ep_combos)
    z_index = _combo_index(z_combos)
    n_eval = eval_idx.size
    counts = np.zeros(n_eval, dtype=np.int64)

    ep_sub = episcore_all[:, :, eval_idx]
    z_sub = zscore_all[:, :, eval_idx]
    ep_flags = (ep_sub > cutoff).any(axis=1)
    z_flags = (z_sub > cutoff).any(axis=1)

    for combo in union_combos:
        either = np.zeros(n_eval, dtype=bool)
        ep_i = ep_index.get(combo)
        z_i = z_index.get(combo)
        if ep_i is not None:
            either |= ep_flags[ep_i]
        if z_i is not None:
            either |= z_flags[z_i]
        counts += either.astype(np.int64)
    return counts


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Dir with episcore_grid_search.parquet, zscore_grid_search.parquet, meta.csv",
)
@click.option(
    "--output-base",
    required=True,
    type=click.Path(file_okay=False),
    help="Output base; repeat dirs land under ref_free/",
)
@click.option("--total-repeats", default=10000, show_default=True, type=int)
@click.option("--repeat-start", default=0, show_default=True, type=int)
@click.option("--repeat-end", default=None, type=int, help="Exclusive. Default: total-repeats")
@click.option(
    "--ref-n",
    default=50,
    show_default=True,
    type=int,
    help="Reference samples drawn per repeat from dev Normal pool",
)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option(
    "--cutoff",
    default=DEFAULT_CUTOFF,
    show_default=True,
    type=float,
    help="Abnormality threshold on per-chr episcore/zscore",
)
@click.option(
    "--min-ff",
    default=0.0,
    show_default=True,
    type=float,
    help="Require ff_before_mq > min_ff for samples in the universe",
)
@click.option(
    "--combo-mode",
    default="all",
    show_default=True,
    type=click.Choice(["all", "fixed"]),
    help="all: every parquet combo; fixed: one episcore + one zscore combo",
)
@click.option("--ep-threshold", default=None, type=float, help="Fixed episcore threshold")
@click.option("--ep-recall", default=None, type=float, help="Fixed episcore recall")
@click.option("--z-threshold", default=None, type=float, help="Fixed zscore threshold")
@click.option("--z-recall", default=None, type=float, help="Fixed zscore recall")
def main(
    input_dir: str,
    output_base: str,
    total_repeats: int,
    repeat_start: int,
    repeat_end: Optional[int],
    ref_n: int,
    seed: int,
    cutoff: float,
    min_ff: float,
    combo_mode: str,
    ep_threshold: Optional[float],
    ep_recall: Optional[float],
    z_threshold: Optional[float],
    z_recall: Optional[float],
) -> None:
    """Run reference-free random-ref episcore/zscore abnormality sweep."""
    input_path = Path(input_dir)
    out_root = Path(output_base) / "ref_free"
    out_root.mkdir(parents=True, exist_ok=True)

    if repeat_end is None:
        repeat_end = total_repeats
    if repeat_start < 0 or repeat_end > total_repeats or repeat_end <= repeat_start:
        raise click.ClickException(
            f"Repeat slice [{repeat_start}, {repeat_end}) must lie within [0, {total_repeats})"
        )

    console.rule("[bold blue]Reference-free episcore/zscore sweep")
    console.print(f"  Input dir      : {input_path}")
    console.print(f"  Output root    : {out_root}")
    console.print(f"  Repeat range   : [{repeat_start}, {repeat_end}) of {total_repeats}")
    use_fixed = combo_mode == "fixed"
    if use_fixed:
        missing = [
            name
            for name, val in (
                ("ep-threshold", ep_threshold),
                ("ep-recall", ep_recall),
                ("z-threshold", z_threshold),
                ("z-recall", z_recall),
            )
            if val is None
        ]
        if missing:
            raise click.ClickException(
                f"--combo-mode fixed requires: {', '.join(missing)}"
            )
    assert ep_threshold is not None or not use_fixed
    assert ep_recall is not None or not use_fixed
    assert z_threshold is not None or not use_fixed
    assert z_recall is not None or not use_fixed

    console.print(f"  ref-n / seed   : {ref_n} / {seed}")
    console.print(f"  combo-mode     : {combo_mode}")
    if use_fixed:
        console.print(f"  episcore combo : threshold={ep_threshold}, recall={ep_recall}")
        console.print(f"  zscore combo   : threshold={z_threshold}, recall={z_recall}")
    console.print(f"  cutoff         : {cutoff}")

    meta = pd.read_csv(input_path / "meta.csv")
    for col in ("sample", "set", "label", "ff_before_mq"):
        if col not in meta.columns:
            raise click.ClickException(f"meta.csv missing column: {col}")
    meta = meta.drop_duplicates("sample", keep="first").copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")

    console.print("[cyan]Loading parquets ...[/cyan]")
    ep_df = pd.read_parquet(input_path / "episcore_grid_search.parquet")
    z_df = pd.read_parquet(input_path / "zscore_grid_search.parquet")

    ep_samples = set(ep_df["sample"].astype(str).unique())
    z_samples = set(z_df["sample"].astype(str).unique())
    meta_samples = set(meta["sample"])
    ff_pass = set(meta.loc[meta["ff_before_mq"] > min_ff, "sample"].astype(str))
    universe = sorted(meta_samples & ep_samples & z_samples & ff_pass)
    if not universe:
        raise click.ClickException("No samples remain after filters")

    sample_index = {s: i for i, s in enumerate(universe)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}

    meta_idx = meta.set_index("sample").reindex(universe)
    set_arr = meta_idx["set"].astype(str).to_numpy()
    label_arr = meta_idx["label"].astype(str).to_numpy()
    ff_arr = pd.to_numeric(meta_idx["ff_before_mq"], errors="coerce").to_numpy()

    is_dev_normal = (set_arr == "dev") & (label_arr == "Normal")
    is_dev_trisomy = (set_arr == "dev") & np.char.startswith(label_arr.astype(str), "T")
    is_test = set_arr == "test"
    eval_mask = is_dev_trisomy | is_test
    ref_pool_idx = np.flatnonzero(is_dev_normal)
    eval_idx = np.flatnonzero(eval_mask)

    if ref_pool_idx.size < ref_n:
        raise click.ClickException(
            f"Need at least {ref_n} dev Normal samples, found {ref_pool_idx.size}"
        )
    if eval_idx.size == 0:
        raise click.ClickException("No evaluation samples (dev trisomy + test)")

    if use_fixed:
        assert ep_threshold is not None and ep_recall is not None
        assert z_threshold is not None and z_recall is not None
        ep_arrays, z_array = _load_fixed_combo_arrays(
            ep_df,
            z_df,
            ep_threshold,
            ep_recall,
            z_threshold,
            z_recall,
            sample_index,
            chr_index,
        )
        ep_combos: List[Combo] = [(ep_threshold, ep_recall)]
        z_combos = [(z_threshold, z_recall)]
        union_combos = ep_combos
        z_array_all = None
    else:
        ep_combos, ep_arrays = _build_dense(
            ep_df,
            ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
            sample_index,
            chr_index,
        )
        z_combos, z_arrays = _build_dense(z_df, ["percentage"], sample_index, chr_index)
        z_array_all = z_arrays[0]
        union_combos = sorted(set(ep_combos) | set(z_combos))

    console.print(f"  universe samples : {len(universe)}")
    console.print(f"  dev Normal pool  : {ref_pool_idx.size}")
    console.print(
        f"  eval samples     : {eval_idx.size} "
        f"(dev trisomy={int(is_dev_trisomy.sum())}, test={int(is_test.sum())})"
    )
    console.print(f"  episcore combos  : {len(ep_combos)}")
    console.print(f"  zscore combos    : {len(z_combos)}")
    console.print(f"  union combos     : {len(union_combos)}")

    if repeat_start == 0:
        eval_info = pd.DataFrame(
            {
                "sample": [universe[i] for i in eval_idx],
                "set": set_arr[eval_idx],
                "label": label_arr[eval_idx],
                "ff_before_mq": ff_arr[eval_idx],
            }
        )
        eval_info.to_csv(out_root / "eval_samples.tsv", sep="\t", index=False)
        run_config = {
            "combo_mode": combo_mode,
            "ref_n": ref_n,
            "cutoff": cutoff,
            "total_repeats": total_repeats,
            "seed": seed,
            "n_ep_combos": len(ep_combos),
            "n_z_combos": len(z_combos),
            "n_union_combos": len(union_combos),
            "episcore_denominator": len(ep_combos) * total_repeats,
            "zscore_denominator": len(z_combos) * total_repeats,
            "either_denominator": len(union_combos) * total_repeats,
        }
        if use_fixed:
            run_config.update(
                {
                    "ep_threshold": ep_threshold,
                    "ep_recall": ep_recall,
                    "z_threshold": z_threshold,
                    "z_recall": z_recall,
                }
            )
        (out_root / "run_config.json").write_text(json.dumps(run_config, indent=2) + "\n")
        console.print(f"[green]OK[/green] Wrote {out_root / 'eval_samples.tsv'}")

    rng = np.random.default_rng(seed)
    ref_draws = generate_unique_ref_draws(
        pool_size=ref_pool_idx.size,
        ref_n=ref_n,
        n_repeats=total_repeats,
        rng=rng,
    )

    n_eval = eval_idx.size
    ep_counts = np.zeros(n_eval, dtype=np.int64)
    z_counts = np.zeros(n_eval, dtype=np.int64)
    either_counts = np.zeros(n_eval, dtype=np.int64)
    manifest_rows: List[Dict[str, object]] = []

    for repeat_index in range(repeat_start, repeat_end):
        local_ref = ref_pool_idx[ref_draws[repeat_index]]
        if use_fixed:
            episcore = compute_episcore(
                np.expand_dims(ep_arrays[0], 0),
                np.expand_dims(ep_arrays[1], 0),
                np.expand_dims(ep_arrays[2], 0),
                np.expand_dims(ep_arrays[3], 0),
                local_ref,
            )[0]
            zscore = compute_zscore(np.expand_dims(z_array, 0), local_ref)[0]
            ep_step = _flag_abnormal(episcore, eval_idx, cutoff).astype(np.int64)
            z_step = _flag_abnormal(zscore, eval_idx, cutoff).astype(np.int64)
            either_step = (ep_step.astype(bool) | z_step.astype(bool)).astype(np.int64)
        else:
            assert z_array_all is not None
            episcore_all = compute_episcore(
                ep_arrays[0], ep_arrays[1], ep_arrays[2], ep_arrays[3], local_ref
            )
            zscore_all = compute_zscore(z_array_all, local_ref)
            ep_step = _accumulate_combo_flags(episcore_all, eval_idx, cutoff)
            z_step = _accumulate_combo_flags(zscore_all, eval_idx, cutoff)
            either_step = _accumulate_either_union(
                episcore_all,
                zscore_all,
                eval_idx,
                cutoff,
                ep_combos,
                z_combos,
                union_combos,
            )
        ep_counts += ep_step
        z_counts += z_step
        either_counts += either_step

        manifest_rows.append(
            {
                "repeat_index": repeat_index,
                "n_ref": int(local_ref.size),
                "n_eval_episcore_abnormal_total": int(ep_step.sum()),
                "n_eval_zscore_abnormal_total": int(z_step.sum()),
                "n_eval_either_abnormal_total": int(either_step.sum()),
            }
        )
        done = repeat_index - repeat_start + 1
        if done % 20 == 0 or done == repeat_end - repeat_start:
            console.print(f"  completed repeat {repeat_index + 1}/{repeat_end}")

    slice_df = pd.DataFrame(
        {
            "eval_pos": np.arange(n_eval, dtype=np.int64),
            "episcore_abnormal_count": ep_counts,
            "zscore_abnormal_count": z_counts,
            "either_abnormal_count": either_counts,
        }
    )
    slice_path = out_root / f"abnormality_counts_{repeat_start}_{repeat_end}.tsv"
    slice_df.to_csv(slice_path, sep="\t", index=False)

    manifest_path = out_root / f"manifest_{repeat_start}_{repeat_end}.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    console.print(f"[green]Done[/green] {repeat_end - repeat_start} repeats")
    console.print(f"  -> {slice_path}")
    console.print(f"  -> {manifest_path}")


if __name__ == "__main__":
    main()
