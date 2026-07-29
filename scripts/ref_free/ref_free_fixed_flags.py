#!/usr/bin/env python3
"""Fixed-combo ref_free that stores per-repeat abnormality flags (dev/test only).

Fast AUC-vs-repeats design
--------------------------
Run 1× with ``total_repeats`` (e.g. 1e6), writing compact uint8 flag shards
``flags_{start}_{end}.npz`` with arrays shape ``(n_repeats_shard, n_eval)`` for
episcore / zscore / ezscore@ez-cutoff. Downstream bootstrap over those flags
gives the AUC learning curve + variation band with **no extra ref draws**.
"""

from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
from rich.console import Console

SCRIPT_DIR = Path(__file__).resolve().parent
REF40_DIR = SCRIPT_DIR.parent / "ref_explore_plus_grid_search"
if str(REF40_DIR) not in sys.path:
    sys.path.insert(0, str(REF40_DIR))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from grid_coverage import assert_dense_coverage, assert_table_coverage  # noqa: E402
from grid_search_ref40 import (  # noqa: E402
    CHR_LIST,
    _build_dense,
    compute_episcore,
    compute_zscore,
)
from ref_free_ezscore import (  # noqa: E402
    DEFAULT_REF_N,
    _compute_ezscore,
    _flag_abnormal,
    _generate_half_partitions,
    _load_fixed_combo_arrays,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)
console = Console()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--output-base", required=True, type=click.Path(file_okay=False))
@click.option("--total-repeats", default=1_000_000, show_default=True, type=int)
@click.option("--repeat-start", default=0, show_default=True, type=int)
@click.option("--repeat-end", default=None, type=int)
@click.option("--ref-n", default=DEFAULT_REF_N, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--cutoff", default=3.0, show_default=True, type=float)
@click.option("--ez-cutoff", default=3.0, show_default=True, type=float,
              help="Single ezscore cutoff used for flag storage / AUC curve")
@click.option("--ep-threshold", default=0.5, show_default=True, type=float)
@click.option("--ep-recall", default=0.65, show_default=True, type=float)
@click.option("--z-threshold", default=0.85, show_default=True, type=float)
@click.option("--z-recall", default=0.95, show_default=True, type=float)
def main(
    input_dir: str,
    output_base: str,
    total_repeats: int,
    repeat_start: int,
    repeat_end: Optional[int],
    ref_n: int,
    seed: int,
    cutoff: float,
    ez_cutoff: float,
    ep_threshold: float,
    ep_recall: float,
    z_threshold: float,
    z_recall: float,
) -> None:
    input_path = Path(input_dir)
    out_root = Path(output_base) / "fixed_flags"
    out_root.mkdir(parents=True, exist_ok=True)
    if repeat_end is None:
        repeat_end = total_repeats
    if repeat_start < 0 or repeat_end > total_repeats or repeat_end <= repeat_start:
        raise click.ClickException(
            f"Repeat slice [{repeat_start}, {repeat_end}) invalid for total={total_repeats}"
        )

    n_shard = repeat_end - repeat_start
    console.rule("[bold blue]Fixed-combo per-repeat flags (dev/test)")
    console.print(f"  repeats [{repeat_start}, {repeat_end}) / {total_repeats}")
    console.print(f"  ez cutoff for flags: {ez_cutoff:g}")

    meta = pd.read_csv(input_path / "meta.csv").drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)
    meta["ff_before_mq"] = pd.to_numeric(meta["ff_before_mq"], errors="coerce")

    ep_df = pd.read_parquet(input_path / "episcore_grid_search.parquet")
    z_df = pd.read_parquet(input_path / "zscore_grid_search.parquet")
    ep_samples = set(ep_df["sample"].astype(str).unique())
    z_samples = set(z_df["sample"].astype(str).unique())
    universe = sorted(set(meta["sample"]) & ep_samples & z_samples)
    sample_index = {s: i for i, s in enumerate(universe)}
    chr_index = {c: i for i, c in enumerate(CHR_LIST)}

    meta_idx = meta.set_index("sample").reindex(universe)
    set_arr = meta_idx["set"].astype(str).to_numpy()
    label_arr = meta_idx["label"].astype(str).to_numpy()
    ff_arr = pd.to_numeric(meta_idx["ff_before_mq"], errors="coerce").to_numpy()

    is_trisomy = np.array([bool(re.match(r"^T\d", s)) for s in label_arr])
    is_normal = label_arr == "Normal"
    is_dev_normal = (set_arr == "dev") & is_normal
    is_dev_trisomy = (set_arr == "dev") & is_trisomy
    is_test = set_arr == "test"
    # Explicitly exclude val
    eval_mask = is_dev_trisomy | is_test
    ref_pool_idx = np.flatnonzero(is_dev_normal)
    eval_idx = np.flatnonzero(eval_mask)

    if ref_pool_idx.size < 2 * ref_n:
        raise click.ClickException(f"Need >= {2*ref_n} dev Normal, found {ref_pool_idx.size}")
    if eval_idx.size == 0:
        raise click.ClickException("No eval samples")

    ep_arrays, z_array = _load_fixed_combo_arrays(
        ep_df, z_df, ep_threshold, ep_recall, z_threshold, z_recall,
        sample_index, chr_index,
    )
    ep_dense = np.expand_dims(ep_arrays[0], 0)
    z_dense = np.expand_dims(z_array, 0)
    assert_table_coverage(ep_df, universe, "episcore", [(ep_threshold, ep_recall)])
    assert_table_coverage(z_df, universe, "zscore", [(z_threshold, z_recall)])
    assert_dense_coverage(ep_dense, universe, [(ep_threshold, ep_recall)], "episcore")
    assert_dense_coverage(z_dense, universe, [(z_threshold, z_recall)], "zscore")

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
        cfg = {
            "mode": "fixed_flags",
            "total_repeats": total_repeats,
            "ref_n": ref_n,
            "seed": seed,
            "cutoff": cutoff,
            "ez_cutoff": ez_cutoff,
            "ep_threshold": ep_threshold,
            "ep_recall": ep_recall,
            "z_threshold": z_threshold,
            "z_recall": z_recall,
            "n_eval": int(eval_idx.size),
            "n_ref_pool": int(ref_pool_idx.size),
            "eval_sets": ["dev_trisomy", "test"],
        }
        (out_root / "run_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
        console.print(f"[green]OK[/green] eval n={eval_idx.size} (dev T + test, no val)")

    rng = np.random.default_rng(seed)
    ref_draws, ez_draws = _generate_half_partitions(
        pool_size=ref_pool_idx.size, half=ref_n, n_repeats=total_repeats, rng=rng
    )

    n_eval = eval_idx.size
    flags_ep = np.zeros((n_shard, n_eval), dtype=np.uint8)
    flags_z = np.zeros((n_shard, n_eval), dtype=np.uint8)
    flags_ez = np.zeros((n_shard, n_eval), dtype=np.uint8)

    for local_i, repeat_index in enumerate(range(repeat_start, repeat_end)):
        ref_idx = ref_pool_idx[ref_draws[repeat_index]]
        ez_ref_idx = ref_pool_idx[ez_draws[repeat_index]]
        episcore = compute_episcore(
            np.expand_dims(ep_arrays[0], 0),
            np.expand_dims(ep_arrays[1], 0),
            np.expand_dims(ep_arrays[2], 0),
            np.expand_dims(ep_arrays[3], 0),
            ref_idx,
        )[0]
        zscore = compute_zscore(np.expand_dims(z_array, 0), ref_idx)[0]
        ez = _compute_ezscore(episcore, zscore, ez_ref_idx)
        flags_ep[local_i] = _flag_abnormal(episcore, eval_idx, cutoff).astype(np.uint8)
        flags_z[local_i] = _flag_abnormal(zscore, eval_idx, cutoff).astype(np.uint8)
        flags_ez[local_i] = _flag_abnormal(ez, eval_idx, ez_cutoff).astype(np.uint8)
        if (local_i + 1) % 500 == 0 or local_i + 1 == n_shard:
            console.print(f"  completed {repeat_index + 1}/{repeat_end}")

    out = out_root / f"flags_{repeat_start}_{repeat_end}.npz"
    np.savez_compressed(
        out,
        flags_ep=flags_ep,
        flags_z=flags_z,
        flags_ez=flags_ez,
        repeat_start=np.asarray([repeat_start], dtype=np.int64),
        repeat_end=np.asarray([repeat_end], dtype=np.int64),
        ez_cutoff=np.asarray([ez_cutoff], dtype=np.float64),
    )
    console.print(f"[green]Done[/green] {n_shard} repeats -> {out}")


if __name__ == "__main__":
    main()
