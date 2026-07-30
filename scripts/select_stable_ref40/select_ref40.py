#!/usr/bin/env python3
"""
Select a Normal+dev ref_40 whose per-chr mean/std match early_ref_17, while
minimizing ezscore pred_label changes (cutoff 4.5) versus the early_ref baseline.

Inputs (under --input-dir, from build_source_tables.py):
    meta.csv or --meta-csv
    beta.csv
    percentage.csv

Outputs (under --output-dir):
    ref40_samples.txt
    ref40_samples.tsv
    baseline_score.tsv / ref40_score.tsv   (all samples with scores)
    reference_meanstd_compare.tsv
    pred_label_compare.tsv
    selection_summary.json
    search_top_candidates.tsv
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
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "reference_explore"))
sys.path.insert(0, str(SCRIPT_DIR))

from calc_zscore_episcore_ezscore import (  # noqa: E402
    DEFAULT_EZSCORE_REF_SAMPLES,
    _compute_episcore,
    _load_percentage_matrix,
    _parse_chr_spec,
    _ref_mean_std,
    _stack_per_chr,
    _zscore_from_ref,
    build_ezscore_ref_mask,
    load_ezscore_ref_samples,
)
from pred_label_utils import (  # noqa: E402
    STRONG_CUTOFF,
    assign_pred_labels_matrix,
    format_comma_scores,
)

console = Console()
CHR_LIST = _parse_chr_spec("1-22")


def _load_merged(input_dir: Path, meta_csv: Path) -> Tuple[pd.DataFrame, Path]:
    beta = pd.read_csv(input_dir / "beta.csv").drop_duplicates("sample", keep="first")
    beta["sample"] = beta["sample"].astype(str)
    meta = pd.read_csv(meta_csv).drop_duplicates("sample", keep="first")
    meta["sample"] = meta["sample"].astype(str)
    merged = meta.merge(beta, on="sample", how="inner", suffixes=("", "_beta"))
    pct_path = input_dir / "percentage.csv"
    if not pct_path.is_file():
        raise FileNotFoundError(pct_path)
    return merged, pct_path


def _prepare_arrays(
    merged: pd.DataFrame, pct_path: Path
) -> Dict[str, np.ndarray]:
    hypo_z = _stack_per_chr(merged, CHR_LIST, "hypo_z_intra")
    hyper_z = _stack_per_chr(merged, CHR_LIST, "hyper_z_intra")
    hypo_c = _stack_per_chr(merged, CHR_LIST, "hypo_cpgs_count", dtype=np.int64)
    hyper_c = _stack_per_chr(merged, CHR_LIST, "hyper_cpgs_count", dtype=np.int64)
    pct = _load_percentage_matrix(pct_path, merged["sample"].tolist(), CHR_LIST)
    return {
        "hypo_z": hypo_z,
        "hyper_z": hyper_z,
        "hypo_c": hypo_c,
        "hyper_c": hyper_c,
        "pct": pct,
    }


def _mean_std_distance(
    arrays: Dict[str, np.ndarray],
    ref_idx: np.ndarray,
    target: Dict[str, Tuple[np.ndarray, np.ndarray]],
) -> float:
    """Normalized L2 distance of mean/std vs early_ref target (hypo/hyper/pct)."""
    dist = 0.0
    for key, arr in (
        ("hypo_z", arrays["hypo_z"]),
        ("hyper_z", arrays["hyper_z"]),
        ("pct", arrays["pct"]),
    ):
        mu, sd = _ref_mean_std(arr[ref_idx])
        t_mu, t_sd = target[key]
        # scale by target sd (or 1) so chromosomes are comparable
        scale_mu = np.where(np.abs(t_sd) > 1e-8, np.abs(t_sd), 1.0)
        scale_sd = np.where(np.abs(t_sd) > 1e-8, np.abs(t_sd), 1.0)
        dist += float(np.nansum(((mu - t_mu) / scale_mu) ** 2))
        dist += float(np.nansum(((sd - t_sd) / scale_sd) ** 2))
    return dist


def _compute_all_scores(
    arrays: Dict[str, np.ndarray],
    ref_idx: np.ndarray,
    ez_idx: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ref_hypo_mean, ref_hypo_std = _ref_mean_std(arrays["hypo_z"][ref_idx])
    ref_hyper_mean, ref_hyper_std = _ref_mean_std(arrays["hyper_z"][ref_idx])
    episcore = _compute_episcore(
        arrays["hypo_z"],
        arrays["hyper_z"],
        arrays["hypo_c"],
        arrays["hyper_c"],
        ref_hypo_mean,
        ref_hypo_std,
        ref_hyper_mean,
        ref_hyper_std,
    )
    pct_mean, pct_std = _ref_mean_std(arrays["pct"][ref_idx])
    zscore = _zscore_from_ref(arrays["pct"], pct_mean, pct_std)
    combined = zscore + episcore
    ez_mean, ez_std = _ref_mean_std(combined[ez_idx])
    ezscore = _zscore_from_ref(combined, ez_mean, ez_std)
    return zscore, episcore, ezscore


def _pred_masks(ez: np.ndarray, strong: float = STRONG_CUTOFF, gray: float = 3.0):
    strong_m = ez > strong
    gray_m = (ez >= gray) & (ez <= strong)
    return strong_m, gray_m


def _pred_change_count_masks(
    base_strong: np.ndarray,
    base_gray: np.ndarray,
    new_strong: np.ndarray,
    new_gray: np.ndarray,
    compare_mask: np.ndarray,
) -> int:
    changed = ((base_strong != new_strong) | (base_gray != new_gray)).any(axis=1)
    return int(np.sum(changed & compare_mask))


def _pred_change_count(
    baseline_pred: np.ndarray,
    new_pred: np.ndarray,
    compare_mask: np.ndarray,
) -> int:
    return int(np.sum((baseline_pred != new_pred) & compare_mask))


def _score_table(
    samples: Sequence[str],
    zscore: np.ndarray,
    episcore: np.ndarray,
    ezscore: np.ndarray,
) -> pd.DataFrame:
    data = {"sample": list(samples)}
    for i, chr_name in enumerate(CHR_LIST):
        num = chr_name.removeprefix("chr")
        data[f"zscore_chr{num}"] = zscore[:, i]
        data[f"episcore_chr{num}"] = episcore[:, i]
        data[f"ezscore_chr{num}"] = ezscore[:, i]
    df = pd.DataFrame(data)
    df["pred_label"] = assign_pred_labels_matrix(ezscore, strong_cutoff=STRONG_CUTOFF)
    df["rc_zscores"] = [format_comma_scores(r) for r in zscore]
    df["beta_zscores"] = [format_comma_scores(r) for r in episcore]
    df["final_zscores"] = [format_comma_scores(r) for r in ezscore]
    return df


def _meanstd_compare_table(
    arrays: Dict[str, np.ndarray],
    early_idx: np.ndarray,
    ref40_idx: np.ndarray,
) -> pd.DataFrame:
    rows = []
    for key, label in (
        ("hypo_z", "episcore_hypo_z_intra"),
        ("hyper_z", "episcore_hyper_z_intra"),
        ("pct", "zscore_percentage"),
    ):
        e_mu, e_sd = _ref_mean_std(arrays[key][early_idx])
        r_mu, r_sd = _ref_mean_std(arrays[key][ref40_idx])
        for i, chr_name in enumerate(CHR_LIST):
            rows.append(
                {
                    "feature": label,
                    "chr": chr_name,
                    "early_ref_mean": e_mu[i],
                    "early_ref_std": e_sd[i],
                    "ref40_mean": r_mu[i],
                    "ref40_std": r_sd[i],
                    "delta_mean": r_mu[i] - e_mu[i],
                    "delta_std": r_sd[i] - e_sd[i],
                }
            )
    return pd.DataFrame(rows)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option(
    "--meta-csv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Default: <input-dir>/meta.csv",
)
@click.option("--output-dir", required=True, type=click.Path(file_okay=False))
@click.option("--ref-n", default=40, show_default=True, type=int)
@click.option("--n-random", default=20000, show_default=True, type=int)
@click.option("--n-swap-rounds", default=2000, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
@click.option(
    "--ezscore-ref-samples",
    default=str(DEFAULT_EZSCORE_REF_SAMPLES),
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--exclude-sample",
    default="PTAY0586P8S1",
    show_default=True,
    help="Optional sample never used as reference (empty string to disable)",
)
@click.option(
    "--keep-early-ref-seed/--no-keep-early-ref-seed",
    default=True,
    show_default=True,
    help="Also evaluate early_ref + (ref_n-17) fillers and use as search seed",
)
def main(
    input_dir: str,
    meta_csv: Optional[str],
    output_dir: str,
    ref_n: int,
    n_random: int,
    n_swap_rounds: int,
    seed: int,
    ezscore_ref_samples: str,
    exclude_sample: str,
    keep_early_ref_seed: bool,
) -> None:
    """Search for a stable ref_40 and write comparison artifacts."""
    in_dir = Path(input_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_path = Path(meta_csv) if meta_csv else in_dir / "meta.csv"
    if not meta_path.is_file():
        raise click.ClickException(f"Missing meta csv: {meta_path}")

    console.rule("[bold blue]Select stable ref_40")
    merged, pct_path = _load_merged(in_dir, meta_path)
    # Require finite percentage for candidates / scoring completeness
    pct_samples = set(
        pd.read_csv(pct_path, sep="\t", usecols=["sample"])["sample"].astype(str)
    )
    arrays = _prepare_arrays(merged, pct_path)
    ez_mask = build_ezscore_ref_mask(merged, samples_file=Path(ezscore_ref_samples))
    ez_idx = np.flatnonzero(ez_mask)
    console.print(f"  Samples merged : {len(merged)}")
    console.print(f"  Ezscore refs   : {ez_idx.size}")

    early_mask = merged["ref_type"].astype(str) == "early_ref"
    early_idx = np.flatnonzero(early_mask.to_numpy())
    if early_idx.size == 0:
        raise click.ClickException("No early_ref samples in meta")

    pool_mask = (merged["set"].astype(str) == "dev") & (
        merged["label"].astype(str) == "Normal"
    )
    if exclude_sample:
        pool_mask &= merged["sample"].astype(str) != exclude_sample
    # must have percentage coverage
    pool_mask &= merged["sample"].astype(str).isin(pct_samples)
    pool_idx = np.flatnonzero(pool_mask.to_numpy())
    if pool_idx.size < ref_n:
        raise click.ClickException(
            f"Candidate pool {pool_idx.size} < ref_n {ref_n}"
        )
    console.print(
        f"  early_ref      : {early_idx.size} | pool Normal+dev: {pool_idx.size}"
    )

    target = {
        "hypo_z": _ref_mean_std(arrays["hypo_z"][early_idx]),
        "hyper_z": _ref_mean_std(arrays["hyper_z"][early_idx]),
        "pct": _ref_mean_std(arrays["pct"][early_idx]),
    }

    # Baseline scores with early_ref
    z0, e0, ez0 = _compute_all_scores(arrays, early_idx, ez_idx)
    base_strong, base_gray = _pred_masks(ez0)
    baseline_pred = np.array(assign_pred_labels_matrix(ez0, strong_cutoff=STRONG_CUTOFF))
    baseline_df = _score_table(merged["sample"].tolist(), z0, e0, ez0)
    baseline_df.to_csv(out_dir / "baseline_score.tsv", sep="\t", index=False)

    # Samples compared for pred_label stability: have percentage, not in early_ref
    # (early_ref historically have NA pred_label), and later also exclude new ref.
    has_pct = merged["sample"].astype(str).isin(pct_samples).to_numpy()
    base_compare = has_pct & (~early_mask.to_numpy())

    rng = np.random.default_rng(seed)
    candidates: List[dict] = []

    def evaluate(ref_global_idx: np.ndarray, tag: str, *, keep_scores: bool = False) -> dict:
        ref_global_idx = np.unique(np.asarray(ref_global_idx, dtype=np.int64))
        if ref_global_idx.size != ref_n:
            raise ValueError(f"{tag}: expected {ref_n}, got {ref_global_idx.size}")
        dist = _mean_std_distance(arrays, ref_global_idx, target)
        z, e, ez = _compute_all_scores(arrays, ref_global_idx, ez_idx)
        new_strong, new_gray = _pred_masks(ez)
        ref_mask = np.zeros(len(merged), dtype=bool)
        ref_mask[ref_global_idx] = True
        compare = base_compare & (~ref_mask)
        n_changed = _pred_change_count_masks(
            base_strong, base_gray, new_strong, new_gray, compare
        )
        n_compare = int(compare.sum())
        out = {
            "tag": tag,
            "ref_idx": ref_global_idx.copy(),
            "dist": float(dist),
            "n_changed": int(n_changed),
            "n_compare": int(n_compare),
            "reference_list": ",".join(
                merged.iloc[ref_global_idx]["sample"].astype(str).tolist()
            ),
        }
        if keep_scores:
            pred = np.array(assign_pred_labels_matrix(ez, strong_cutoff=STRONG_CUTOFF))
            out.update({"z": z, "e": e, "ez": ez, "pred": pred})
        return out

    best: Optional[dict] = None

    def consider(res: dict) -> None:
        nonlocal best
        light = {
            "tag": res["tag"],
            "dist": res["dist"],
            "n_changed": res["n_changed"],
            "n_compare": res["n_compare"],
            "reference_list": res["reference_list"],
            "ref_idx": res["ref_idx"],
        }
        candidates.append(light)
        if best is None or (res["n_changed"], res["dist"]) < (
            best["n_changed"],
            best["dist"],
        ):
            # Recompute once with scores retained for the running best
            best = evaluate(res["ref_idx"], res["tag"], keep_scores=True)

    # Seed: keep all early_ref + best fillers by mean/std distance
    early_in_pool = np.intersect1d(early_idx, pool_idx)
    fillers_pool = np.setdiff1d(pool_idx, early_in_pool)
    if keep_early_ref_seed and early_in_pool.size <= ref_n:
        need = ref_n - int(early_in_pool.size)
        console.print(f"Evaluating early_ref seed + {need} fillers ...")
        # random filler draws
        n_seed_draws = min(2000, n_random)
        for i in range(n_seed_draws):
            fill = rng.choice(fillers_pool, size=need, replace=False)
            ref = np.sort(np.concatenate([early_in_pool, fill]))
            consider(evaluate(ref, f"early_seed_{i}"))
        # also pure distance-greedy fillers
        remaining = list(map(int, fillers_pool))
        chosen = list(map(int, early_in_pool))
        while len(chosen) < ref_n and remaining:
            best_s, best_d = None, np.inf
            for s in remaining:
                trial = np.array(chosen + [s], dtype=np.int64)
                # temporary pad: if size < ref_n, distance still informative on partial set
                d = _mean_std_distance(arrays, trial, target)
                if d < best_d:
                    best_d, best_s = d, s
            chosen.append(best_s)
            remaining.remove(best_s)
        consider(evaluate(np.array(sorted(chosen), dtype=np.int64), "early_greedy"))

    # Random search over full pool
    console.print(f"Random search: {n_random} draws ...")
    seen = set()
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("random", total=n_random)
        for i in range(n_random):
            draw = np.sort(rng.choice(pool_idx, size=ref_n, replace=False))
            key = draw.tobytes()
            if key in seen:
                progress.advance(task)
                continue
            seen.add(key)
            consider(evaluate(draw, f"random_{i}"))
            progress.advance(task)

    # Local swap optimization on current best
    if best is None:
        raise click.ClickException("No candidates evaluated")
    console.print(
        f"Local swaps on best (n_changed={best['n_changed']}, dist={best['dist']:.4f}) ..."
    )
    cur_idx = best["ref_idx"].copy()
    cur_changed, cur_dist = best["n_changed"], best["dist"]
    for i in range(n_swap_rounds):
        in_ref = set(map(int, cur_idx))
        out_cand = [int(x) for x in pool_idx if int(x) not in in_ref]
        if not out_cand:
            break
        drop = int(rng.choice(cur_idx))
        add = int(rng.choice(out_cand))
        trial = np.sort(
            np.array([x for x in cur_idx if x != drop] + [add], dtype=np.int64)
        )
        res = evaluate(trial, f"swap_{i}", keep_scores=False)
        consider(res)
        if (res["n_changed"], res["dist"]) < (cur_changed, cur_dist):
            cur_idx = trial
            cur_changed, cur_dist = res["n_changed"], res["dist"]

    candidates.sort(key=lambda c: (c["n_changed"], c["dist"]))
    # Ensure best retains full score matrices
    if "z" not in best:
        best = evaluate(best["ref_idx"], best["tag"], keep_scores=True)
    # If a lighter candidate beat best on sort key, refresh
    top_light = candidates[0]
    if (top_light["n_changed"], top_light["dist"]) < (best["n_changed"], best["dist"]) or (
        top_light["ref_idx"].tobytes() != best["ref_idx"].tobytes()
        and (top_light["n_changed"], top_light["dist"])
        <= (best["n_changed"], best["dist"])
    ):
        best = evaluate(top_light["ref_idx"], top_light["tag"], keep_scores=True)

    # Write outputs
    ref_samples = merged.iloc[best["ref_idx"]]["sample"].astype(str).tolist()
    (out_dir / "ref40_samples.txt").write_text("\n".join(ref_samples) + "\n")
    pd.DataFrame({"sample": ref_samples}).to_csv(
        out_dir / "ref40_samples.tsv", sep="\t", index=False
    )

    score_df = _score_table(merged["sample"].tolist(), best["z"], best["e"], best["ez"])
    score_df.to_csv(out_dir / "ref40_score.tsv", sep="\t", index=False)

    meanstd = _meanstd_compare_table(arrays, early_idx, best["ref_idx"])
    meanstd.to_csv(out_dir / "reference_meanstd_compare.tsv", sep="\t", index=False)

    ref_mask = np.zeros(len(merged), dtype=bool)
    ref_mask[best["ref_idx"]] = True
    compare = base_compare & (~ref_mask)
    cmp_df = pd.DataFrame(
        {
            "sample": merged["sample"].astype(str),
            "label": merged["label"].astype(str),
            "set": merged["set"].astype(str),
            "old_ref_type": merged["ref_type"].astype(str),
            "in_early_ref": early_mask.to_numpy(),
            "in_ref40": ref_mask,
            "compared": compare,
            "pred_label_baseline": baseline_pred,
            "pred_label_ref40": best["pred"],
            "changed": baseline_pred != best["pred"],
            "meta_pred_label": merged["pred_label"].astype(str)
            if "pred_label" in merged.columns
            else "",
        }
    )
    cmp_df.to_csv(out_dir / "pred_label_compare.tsv", sep="\t", index=False)

    top = pd.DataFrame(
        [
            {
                "rank": i + 1,
                "tag": c["tag"],
                "n_changed": c["n_changed"],
                "n_compare": c["n_compare"],
                "meanstd_dist": c["dist"],
                "reference_list": c["reference_list"],
            }
            for i, c in enumerate(candidates[:50])
        ]
    )
    top.to_csv(out_dir / "search_top_candidates.tsv", sep="\t", index=False)

    overlap_early = len(set(ref_samples) & set(merged.iloc[early_idx]["sample"].astype(str)))
    summary = {
        "ref_n": ref_n,
        "n_random": n_random,
        "n_swap_rounds": n_swap_rounds,
        "seed": seed,
        "exclude_sample": exclude_sample,
        "pool_size": int(pool_idx.size),
        "early_ref_n": int(early_idx.size),
        "overlap_with_early_ref": overlap_early,
        "best_tag": best["tag"],
        "best_n_changed": int(best["n_changed"]),
        "best_n_compare": int(best["n_compare"]),
        "best_meanstd_dist": float(best["dist"]),
        "ezscore_cutoff": STRONG_CUTOFF,
        "ezscore_ref_n": int(ez_idx.size),
        "ref40_samples": ref_samples,
        "mean_abs_delta_mean": {
            feat: float(
                meanstd.loc[meanstd["feature"] == feat, "delta_mean"].abs().mean()
            )
            for feat in meanstd["feature"].unique()
        },
        "mean_abs_delta_std": {
            feat: float(
                meanstd.loc[meanstd["feature"] == feat, "delta_std"].abs().mean()
            )
            for feat in meanstd["feature"].unique()
        },
    }
    (out_dir / "selection_summary.json").write_text(json.dumps(summary, indent=2))

    console.print(f"[green]OK[/green] Best ref_40: n_changed={best['n_changed']} / {best['n_compare']}")
    console.print(f"  meanstd_dist={best['dist']:.6f} overlap_early={overlap_early}/{early_idx.size}")
    console.print(f"  Wrote outputs under {out_dir}")
    console.rule("[bold green]Done")


if __name__ == "__main__":
    try:
        main(standalone_mode=False)
    except click.ClickException as exc:
        console.print(f"[bold red]Error:[/bold red] {exc}")
        sys.exit(1)
