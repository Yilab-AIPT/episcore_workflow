#!/usr/bin/env python3
"""Search subsets of cartesian ez pairs that improve Normal/Trisomy separation.

Uses aggregated ``pair_abnormal_count.npz`` from a filtered ``ref_free_ezscore``
run (``--store-pair-counts``).

Strategy
--------
1. Score each pair alone (AUC on ff>=1% eval samples at a target ez cutoff).
2. Sweep top-K unions (K in a ladder).
3. Greedy forward selection: repeatedly add the pair that most improves AUC.
4. Write best subset + reweighted signal ratios for that subset.
"""

from __future__ import annotations

import json
from pathlib import Path

import click
import numpy as np
import pandas as pd
from rich.console import Console

from separation import is_trisomy_label, roc_auc, separation_index

console = Console()


def _ez_ratio_col(cutoff: float) -> str:
    return f"ezscore_signal_ratio_{cutoff:g}"


def _ff_eval_mask(df: pd.DataFrame, ff_min: float) -> np.ndarray:
    ff = pd.to_numeric(df["ff_before_mq"], errors="coerce").to_numpy()
    labels = df["label"].astype(str)
    keep = (ff >= ff_min) & (
        labels.eq("Normal") | labels.map(is_trisomy_label).to_numpy()
    )
    # exclude val from subset search (tune on eval only)
    if "set" in df.columns:
        keep &= df["set"].astype(str).ne("val").to_numpy()
    return keep.to_numpy() if hasattr(keep, "to_numpy") else np.asarray(keep)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--output-base", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--ff-min", default=0.01, show_default=True, type=float)
@click.option("--target-ez-cutoff", default=3.0, show_default=True, type=float)
@click.option("--max-greedy", default=40, show_default=True, type=int)
@click.option(
    "--top-k-ladder",
    default="1,5,10,20,40,80,160",
    show_default=True,
    help="Comma-separated K values for top-K union sweep",
)
def main(
    output_base: str,
    ff_min: float,
    target_ez_cutoff: float,
    max_greedy: int,
    top_k_ladder: str,
) -> None:
    out_root = Path(output_base) / "ref_free_ezscore"
    pair_path = out_root / "pair_abnormal_count.npz"
    scores_path = out_root / "abnormality_signal_ratio.tsv"
    config_path = out_root / "run_config.json"
    if not pair_path.is_file():
        raise click.ClickException(
            f"Missing {pair_path}; rerun with --store-pair-counts"
        )
    config = json.loads(config_path.read_text())
    df = pd.read_csv(scores_path, sep="\t")
    data = np.load(pair_path)
    pair_counts = data["pair_abnormal_count"]  # [n_pairs, n_cutoffs, n_eval]
    ez_cutoffs = [float(x) for x in data["ez_cutoffs"]]
    n_repeats = int(data["n_repeats"][0]) if "n_repeats" in data.files else int(
        config["total_repeats"]
    )
    if target_ez_cutoff not in ez_cutoffs:
        raise click.ClickException(
            f"target-ez-cutoff {target_ez_cutoff} not in {ez_cutoffs}"
        )
    c_idx = ez_cutoffs.index(target_ez_cutoff)
    n_pairs = pair_counts.shape[0]
    if n_pairs != int(config["n_ez_combos"]):
        console.print(
            f"[yellow]pair n={n_pairs} vs config n_ez_combos="
            f"{config['n_ez_combos']}[/yellow]"
        )

    mask = _ff_eval_mask(df, ff_min)
    y = df.loc[mask, "label"].map(is_trisomy_label).to_numpy()
    if y.sum() == 0 or (~y).sum() == 0:
        raise click.ClickException("Need both Normal and Trisomy with ff>=ff_min on eval")

    # per-pair ratios on masked samples
    # pair_counts[:, c_idx, :] / n_repeats  (single pair denom = 1 * n_repeats)
    single_ratios = pair_counts[:, c_idx, :][:, mask].astype(np.float64) / float(n_repeats)
    single_auc = np.array([roc_auc(single_ratios[i], y) for i in range(n_pairs)])
    order = np.argsort(-np.nan_to_num(single_auc, nan=-1.0))

    def auc_for_subset(idxs: list[int]) -> float:
        if not idxs:
            return float("nan")
        tot = pair_counts[idxs][:, c_idx, :][:, mask].sum(axis=0).astype(np.float64)
        ratios = tot / (float(len(idxs)) * float(n_repeats))
        return roc_auc(ratios, y)

    baseline_all = auc_for_subset(list(range(n_pairs)))
    trials = [
        {
            "name": "all_pairs",
            "n_pairs": n_pairs,
            "sep": baseline_all,
            "pair_indices": list(range(n_pairs)),
        }
    ]

    ladder = [int(x) for x in top_k_ladder.split(",") if x.strip()]
    for k in ladder:
        k = min(k, n_pairs)
        idxs = order[:k].tolist()
        trials.append(
            {
                "name": f"top_{k}",
                "n_pairs": k,
                "sep": auc_for_subset(idxs),
                "pair_indices": idxs,
            }
        )

    # greedy forward
    selected: list[int] = []
    remaining = set(range(n_pairs))
    greedy_trace = []
    for step in range(min(max_greedy, n_pairs)):
        best_i, best_auc = None, -1.0
        for i in remaining:
            auc = auc_for_subset(selected + [i])
            if np.isfinite(auc) and auc > best_auc:
                best_auc, best_i = auc, i
        if best_i is None:
            break
        selected.append(best_i)
        remaining.remove(best_i)
        greedy_trace.append({"step": step + 1, "added": best_i, "sep": best_auc})
    if selected:
        trials.append(
            {
                "name": f"greedy_{len(selected)}",
                "n_pairs": len(selected),
                "sep": auc_for_subset(selected),
                "pair_indices": list(selected),
            }
        )

    trials_sorted = sorted(
        trials, key=lambda t: (-(t["sep"] if np.isfinite(t["sep"]) else -1), t["n_pairs"])
    )
    best = trials_sorted[0]
    console.print(
        f"[green]Best subset[/green] {best['name']} sep={best['sep']:.4f} "
        f"n_pairs={best['n_pairs']} (all_pairs sep={baseline_all:.4f})"
    )

    # Materialize signal ratios for best subset across all ez cutoffs
    idxs = best["pair_indices"]
    denom = float(len(idxs) * n_repeats)
    result = df.copy()
    for j, c in enumerate(ez_cutoffs):
        tot = pair_counts[idxs][:, j, :].sum(axis=0).astype(np.float64)
        result[_ez_ratio_col(c)] = tot / denom
        result[f"ezscore_abnormal_count_{c:g}"] = tot.astype(np.int64)
    if 3.0 in ez_cutoffs:
        result["ezscore_signal_ratio"] = result[_ez_ratio_col(3.0)]
        result["ezscore_abnormal_count"] = result[f"ezscore_abnormal_count_3"]

    best_path = out_root / "abnormality_signal_ratio_best_subset.tsv"
    result.to_csv(best_path, sep="\t", index=False, float_format="%.6f")

    # Separation summary for best subset
    is_val = result["set"].astype(str).eq("val") if "set" in result.columns else np.zeros(len(result), dtype=bool)
    sep_eval = separation_index(result[~is_val], _ez_ratio_col(target_ez_cutoff), ff_min=ff_min)
    sep_val = (
        separation_index(result[is_val], _ez_ratio_col(target_ez_cutoff), ff_min=ff_min)
        if is_val.any()
        else {}
    )

    report = {
        "target_ez_cutoff": target_ez_cutoff,
        "ff_min": ff_min,
        "n_repeats": n_repeats,
        "n_pairs_total": n_pairs,
        "baseline_all_sep": baseline_all,
        "best": best,
        "trials": [
            {k: v for k, v in t.items() if k != "pair_indices"}
            | {"n_pair_indices": len(t["pair_indices"])}
            for t in trials_sorted
        ],
        "greedy_trace": greedy_trace,
        "separation_eval_best": sep_eval,
        "separation_val_best": sep_val,
        "best_subset_tsv": str(best_path),
    }
    # keep full best indices
    report["best"] = best
    report_path = out_root / "ez_pair_subset_search.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    # human-readable ranking of single pairs
    rank_df = pd.DataFrame(
        {
            "pair_index": np.arange(n_pairs),
            "sep_auc": single_auc,
        }
    ).sort_values("sep_auc", ascending=False)
    if "ez_pairs" in config and "ep_combos" in config and "z_combos" in config:
        rows = []
        for i in range(n_pairs):
            ei, zi = config["ez_pairs"][i]
            ep = config["ep_combos"][ei]
            zc = config["z_combos"][zi]
            rows.append(
                {
                    "pair_index": i,
                    "ep_threshold": ep[0],
                    "ep_recall": ep[1],
                    "z_threshold": zc[0],
                    "z_recall": zc[1],
                    "sep_auc": float(single_auc[i]),
                }
            )
        rank_df = pd.DataFrame(rows).sort_values("sep_auc", ascending=False)
    rank_df.to_csv(out_root / "ez_pair_single_auc.tsv", sep="\t", index=False, float_format="%.6f")

    console.print(f"[green]OK[/green] Wrote {report_path}")
    console.print(f"  -> {best_path}")


if __name__ == "__main__":
    main()
