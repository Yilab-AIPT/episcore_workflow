#!/usr/bin/env python3
"""Render summary plots for select_stable_ref40 (no Jupyter required)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260730-stable_ref40")
CHR = [f"chr{i}" for i in range(1, 23)]


def main() -> None:
    meta = pd.read_csv(OUT / "meta.csv")
    beta = pd.read_csv(OUT / "beta.csv")
    pct = pd.read_csv(OUT / "percentage.csv", sep="\t")
    fetch = pd.read_csv(OUT / "zscore_fetch_report.tsv", sep="\t")
    meanstd = pd.read_csv(OUT / "reference_meanstd_compare.tsv", sep="\t")
    cmp = pd.read_csv(OUT / "pred_label_compare.tsv", sep="\t")
    summary = json.loads((OUT / "selection_summary.json").read_text())
    base = pd.read_csv(OUT / "baseline_score.tsv", sep="\t")
    new = pd.read_csv(OUT / "ref40_score.tsv", sep="\t")

    # episcore early_ref meanstd
    early_samples = set(meta.loc[meta.ref_type == "early_ref", "sample"].astype(str))
    rows = []
    for _, r in beta[beta["sample"].astype(str).isin(early_samples)].iterrows():
        for c in CHR:
            rows.append(
                {
                    "chr": c,
                    "hypo_z_intra": r[f"{c}_hypo_z_intra"],
                    "hyper_z_intra": r[f"{c}_hyper_z_intra"],
                    "hypo_cpgs_count": r[f"{c}_hypo_cpgs_count"],
                    "hyper_cpgs_count": r[f"{c}_hyper_cpgs_count"],
                }
            )
    early = pd.DataFrame(rows)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for ax, col, title in zip(
        axes.ravel(),
        ["hypo_z_intra", "hyper_z_intra", "hypo_cpgs_count", "hyper_cpgs_count"],
        ["hypo z_intra", "hyper z_intra", "hypo CpG count", "hyper CpG count"],
    ):
        stats = early.groupby("chr")[col].agg(["mean", "std"]).reindex(CHR)
        ax.errorbar(range(22), stats["mean"], yerr=stats["std"], fmt="o-", capsize=2)
        ax.set_title(f"early_ref (n=17): {title}")
        ax.set_xticks(range(22))
        ax.set_xticklabels([c.replace("chr", "") for c in CHR], fontsize=8)
    fig.suptitle("Episcore source — early_ref mean±std by chromosome")
    fig.tight_layout()
    fig.savefig(OUT / "plot_episcore_early_ref_meanstd.png", dpi=150)
    plt.close(fig)

    # zscore percentage early_ref
    pct_m = pct.merge(meta[["sample", "ref_type"]], on="sample", how="left")
    early_pct = pct_m[pct_m.ref_type == "early_ref"]
    stats = early_pct.groupby("chr")["percentage"].agg(["mean", "std"]).reindex(CHR)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.errorbar(range(22), stats["mean"], yerr=stats["std"], fmt="o-", capsize=2, color="#4C78A8")
    ax.set_xticks(range(22))
    ax.set_xticklabels([c.replace("chr", "") for c in CHR])
    ax.set_xlabel("chr")
    ax.set_ylabel("percentage (%)")
    ax.set_title("Zscore source — early_ref mean±std percentage (cutoff=0.85)")
    fig.tight_layout()
    fig.savefig(OUT / "plot_zscore_early_ref_percentage_meanstd.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    fetch["status"].value_counts().plot(kind="bar", ax=ax, color="#F58518")
    ax.set_title("Missing zscore samples — fetch result")
    ax.set_ylabel("n samples")
    fig.tight_layout()
    fig.savefig(OUT / "plot_zscore_fetch_status.png", dpi=150)
    plt.close(fig)

    # meanstd compare
    features = [
        "episcore_hypo_z_intra",
        "episcore_hyper_z_intra",
        "zscore_percentage",
    ]
    fig, axes = plt.subplots(len(features), 2, figsize=(12, 10), sharex=True)
    x = np.arange(22)
    for i, feat in enumerate(features):
        sub = meanstd[meanstd.feature == feat].set_index("chr").reindex(CHR)
        axes[i, 0].plot(x, sub["early_ref_mean"], "o-", label="early_ref_17", color="#4C78A8")
        axes[i, 0].plot(x, sub["ref40_mean"], "s--", label="ref_40", color="#F58518")
        axes[i, 0].set_ylabel("mean")
        axes[i, 0].set_title(f"{feat} mean")
        axes[i, 0].legend(fontsize=8)
        axes[i, 1].plot(x, sub["early_ref_std"], "o-", label="early_ref_17", color="#4C78A8")
        axes[i, 1].plot(x, sub["ref40_std"], "s--", label="ref_40", color="#F58518")
        axes[i, 1].set_ylabel("std")
        axes[i, 1].set_title(f"{feat} std")
        axes[i, 1].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("chr", "") for c in CHR], fontsize=8)
        ax.set_xlabel("chr")
    fig.suptitle("Reference mean/std: early_ref_17 vs selected ref_40")
    fig.tight_layout()
    fig.savefig(OUT / "plot_ref40_vs_early_meanstd.png", dpi=150)
    plt.close(fig)

    compared = cmp[cmp.compared]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(
        ["unchanged", "changed"],
        [(~compared.changed).sum(), compared.changed.sum()],
        color=["#54A24B", "#E45756"],
    )
    ax.set_ylabel("n samples")
    ax.set_title(f'Ezscore pred_label stability (cutoff={summary["ezscore_cutoff"]})')
    fig.tight_layout()
    fig.savefig(OUT / "plot_pred_label_stability.png", dpi=150)
    plt.close(fig)

    m = base.merge(new, on="sample", suffixes=("_base", "_new"))
    m = m.merge(cmp[["sample", "compared"]], on="sample")
    eval_m = m[m.compared]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, prefix, title in zip(
        axes,
        ["episcore", "zscore", "ezscore"],
        ["episcore (beta)", "zscore (rc)", "ezscore (final)"],
    ):
        cols_b = [f"{prefix}_chr{i}_base" for i in range(1, 23)]
        cols_n = [f"{prefix}_chr{i}_new" for i in range(1, 23)]
        xb = eval_m[cols_b].to_numpy().ravel()
        xn = eval_m[cols_n].to_numpy().ravel()
        mask = np.isfinite(xb) & np.isfinite(xn)
        ax.scatter(xb[mask], xn[mask], s=4, alpha=0.15, color="#4C78A8")
        lim = np.nanpercentile(np.concatenate([xb[mask], xn[mask]]), [0.5, 99.5])
        ax.plot(lim, lim, "--", color="gray", lw=1)
        r = np.corrcoef(xb[mask], xn[mask])[0, 1]
        ax.set_title(f"{title}\nPearson r={r:.4f}")
        ax.set_xlabel("early_ref_17")
        ax.set_ylabel("ref_40")
    fig.suptitle("Per-chr scores on shared analyze samples")
    fig.tight_layout()
    fig.savefig(OUT / "plot_score_correlation_baseline_vs_ref40.png", dpi=150)
    plt.close(fig)

    print(f"Wrote plots under {OUT}")


if __name__ == "__main__":
    main()
