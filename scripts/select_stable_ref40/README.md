# Select stable ref_40

Offline toolkit to (1) summarize episcore / zscore calculation sources, (2) fetch
missing zscore percentages, and (3) choose a Normal+dev `ref_40` that preserves
early_ref mean/std and ezscore `pred_label` (cutoff 4.5).

## Inputs

| Input | Path |
|-------|------|
| Episcore samplesheet | `/lustre1/cqyi/syfan/nipt_article_plot/episcore_result_samplesheet.csv` |
| Percentage (partial) | `/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260607-ref_40/percentage.csv` |
| Meta samplesheet | `/lustre1/cqyi/syfan/nipt_article_plot/temporary_updated_samplesheet.csv` |
| Ezscore refs (n=25) | `/lustre1/cqyi/myli/bert/analysis_nipt/multiomics/chr_stats_reference_samples.txt` |

Combos:

- Episcore: early.config recall=0.65, threshold=0.5
- Zscore: recall=0.95, cutoff=0.85

## Run

```bash
bash /lustre1/cqyi/AIPT_2.0/workflow/episcore/scripts/select_stable_ref40/run_all.sh

# or step-wise via singularity + common_tools.sif:
# build_source_tables.py → select_ref40.py → write_updated_meta.py → render_plots.py
```

Useful knobs for `select_ref40.py`:

- `--n-random` / `--n-swap-rounds` / `--seed`
- `--exclude-sample ''` to allow all Normal+dev (default historically excluded `PTAY0586P8S1`)

## Outputs

Directory: `/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260730-stable_ref40/`

| File | Description |
|------|-------------|
| `beta.csv` | Episcore source table (all meta samples) |
| `percentage.csv` | Zscore percentages at cutoff=0.85 |
| `zscore_fetch_report.tsv` | Per-missing-sample fetch status |
| `ref40_samples.txt` | Selected ref_40 list |
| `baseline_score.tsv` / `ref40_score.tsv` | Recalculated scores |
| `reference_meanstd_compare.tsv` | early_ref vs ref_40 mean/std |
| `pred_label_compare.tsv` | Pred-label diffs (ezscore cutoff 4.5) |
| `selection_summary.json` | Search metrics |
| `temporary_updated_samplesheet_ref40.csv` | Updated meta (`ref_type` + `*_zscores` + `pred_label`) |
| `plot_*.png` | Summary figures |

Also copied: `/lustre1/cqyi/syfan/nipt_article_plot/temporary_updated_samplesheet_ref40.csv`

## Notebooks

- `notebooks/aipt_2.0/summarize_score_sources.ipynb`
- `notebooks/aipt_2.0/select_stable_ref40.ipynb`
