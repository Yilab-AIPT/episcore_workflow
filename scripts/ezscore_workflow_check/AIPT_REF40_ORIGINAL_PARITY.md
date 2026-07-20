# aipt_ref_40 parity with original Snakemake ezscore workflow

This document records changes made so **`params.step=aipt_ref_40`** reproduces the same
episcore / zscore / ezscore as the original Snakemake stack (`dna_5mc_analysis_pipeline`
+ `workflow/zscore` + ezscore normalization), under the **fixed early combo**:

| score    | read-prob threshold | recall |
|----------|---------------------|--------|
| episcore | 0.5                 | 0.65   |
| zscore   | 0.85                | 0.95   |

Use this when adapting **other launch scripts** (batch jobs, offline checks, notebooks)
that run `aipt_ref_40`. The reference implementation lives in this directory:
`run_aipt_ref40_check.slurm`, `build_original_combo_assets.py`, `compare_with_meta.py`.

---

## 1. Problem that was fixed

Before these changes, `aipt_ref_40` episcore **did not match NIPT / original Snakemake**
because depth filtering used **target-only** coverage (`meth + unmeth` from the target
bedGraph). The original pipeline filters on **`raw_total = target + background`** read
counts at the split threshold (`0.5`).

Everything else (hypo/hyper z_intra → s_inter formula, zscore percentage logic, ezscore
normalization) was already aligned when given the correct reference assets.

---

## 2. Episcore: raw_total depth filter (main code change)

### Original (NIPT / Snakemake)

1. Split clean BAM at **`threshold=0.5`** → **target** (`prob ≥ 0.5`) + **background** (`prob < 0.5`).
2. MethylDackel on **both** BAMs.
3. Merge bedGraphs → `raw_total_count = target_meth + target_unmeth + background_meth + background_unmeth`.
4. Keep CpGs with **`raw_total_count > beta_depth_threshold`** (default **30**).
5. Beta for episcore still comes from **target** meth/unmeth only.

### aipt_ref_40 (after fix)

Same logic, wired through flexible per-threshold episcore:

| step | module / script |
|------|-----------------|
| Split target BAMs for all episcore combo thresholds **plus** background at `params.threshold` | `SPLIT_BAM_BY_THRESHOLDS` + `split_reads_by_thresholds.py` |
| MethylDackel on all target-threshold BAMs | `METHYLDACKEL_TARGET` in `subworkflows/local/aipt_episcore.nf` |
| MethylDackel on background BAM at `params.threshold` | `METHYLDACKEL_BACKGROUND` |
| Build depth-passing CpG list from thr=`params.threshold` target+bg bedGraphs | **`FILTER_CPG_BY_DEPTH`** + **`bin/filter_cpg_by_depth.py`** |
| Episcore per chromosome using combo-specific target bedGraph + recall list **intersected with depth list** | `CALC_EPISCORE_FLEXIBLE` + `bin/calc_episcore_flexible.py` |

### New / modified files

```
bin/filter_cpg_by_depth.py              # NEW — outer-merge target+bg, raw_total > depth
bin/split_reads_by_thresholds.py        # --background-threshold (emit background_thr_{t}.txt)
bin/calc_episcore_flexible.py           # --depth-filtered-cpgs (replaces --depth target-only filter)
modules/local/filter_cpg_by_depth/main.nf   # NEW process
modules/local/split_bam_by_thresholds/main.nf # +background_threshold input, background_bams output
modules/local/calc_episcore_flexible/main.nf  # takes depth_filtered_cpgs, no longer passes --depth
subworkflows/local/aipt_episcore.nf     # full wiring (see below)
conf/alioth_slurm.config                # FILTER_CPG_BY_DEPTH container/resources
conf/alioth_local.config                # same
```

**Important:** `bin/filter_cpg_by_depth.py` must be **executable** (`chmod +x`). Nextflow
invokes it directly; missing `+x` causes exit **126 Permission denied**.

### Workflow wiring (`subworkflows/local/aipt_episcore.nf`)

1. **`depth_thr`** = `params.threshold` formatted as `%g` (typically `"0.5"`).
2. Ensure `depth_thr` is in the split threshold list even if absent from
   `best_combo_episcore.csv` (needed for depth bedGraph pair).
3. `SPLIT_BAM_BY_THRESHOLDS(ch_samplesheet, split_thresholds_str, params.threshold)`.
4. Target path unchanged: index → MethylDackel → `groupTuple` bedGraphs per sample.
5. Background path: optional `background_bams` → index → MethylDackel → join with
   **target bedGraph at `depth_thr` only** → `FILTER_CPG_BY_DEPTH`.
6. `CALC_EPISCORE_FLEXIBLE(ch_bedgraphs.join(filtered_cpgs), ...)` — depth list is
   **shared across all combo thresholds** for that sample (same as original: one
   depth filter at split threshold 0.5).

### Params that must match original

```groovy
// conf/aipt_ref_40.config (or your profile)
threshold              = 0.5    // BAM split for depth filter AND SNP-FF when FF enabled
beta_depth_threshold   = 30    // raw_total > 30
cpg_recall_dir         = .../assets/220k_cpg_recall_list   // 220k_cpg_recall_{recall}.txt
```

For a **fixed-combo** reproduction run, point `--grid_search_result` at assets where
**every chromosome** has the same `(threshold, recall)` — see §3.

---

## 3. Reference assets: match original mean/std (not ref-40 grid-search output)

`aipt_ref_40` normally consumes **ref-40 grid-search** outputs
(`best_combo_episcore.csv`, `best_reference_matrix.tsv`, …). To match the **original
Snakemake combo**, build a synthetic asset directory with
`scripts/ezscore_workflow_check/build_original_combo_assets.py`:

```bash
python3 build_original_combo_assets.py --output-dir /path/to/grid_search_assets
```

This writes:

| file | source / note |
|------|----------------|
| `best_combo_episcore.csv` | fixed **0.5 / 0.65** all chr |
| `best_combo_zscore.csv` | fixed **0.85 / 0.95** all chr |
| `best_reference_matrix.tsv` | episcore stats from `assets/early_reference_beta_zscore.tsv` (**ddof=0**); zscore stats from production `igtc_220k.early17...reference.csv` (**ddof=1**, matches `mq_zscore_analyzer.py` / polars) |
| `best_ezscore_ref_20_matrix.tsv` | from `chr_stats.csv` (mu/sigma of episcore+zscore) |
| `best_sample_scores_recalc_ezscore.tsv` | empty stub so `hasEzscoreAssets()` is true and plots run |
| `LOGIC_NOTES.md` | short parity notes (regenerated each build) |

Pass to Nextflow:

```bash
nextflow run .../main.nf \
  -profile aipt_ref_40,alioth_slurm,singularity \
  --step aipt_ref_40 \
  --grid_search_result /path/to/grid_search_assets \
  ...
```

**Do not** use the default `grid_search_result` path in `conf/aipt_ref_40.config` when
the goal is original-combo parity — that directory holds **per-chr flexible combos** from
ref-40 grid search, not the fixed early production combo.

---

## 4. Comparison-only runs: skip SNP-FF (`skip_ff`)

Meta ground-truth scores (`beta_zscores`, `rc_zscores`, `final_zscores`) do **not**
require fetal fraction. For parity checks, run with:

```bash
--skip_ff true
```

Changes tied to this flag:

| file | change |
|------|--------|
| `nextflow.config`, `conf/aipt_ref_40.config` | `params.skip_ff` |
| `workflows/aipt_ref_40.nf` | skips `AIPT_FF`; passes empty FF channel to merge |
| `modules/local/merge_scores/main.nf` | passes `--skip-ff` to script |
| `bin/merge_scores.py` | omits `ff_before_mq` column when `--skip-ff` |
| `bin/plot_scores.py` | tolerates missing `ff_before_mq` (uses placeholder x-axis) |

Production scoring should keep **`skip_ff=false`** (default) so `{sample}_scores.tsv`
includes `ff_before_mq`.

---

## 5. Zscore and ezscore (unchanged logic; asset-sensitive)

No new processes were added for zscore/ezscore parity. Alignment depends on:

- **Zscore:** `AIPT_ZSCORE` / `calc_zscore_flexible.py` with `best_combo_zscore.csv`
  at **0.85/0.95** and reference `percentage_mean/std` built with **ddof=1**.
- **Ezscore:** `merge_scores.py` computes `(episcore + zscore - mean) / std` using
  `best_ezscore_ref_20_matrix.tsv` from `chr_stats.csv`.
- **CpG overlap:** pybedtools (original) vs searchsorted on CpG start (Nextflow) —
  should agree for single-base CpG sites.

---

## 6. Checklist for adapting another aipt_ref_40 launch script

1. **Profile:** `-profile aipt_ref_40,alioth_slurm,singularity` (or local equivalent).
2. **Assets:** run `build_original_combo_assets.py` **or** ensure your
   `--grid_search_result` uses the same reference sources and fixed combos above.
3. **Params:** `threshold=0.5`, `beta_depth_threshold=30`, correct `cpg_recall_dir`.
4. **Samplesheet:** `[sample, clean_bam, deconv_res]` — same as NIPT; use
   `PREPARE_INPUTS` merge if multiple rows per sample.
5. **Comparison vs meta:** `--skip_ff true`; read scores from
   `{outdir}/merge_scores/{sample}_scores.tsv`; compare to meta columns
   `beta_zscores` / `rc_zscores` / `final_zscores` (see `compare_with_meta.py`).
6. **Resume:** `-work-dir ... -resume` after code fixes; rerun episcore stages if
   depth-filter logic changed (cached `CALC_EPISCORE_FLEXIBLE` from old target-only
   depth is invalid).
7. **New bin script:** after pulling, `chmod +x bin/filter_cpg_by_depth.py`.

### Minimal Nextflow invocation (parity check)

```bash
nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
  -profile aipt_ref_40,alioth_slurm,singularity \
  -work-dir /path/to/work_pipeline \
  -resume \
  --step aipt_ref_40 \
  --input /path/to/check_samplesheet.csv \
  --outdir /path/to/pipeline \
  --grid_search_result /path/to/grid_search_assets \
  --skip_ff true
```

---

## 7. Known remaining differences (usually negligible)

| topic | original | aipt_ref_40 after fix |
|-------|----------|------------------------|
| Depth filter site set | outer merge can retain sites with background-only coverage passing depth | intersect with target bedGraph + recall list (sites without target methylation never contribute beta) |
| Episcore per-chr threshold | single split at 0.5 | flexible: each chr can use its combo threshold bedGraph, **same** depth list |
| Zscore CpG overlap | pybedtools intersect | start-position searchsorted |
| FF | always computed in production Snakemake path | optional; skipped with `--skip_ff true` for score-only checks |

---

## 8. Files touched (summary for code search)

**Episcore depth parity**

- `subworkflows/local/aipt_episcore.nf`
- `modules/local/split_bam_by_thresholds/main.nf`
- `modules/local/filter_cpg_by_depth/main.nf` *(new)*
- `modules/local/calc_episcore_flexible/main.nf`
- `bin/split_reads_by_thresholds.py`
- `bin/filter_cpg_by_depth.py` *(new)*
- `bin/calc_episcore_flexible.py`
- `conf/alioth_slurm.config`, `conf/alioth_local.config`

**FF-optional / plotting (ezscore check)**

- `workflows/aipt_ref_40.nf`
- `modules/local/merge_scores/main.nf`
- `bin/merge_scores.py`
- `bin/plot_scores.py`
- `nextflow.config`, `conf/aipt_ref_40.config`

**Check harness (this directory)**

- `build_original_combo_assets.py` — builds parity assets + `LOGIC_NOTES.md`
- `prepare_check_samples.py`, `compare_with_meta.py`
- `run_aipt_ref40_check.slurm`, `submit_aipt_ref40_check.sh`

---

## 9. Validation outputs

After a successful check run:

```
{output_dir}/pipeline/merge_scores/{sample}_scores.tsv   # rebuilt scores
{output_dir}/comparison/check_comparison_detail.tsv
{output_dir}/comparison/check_comparison_summary.tsv
```

Notebook: `notebooks/aipt_2.0/ezscore_check.ipynb`.

Meta mapping (`compare_with_meta.py`):

| pipeline column | meta.csv column |
|-----------------|-----------------|
| `episcore`      | `beta_zscores`  |
| `zscore`        | `rc_zscores`    |
| `ezscore`       | `final_zscores` |

Vectors in meta are comma-separated chr1–chr22 values.
