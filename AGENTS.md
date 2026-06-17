# Episcore — Agent Guide

## What this project is

**Episcore** is a Nextflow pipeline for **NIPT trisomy detection** using methylation **episcore** (beta → chromosome-level z-scores) plus **SNP-based fetal fraction (FF)** estimation. It sits in the AIPT 2.0 workflow at `/lustre1/cqyi/AIPT_2.0/workflow/episcore`.

Upstream work (outside this repo) produces **clean BAMs** and **deconvolution read-probability tables**; this pipeline splits reads by deconv threshold, extracts methylation (MethylDackel), computes beta/z-scores against a reference matrix, estimates FF from SNPs, and emits per-sample reports.

## Entry points

| `params.step` | Workflow | Samplesheet columns |
|---------------|----------|---------------------|
| `split_bam` | `NIPT` → `SPLIT_BAM` | `sample`, `clean_bam`, `deconv_res` |
| `beta_zscore` | `NIPT` → `CALC_BETA_ZSCORE`, `ESTIMATE_FF`, `REPORT` | `sample`, `target_bam`, `background_bam` (or `clean_bam`+`deconv_res` for `split_bam`) |
| `grid_search` | `GRID_SEARCH` → `EXTRACT_BETA` | `sample`, `clean_bam`, `deconv_res` |
| `est_ff_from_bam` | `SNP_EST_FF` → `SPLIT_BAM`, `BAM_TO_PILEUP`, `ESTIMATE_FF_HIGHER_PRECISION` | `sample`, `clean_bam`, `deconv_res` |
| `est_ff_from_pileup` | `SNP_EST_FF` → `ESTIMATE_FF_HIGHER_PRECISION` | `sample`, `pileup` |
| `perturbed_res` | `PERTURBED_RES` → `MERGE_DECONV_RES`, `REPLACE_DECONV_PROB`, `SPLIT_BAM_BY_DECONV_RES`, `CALC_BETA_ZSCORE`, `ESTIMATE_FF` | `sample`, `full_name`, `clean_bam`, `perturbed_res`, `original_res` |

`main.nf` routes by `params.step`: NIPT steps use `validateAndParseSamplesheet`; grid search uses `validateAndParseGridSearchParameters`; SNP FF steps (`est_ff_from_bam` / `est_ff_from_pileup`) use `validateAndParseSnpFFSamplesheet`; the methylation-perturbation step (`perturbed_res`) uses `validateAndParsePerturbedResSamplesheet` (`lib/`).

The **NIPT FF** subworkflow (`subworkflows/local/estimate_ff.nf`, used by `beta_zscore` and `perturbed_res`) runs `BAM_TO_PILEUP` then branches on `params.ff_precision`: when set (e.g. `0.001`), it calls `ESTIMATE_FF_HIGHER_PRECISION` (`bin/estimate_ff_with_higher_precision.py`, passes `--ff-precision`); when `null` (default), it calls `SNP_TO_FF` (`bin/estimate_ff.py`, fixed-step grid search).

The **SNP FF** workflow (`workflows/snp_est_ff.nf`, entry `EST_FF`) is a standalone path focused only on SNP-based fetal fraction. It reuses `SPLIT_BAM` + `BAM_TO_PILEUP` to build a pileup (`est_ff_from_bam`) or consumes a pre-computed pileup (`est_ff_from_pileup`), then always runs `ESTIMATE_FF_HIGHER_PRECISION` (requires `--ff-precision`). It does **not** go through `CALC_BETA_ZSCORE` / `REPORT`.

The **methylation-perturbation** workflow (`workflows/perturbed_res.nf`, entry `PERTURB`, `params.step=perturbed_res`) measures how perturbing read methylation status changes downstream episcore/FF. Per sample it (1) merges the (shared) `original_res` files via `MERGE_DECONV_RES`, then for each perturbation condition (`full_name`, `{sample}_*`) (2) overwrites the original `prob_class_1` with the perturbed read probabilities via `REPLACE_DECONV_PROB` (`bin/replace_deconv_prob.py`, left-join/coalesce on read `name`), (3) splits the single clean BAM with `SPLIT_BAM_BY_DECONV_RES` directly — **no `SAMTOOLS_MERGE` / `PICARD_MARKDUPLICATES`** since there is one BAM per sample — and (4) runs `CALC_BETA_ZSCORE` + `ESTIMATE_FF`. It does **not** go through `REPORT`. `meta = [id: full_name, sample: sample]`, so all per-condition outputs are keyed by `full_name`.

## Directory map

```
main.nf                 # Top-level router (MAIN / SUB / EST_FF / PERTURB workflows)
nextflow.config         # params, profiles, workDir, reports
workflows/              # nipt.nf, grid_search.nf, snp_est_ff.nf, perturbed_res.nf
subworkflows/local/     # split_bam, calc_beta_zscore, estimate_ff, report, extract_beta
modules/local/          # Project-specific processes (incl. replace_deconv_prob)
modules/nf-core/        # samtools, picard, methyldackel (vendored)
lib/                    # Groovy samplesheet / grid-search / snp-ff / perturbed-res parsers
bin/                    # Python CLI scripts invoked by processes
conf/                   # Profile-specific params + executor/container config
assets/                 # Reference FASTA, CpG/SNP lists, reference z-score matrices (not in git)
containers/             # Singularity images (gitignored); default Python env: common_tools.sif
scripts/                # Offline grid-search / reference exploration (not part of main.nf)
run_grid_search*.sh     # Batch nextflow launches on lustre
notebooks/              # Analysis notebooks (gitignored)
```

## Typical run (production)

```bash
nextflow run /lustre1/cqyi/AIPT_2.0/workflow/episcore/main.nf \
  -profile early,alioth_slurm,singularity \
  --input /path/to/samplesheet.csv \
  --outdir /lustre1/cqyi/AIPT_2.0/results/episcore_output/<run_id> \
  --step split_bam   # or beta_zscore | est_ff_from_bam | est_ff_from_pileup | perturbed_res
```

- **workDir** (fixed in `nextflow.config`): `/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir`
- **Profiles**: `early`, `early_240k`, `middle`, `filter_size`, `at_analyze`, `at_ref`, `grid_search`, `perturbed_res`, `test`, plus executors `alioth_slurm`, `alioth_local`, `dev`, `singularity`
- Panel variants differ mainly in `cpg_list`, `snp_list`, `reference_beta_zscore_matrix` under `conf/*.config`
- **Containers**: local/Python processes use `containers/common_tools.sif` (default Python environment for `bin/*.py`); `METHYLDACKEL.*` → `methyldackel.sif`, `PICARD_*` → `picard_3.4.0.sif` (see `conf/alioth_slurm.config`)
- **`ff_precision`**: optional on NIPT (`--ff_precision 0.001`); when set, `ESTIMATE_FF` uses higher-precision FF estimation instead of `estimate_ff.py`

## Channel conventions

- **meta**: `[id: sample]` from samplesheet `sample` column (for `perturbed_res`: `[id: full_name, sample: sample]`)
- **split_bam input**: `[meta, clean_bam, deconv_res]`
- **beta_zscore input**: `[meta, target_bam, background_bam]`
- **est_ff_from_bam input**: `[meta, clean_bam, deconv_res]`
- **est_ff_from_pileup input**: `[meta, pileup]`
- **perturbed_res input**: `[meta, clean_bam, perturbed_res, original_res]` (meta carries `id`=full_name + `sample`)
- Processes use `tag "$meta.id"`; outputs often prefixed with `${meta.id}`

## Python ↔ Nextflow contract

- Scripts in `bin/` are called by name from process `script:` blocks (on `PATH` via Nextflow `bin/`).
- CLI: **click** + **rich**; tabular IO: **pandas** / **polars** (see `filter_deconv_res.py`, `merge_deconv_res.py`).
- `estimate_ff.py` imports `FFEstimator` from `bin/FFEstimator.py` (same directory in container).
- `estimate_ff_with_higher_precision.py` (NIPT when `params.ff_precision` is set; always in SNP FF workflow) reuses `FFEstimator` + `load_and_validate_data`/`parse_*` from `estimate_ff.py`; iterative range-narrowing search via `--ff-precision` ← `params.ff_precision`. `--mode-list` ← `params.snp_est_mode`, `--min-raw-depth` ← `params.snp_depth_threshold`, optional `--known-sites` ← `params.snp_list` (shared with `bam_to_pileup.py`, filters the pileup to panel sites before estimation). It is the offline `scripts/ff_decimal/` variant adapted to a single pileup per process.
- `replace_deconv_prob.py` (perturbed_res workflow) overwrites the merged original `prob_class_1` with the perturbed read probabilities via a polars left-join + `coalesce` on read `name` (perturbed reads are a subset of the original), emitting a `name, prob_class_1, insert_size` parquet for `SPLIT_BAM_BY_DECONV_RES`. Uses the streaming engine for the large (~10^8-row) original files.

## Agent do / don't

**Do**

- Match existing patterns: `include { MODULE } from '...'`, subworkflow `take` / `main` / `emit`, `params.*` for config
- Put reusable parsing in `lib/*.nf`; new processes under `modules/local/<name>/main.nf`
- Keep changes minimal; read neighboring modules before editing channel joins
- Use lustre paths for data/results; respect `.gitignore` (no `assets/`, `containers/`, `notebooks/`, `results/`)

**Don't**

- Commit secrets, large binaries, or generated logs under `scripts/**/logs`
- Change `workDir` or cluster queue names without explicit user request
- Duplicate `SPLIT_BAM` merge logic — `extract_beta.nf` mirrors `split_bam.nf`; refactor both if needed

## Related rules

| File | When it applies |
|------|-----------------|
| `.cursor/rules/episcore-project.mdc` | Always |
| `.cursor/rules/nextflow.mdc` | `*.nf` |
| `.cursor/rules/python.mdc` | `bin/*.py`, `scripts/**/*.py` |
| `.cursor/rules/hpc-config.mdc` | `conf/*`, `run_*.sh` |

## Offline analysis

`scripts/` holds batch Python and shell wrappers for parameter sweeps (thresholds, panel sizes, XO). These are **not** invoked by `main.nf`; treat as separate from pipeline changes unless the user asks to integrate them.
