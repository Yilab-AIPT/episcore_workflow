# Episcore — Agent Guide

## What this project is

**Episcore** is a Nextflow pipeline for **NIPT trisomy detection** using methylation **episcore** (beta → chromosome-level z-scores) plus **SNP-based fetal fraction (FF)** estimation. It sits in the AIPT 2.0 workflow at `/lustre1/cqyi/AIPT_2.0/workflow/episcore`.

Upstream work (outside this repo) produces **clean BAMs** and **deconvolution read-probability tables**; this pipeline splits reads by deconv threshold, extracts methylation (MethylDackel), computes beta/z-scores against a reference matrix, estimates FF from SNPs, and emits per-sample reports.

## Entry points

| `params.step` | Workflow | Samplesheet columns |
|---------------|----------|---------------------|
| `split_bam` | `NIPT` → `SPLIT_BAM` | `sample`, `clean_bam`, `deconv_res` |
| `beta_zscore` | `NIPT` → `CALC_BETA_ZSCORE`, `ESTIMATE_FF`, `REPORT` | `sample`, `target_bam`, `background_bam` |
| `grid_search` | `GRID_SEARCH` → `EXTRACT_BETA` | `sample`, `clean_bam`, `deconv_res` |
| `est_ff_from_bam` | `SNP_EST_FF` → `SPLIT_BAM`, `BAM_TO_PILEUP`, `ESTIMATE_FF_HIGHER_PRECISION` | `sample`, `clean_bam`, `deconv_res` |
| `est_ff_from_pileup` | `SNP_EST_FF` → `ESTIMATE_FF_HIGHER_PRECISION` | `sample`, `pileup` |

`main.nf` routes by `params.step`: NIPT steps use `validateAndParseSamplesheet`; grid search uses `validateAndParseGridSearchParameters`; SNP FF steps (`est_ff_from_bam` / `est_ff_from_pileup`) use `validateAndParseSnpFFSamplesheet` (`lib/`).

The **SNP FF** workflow (`workflows/snp_est_ff.nf`, entry `EST_FF`) is a standalone path focused only on SNP-based fetal fraction. It reuses `SPLIT_BAM` + `BAM_TO_PILEUP` to build a pileup (`est_ff_from_bam`) or consumes a pre-computed pileup (`est_ff_from_pileup`), then runs `ESTIMATE_FF_HIGHER_PRECISION` (iterative range-narrowing grid search, `bin/estimate_ff_with_higher_precision.py`). It does **not** go through `CALC_BETA_ZSCORE` / `REPORT`.

## Directory map

```
main.nf                 # Top-level router (MAIN / SUB / EST_FF workflows)
nextflow.config         # params, profiles, workDir, reports
workflows/              # nipt.nf, grid_search.nf, snp_est_ff.nf
subworkflows/local/     # split_bam, calc_beta_zscore, estimate_ff, report, extract_beta
modules/local/          # Project-specific processes
modules/nf-core/        # samtools, picard, methyldackel (vendored)
lib/                    # Groovy samplesheet / grid-search / snp-ff parsers
bin/                    # Python CLI scripts invoked by processes
conf/                   # Profile-specific params + executor/container config
assets/                 # Reference FASTA, CpG/SNP lists, reference z-score matrices (not in git)
containers/             # Singularity images (gitignored)
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
  --step split_bam   # or beta_zscore | est_ff_from_bam | est_ff_from_pileup
```

- **workDir** (fixed in `nextflow.config`): `/lustre1/cqyi/AIPT_2.0/tmp/episcore_workflow_tmp_dir`
- **Profiles**: `early`, `early_240k`, `middle`, `filter_size`, `at_analyze`, `at_ref`, `grid_search`, `test`, plus executors `alioth_slurm`, `alioth_local`, `dev`, `singularity`
- Panel variants differ mainly in `cpg_list`, `snp_list`, `reference_beta_zscore_matrix` under `conf/*.config`

## Channel conventions

- **meta**: `[id: sample]` from samplesheet `sample` column
- **split_bam input**: `[meta, clean_bam, deconv_res]`
- **beta_zscore input**: `[meta, target_bam, background_bam]`
- **est_ff_from_bam input**: `[meta, clean_bam, deconv_res]`
- **est_ff_from_pileup input**: `[meta, pileup]`
- Processes use `tag "$meta.id"`; outputs often prefixed with `${meta.id}`

## Python ↔ Nextflow contract

- Scripts in `bin/` are called by name from process `script:` blocks (on `PATH` via Nextflow `bin/`).
- CLI: **click** + **rich**; tabular IO: **pandas** / **polars** (see `filter_deconv_res.py`, `merge_deconv_res.py`).
- `estimate_ff.py` imports `FFEstimator` from `bin/FFEstimator.py` (same directory in container).
- `estimate_ff_with_higher_precision.py` (SNP FF workflow) reuses `FFEstimator` + `load_and_validate_data`/`parse_*` from `estimate_ff.py`; iterative range-narrowing search via `--ff-precision`. `--mode-list` ← `params.snp_est_mode`, `--min-raw-depth` ← `params.snp_depth_threshold`, optional `--known-sites` ← `params.snp_list` (shared with `bam_to_pileup.py`, filters the pileup to panel sites before estimation). It is the offline `scripts/ff_decimal/` variant adapted to a single pileup per process.

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
