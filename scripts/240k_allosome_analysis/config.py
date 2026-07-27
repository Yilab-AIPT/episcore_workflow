"""Shared paths for 240k allosome analysis (20260720 cohort)."""

from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]  # workflow/episcore

INPUT_DIR = Path(
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260720-240k_early_allosomes_samples"
)
OUTPUT_DIR = Path(
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/20260720-240k_early_allosomes_samples"
)

# Prior 240k cohorts (scatter background)
OLD_SAMPLESHEET_20260416 = Path(
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260416-240k.csv"
)
OLD_SAMPLESHEET_20260507 = Path(
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/20260507-240k_XO_samples.csv"
)
OLD_FF_20260416 = Path(
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/"
    "20260416-early_samples-240krecall0/collect_reports/summary_report.tsv"
)
OLD_FF_20260507 = Path(
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/"
    "20260507-240k_XO_samples-240krecall0/collect_reports/summary_report.tsv"
)
OLD_BETA_20260416 = Path(
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/"
    "20260416-early_samples-240krecall0/extract_beta_value"
)
OLD_BETA_20260507 = Path(
    "/lustre1/cqyi/AIPT_2.0/results/episcore_output/"
    "20260507-240k_XO_samples-240krecall0/extract_beta_value"
)

# Recall CpG lists (240k panel)
CPG_RECALL_DIR = Path(
    "/lustre1/cqyi/AIPT_2.0/data/meta/episcore/"
    "20260525-grid_search_240k_panel_240k_model/recall_list_240k"
)
FULL_240K_CPG_LIST = PROJECT_DIR / "assets" / "cpgs_in_240k_probes.txt"

# Score cutoffs
EPISCORE_THRESHOLD = 0.5  # BAM-split / prob_class_1 for methylation episcore
ZSCORE_CUTOFF = 0.85  # prob_class_1 cutoff for read-count zscore
ZSCORE_MTCOUNT = 1.0

# Derived meta / result locations under INPUT_DIR / OUTPUT_DIR
COHORT_LABELS = INPUT_DIR / "cohort_labels.csv"
MQRES = INPUT_DIR / "mqres.csv"
META = INPUT_DIR / "meta.csv"
EPISCORE_SAMPLES_META = INPUT_DIR / "episcore_samples_meta.csv"
ZSCORE_SAMPLES_META = INPUT_DIR / "zscore_samples_meta.csv"
NF_SAMPLESHEET = INPUT_DIR / "samplesheet_nf.csv"

CHRY_FF_TSV = OUTPUT_DIR / "tables" / "chry_ff.tsv"
EPISCORE_RECALL_DIR = OUTPUT_DIR / "episcore_recall"
ZSCORE_RECALL_DIR = OUTPUT_DIR / "zscore_recall"
EPISCORE_COLLECTED = OUTPUT_DIR / "tables" / "chrX_episcore_vs_recall.tsv"
ZSCORE_COLLECTED = OUTPUT_DIR / "tables" / "chrX_zscore_vs_recall.tsv"
PLOTS_DIR = OUTPUT_DIR / "plots"

SINGULARITY_IMAGE = PROJECT_DIR / "containers" / "common_tools.sif"
MAIN_NF = PROJECT_DIR / "main.nf"
