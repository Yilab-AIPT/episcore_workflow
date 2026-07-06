/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUTS } from '../subworkflows/local/prepare_inputs.nf'
include { AIPT_EPISCORE  } from '../subworkflows/local/aipt_episcore.nf'
include { AIPT_ZSCORE    } from '../subworkflows/local/aipt_zscore.nf'
include { AIPT_FF        } from '../subworkflows/local/aipt_ff.nf'
include { MERGE_SCORES   } from '../modules/local/merge_scores/main.nf'
include { PLOT_SCORES    } from '../modules/local/plot_scores/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Distinct episcore probability thresholds from best_combo_episcore.csv,
// normalised to the canonical (%g-equivalent) string so that the read-name
// files emitted by split_reads_by_thresholds.py match the BAM names.
def parseEpiscoreThresholds(grid_search_result) {
    def combo_file = file("${grid_search_result}/best_combo_episcore.csv")
    if (!combo_file.exists()) {
        error "Missing best_combo_episcore.csv under --grid_search_result: ${combo_file}"
    }
    def thresholds = combo_file.readLines()
        .drop(1)
        .findAll { it?.trim() }
        .collect { it.split(',')[1].trim() }
        .collect { (it as BigDecimal).stripTrailingZeros().toPlainString() }
        .unique()
    return thresholds.join(',')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AIPT_REF_40 {
    take:
    ch_samplesheet  // channel: [meta, clean_bam, deconv_res]

    main:
    def ep_thresholds = parseEpiscoreThresholds(params.grid_search_result)

    // Collapse multi-row samples to one BAM (merged + deduped) and one merged
    // deconv_res per sample before any scoring.
    PREPARE_INPUTS(ch_samplesheet)
    ch_prepared = PREPARE_INPUTS.out.samplesheet

    // Per-chromosome scores
    AIPT_EPISCORE(ch_prepared, ep_thresholds)
    AIPT_ZSCORE(ch_prepared)
    AIPT_FF(ch_prepared)

    // Combine episcore + zscore + ff -> {sample}_scores.tsv (+ ezscore)
    AIPT_EPISCORE.out.episcore
        .join(AIPT_ZSCORE.out.zscore, by: 0)
        .join(AIPT_FF.out.ff, by: 0)
        .set { ch_merge_input }

    MERGE_SCORES(
        ch_merge_input,
        file("${params.grid_search_result}/best_ezscore_ref_20_matrix.tsv")
    )

    // Visualise score distribution per sample
    PLOT_SCORES(
        MERGE_SCORES.out.scores,
        file("${params.grid_search_result}/best_sample_scores_recalc_ezscore.tsv"),
        params.score_cutoff
    )

    emit:
    scores = MERGE_SCORES.out.scores
    plot   = PLOT_SCORES.out.plot
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
