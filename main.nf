/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateAndParseSamplesheet } from './lib/samplesheet_parser.nf'
include { validateAndParseGridSearchParameters } from './lib/grid_search_parameters_parser.nf'
include { NIPT  } from './workflows/nipt'
include { GRID_SEARCH } from './workflows/grid_search'

workflow MAIN {
    take:
    ch_samplesheet

    main:
    NIPT (
        ch_samplesheet
    )
}

workflow SUB {
    take:
    ch_samplesheet

    main:
    GRID_SEARCH (
        ch_samplesheet
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    if (params.step in ['preprocess', 'deconv', 'split_bam', 'beta_zscore', 'rc_zscore']) {
        // Validate and parse samplesheet
        ch_samplesheet = validateAndParseSamplesheet(params.input, params.step)

        //
        // WORKFLOW: Run main workflow
        //
        MAIN (
            ch_samplesheet
        )
    }

    if (params.step in ['grid_search']) {
        // Validate and parse samplesheet
        ch_samplesheet = validateAndParseGridSearchParameters(params.input)

        //
        // WORKFLOW: Run grid search workflow
        //
        SUB (
            ch_samplesheet
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
