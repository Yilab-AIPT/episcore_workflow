/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateAndParseSamplesheet } from './lib/samplesheet_parser.nf'
include { validateAndParseGridSearchParameters } from './lib/grid_search_parameters_parser.nf'
include { validateAndParseSnpFFSamplesheet } from './lib/snp_est_ff_samplesheet_parser.nf'
include { validateAndParsePerturbedResSamplesheet } from './lib/perturbed_res_samplesheet_parser.nf'
include { validateAndParseAiptRef40Samplesheet } from './lib/aipt_ref_40_samplesheet_parser.nf'
include { NIPT  } from './workflows/nipt'
include { GRID_SEARCH } from './workflows/grid_search'
include { SNP_EST_FF } from './workflows/snp_est_ff'
include { PERTURBED_RES } from './workflows/perturbed_res'
include { AIPT_REF_40 } from './workflows/aipt_ref_40'

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

workflow EST_FF {
    take:
    ch_samplesheet

    main:
    SNP_EST_FF (
        ch_samplesheet
    )
}

workflow PERTURB {
    take:
    ch_samplesheet

    main:
    PERTURBED_RES (
        ch_samplesheet
    )
}

workflow REF_40 {
    take:
    ch_samplesheet

    main:
    AIPT_REF_40 (
        ch_samplesheet
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    if (params.step in ['split_bam', 'episcore']) {
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

    if (params.step in ['est_ff_from_bam', 'est_ff_from_pileup']) {
        // Validate and parse samplesheet
        ch_samplesheet = validateAndParseSnpFFSamplesheet(params.input, params.step)

        //
        // WORKFLOW: Run SNP-based fetal fraction estimation workflow
        //
        EST_FF (
            ch_samplesheet
        )
    }

    if (params.step in ['perturbed_res']) {
        // Validate and parse samplesheet
        ch_samplesheet = validateAndParsePerturbedResSamplesheet(params.input, params.step)

        //
        // WORKFLOW: Run methylation-perturbation workflow
        //
        PERTURB (
            ch_samplesheet
        )
    }

    if (params.step in ['aipt_ref_40']) {
        // Validate and parse samplesheet
        ch_samplesheet = validateAndParseAiptRef40Samplesheet(params.input)

        //
        // WORKFLOW: Run AIPT ref-40 episcore/zscore/ezscore scoring workflow
        //
        REF_40 (
            ch_samplesheet
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
