/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPLIT_BAM } from '../subworkflows/local/split_bam.nf'
include { MERGE_CLEAN_BAM } from '../subworkflows/local/merge_clean_bam.nf'
include { BAM_TO_PILEUP } from '../modules/local/bam_to_pileup/main.nf'
include { BAM_TO_PILEUP_SINGLE } from '../modules/local/bam_to_pileup_single/main.nf'
include { ESTIMATE_FF_HIGHER_PRECISION } from '../modules/local/estimate_ff_higher_precision/main.nf'
include { SUMMARIZE_SNP_FF } from '../modules/local/summarize_snp_ff/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SNP_EST_FF {
    take:
    ch_samplesheet // channel:
                   //   est_ff_from_bam                      -> [meta, clean_bam, deconv_res]
                   //   est_ff_from_bam_without_deconv_res   -> [meta, clean_bam]
                   //   est_ff_from_pileup                   -> [meta, pileup]

    main:
    // Resolve a per-sample merged pileup channel [meta, pileup] regardless of
    // the entry point.
    if (params.step == 'est_ff_from_bam') {
        // Split clean bam by deconv result, then build the SNP pileup.
        SPLIT_BAM(ch_samplesheet)
        SPLIT_BAM.out.splitted_bam
            .set { ch_splitted_bam }

        BAM_TO_PILEUP(
            ch_splitted_bam,
            file(params.snp_list)
        )
        BAM_TO_PILEUP.out.pileup
            .set { ch_pileup }
    } else if (params.step == 'est_ff_from_bam_without_deconv_res') {
        // Merge multi-BAM samples, then build pileup from the single clean BAM.
        MERGE_CLEAN_BAM(ch_samplesheet)
        MERGE_CLEAN_BAM.out.merged_bam
            .set { ch_merged_bam }

        BAM_TO_PILEUP_SINGLE(
            ch_merged_bam,
            file(params.snp_list)
        )
        BAM_TO_PILEUP_SINGLE.out.pileup
            .set { ch_pileup }
    } else {
        // est_ff_from_pileup: samplesheet already provides the merged pileup.
        ch_samplesheet
            .set { ch_pileup }
    }

    // Estimate fetal fraction with iterative high-precision grid search,
    // optionally restricting the pileup to the SNP panel (--known-sites).
    ESTIMATE_FF_HIGHER_PRECISION(
        ch_pileup,
        file(params.snp_list),
        params.snp_depth_threshold,
        params.snp_est_mode,
        params.ff_precision
    )
    ESTIMATE_FF_HIGHER_PRECISION.out.ff
        .set { ch_snp_ff }

    // Write a collated samplesheet of per-sample FF output paths for downstream use.
    ch_snp_ff
        .map { meta, ff ->
            def ff_path = "${params.outdir}/estimate_ff_higher_precision/${ff.name}"
            "${meta.id},${ff_path}\n"
        }
        .collectFile(
            name: 'snp_ff.csv',
            storeDir: "${params.outdir}/samplesheet",
            seed: "sample,ff\n"
        )

    // Summarize per-sample FF estimates into one table.
    SUMMARIZE_SNP_FF(
        ch_snp_ff.map { meta, ff -> ff }.collect()
    )

    emit:
    snp_pileup  = ch_pileup
    snp_ff      = ch_snp_ff
    ff_summary  = SUMMARIZE_SNP_FF.out.summary
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
