/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPLIT_BAM } from '../subworkflows/local/split_bam.nf'
include { CALC_BETA_ZSCORE } from '../subworkflows/local/calc_beta_zscore.nf'
include { ESTIMATE_FF } from '../subworkflows/local/estimate_ff.nf'
include { REPORT } from '../subworkflows/local/report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NIPT {
    take:
    ch_samplesheet // channel: samplesheet read in from --input, with columns ['sample', 'tissue', 'bam']

    main:
    // Preprocess raw seqeuncing data, from fastq to clean bam
    if (params.step == 'preprocess') {
        // Return ch_clean_bam
    } else {
        ch_clean_bam = channel.empty().mix(ch_samplesheet)
    }

    // Deconvolute the clean bam to get read probability using MethylQueen
    if (params.step in ['preprocess', 'deconv']) {
        // Return ch_deconv_res
    } else {
        ch_deconv_res = channel.empty().mix(ch_samplesheet)
    }

    // Split clean bam by deconv result
    if (params.step in ['preprocess', 'deconv', 'split_bam']) {
        // Split clean bam by deconv result
        SPLIT_BAM(ch_deconv_res)
        SPLIT_BAM.out.splitted_bam
            .set { ch_splitted_bam }

        ch_splitted_bam
            .map {meta, target_bam, background_bam ->
                def target_bam_path = "${params.outdir}/split_bam_by_deconv_res/${target_bam.name}"
                def background_bam_path = "${params.outdir}/split_bam_by_deconv_res/${background_bam.name}"
                "${meta.id},${target_bam_path},${background_bam_path}\n"
            }
            .collectFile(
                name: 'splitted_bam.csv',
                storeDir: "${params.outdir}/samplesheet",
                seed: "sample,target_bam,background_bam\n"
            )
    } else {
        ch_splitted_bam = channel.empty().mix(ch_samplesheet)
    }

    // Calculate beta-zscore and SNP-est-ff
    if (params.step in ['preprocess', 'deconv', 'split_bam', 'beta_zscore']) {
        // Calculate beta-zscore
        CALC_BETA_ZSCORE(ch_splitted_bam)
        CALC_BETA_ZSCORE.out.zscore
            .set { ch_zscore }
        CALC_BETA_ZSCORE.out.beta_value
            .set { ch_beta_value }

        // Estimate fetal fraction
        ESTIMATE_FF(ch_splitted_bam)
        ESTIMATE_FF.out.snp_pileup
            .set { ch_snp_pileup }
        ESTIMATE_FF.out.snp_ff
            .set { ch_snp_ff }
    } else {
        ch_zscore = channel.empty()
        ch_beta_value = channel.empty()
        ch_snp_pileup = channel.empty()
        ch_snp_ff = channel.empty()
    }

    // Calculate rc-zscore
    if (params.step in ['preprocess', 'deconv', 'rc_zscore']) {
        // Return ch_rc_zscore
    } else {
        ch_rc_zscore = channel.empty()
    }

    // Generate final report
    // For early version: only focusing on beta_zscore results
    // Interfaces for ch_clean_bam/ch_deconv_res/ch_splitted_bam/ch_rc_zscore are left for future development
    if (params.step in ['preprocess', 'deconv', 'split_bam', 'beta_zscore']) {
        REPORT(
            ch_zscore,
            ch_beta_value,
            ch_snp_pileup,
            ch_snp_ff
        )
        REPORT.out.report
            .set { ch_report }
    } else {
        ch_report = channel.empty()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
