/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPLIT_BAM } from '../subworkflows/local/split_bam.nf'
include { CALC_EPISCORE } from '../subworkflows/local/calc_episcore.nf'
include { ESTIMATE_FF } from '../subworkflows/local/estimate_ff.nf'
include { REPORT } from '../subworkflows/local/report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NIPT {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // Split clean bam by deconv result (preprocess/deconv done upstream)
    if (params.step == 'split_bam') {
        SPLIT_BAM(ch_samplesheet)
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

    // Calculate episcore and SNP-est-ff
    if (params.step in ['split_bam', 'episcore']) {
        CALC_EPISCORE(ch_splitted_bam)
        CALC_EPISCORE.out.episcore
            .set { ch_episcore }
        CALC_EPISCORE.out.beta_value
            .set { ch_beta_value }

        ESTIMATE_FF(ch_splitted_bam)
        ESTIMATE_FF.out.snp_pileup
            .set { ch_snp_pileup }
        ESTIMATE_FF.out.snp_ff
            .set { ch_snp_ff }
    } else {
        ch_episcore = channel.empty()
        ch_beta_value = channel.empty()
        ch_snp_pileup = channel.empty()
        ch_snp_ff = channel.empty()
    }

    // Generate final report
    if (params.step in ['split_bam', 'episcore']) {
        REPORT(
            ch_episcore,
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
