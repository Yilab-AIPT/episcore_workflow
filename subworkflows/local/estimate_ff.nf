//
// Estimate fetal fraction using bam
//

include { BAM_TO_PILEUP } from '../../modules/local/bam_to_pileup/main.nf'
include { SNP_TO_FF } from '../../modules/local/snp_to_ff/main.nf'

workflow ESTIMATE_FF {
    take:
    ch_samplesheet // channel: samplesheet with columns ['sample', 'target_bam', 'background_bam']

    main:
    BAM_TO_PILEUP(
        ch_samplesheet,
        file(params.snp_list)
    )
    BAM_TO_PILEUP.out.pileup
        .set { ch_pileup }

    SNP_TO_FF(
        ch_pileup,
        params.snp_depth_threshold,
        params.snp_est_mode
    )
    SNP_TO_FF.out.ff
        .set { ch_ff }

    emit:
    snp_pileup = ch_pileup
    snp_ff = ch_ff
}