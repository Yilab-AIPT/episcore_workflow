//
// AIPT ref-40 fetal-fraction (ff_before_mq) from the clean BAM.
//
// Reuses the SNP pileup + higher-precision FF estimator. The clean BAM is fed
// as both the "target" and "background" pileup inputs (FF is estimated in cfDNA
// mode from read-count ratios, which are invariant to the resulting doubling),
// so ff_before_mq matches the whole-sample cfDNA estimate.
//

include { BAM_TO_PILEUP } from '../../modules/local/bam_to_pileup/main.nf'
include { ESTIMATE_FF_HIGHER_PRECISION } from '../../modules/local/estimate_ff_higher_precision/main.nf'

workflow AIPT_FF {
    take:
    ch_samplesheet      // channel: [meta, clean_bam, deconv_res]

    main:
    ch_samplesheet
        .map { meta, clean_bam, deconv_res -> [meta, clean_bam, clean_bam] }
        .set { ch_bam }

    BAM_TO_PILEUP(ch_bam, file(params.snp_list))

    ESTIMATE_FF_HIGHER_PRECISION(
        BAM_TO_PILEUP.out.pileup,
        file(params.snp_list),
        params.snp_depth_threshold,
        params.snp_est_mode,
        params.ff_precision
    )

    emit:
    ff = ESTIMATE_FF_HIGHER_PRECISION.out.ff
}
