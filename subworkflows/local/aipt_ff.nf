//
// AIPT ref-40 SNP-based fetal fraction (ff_before_mq + ff_after_mq).
//
// Splits the clean BAM by deconv probabilities, builds target/background pileups,
// then runs the higher-precision FF estimator. ff_before_mq is the cfDNA baseline;
// ff_after_mq is the enriched estimate from model-selected reads.
//
// Uses SPLIT_BAM_BY_DECONV_RES directly (not full SPLIT_BAM): PREPARE_INPUTS
// upstream already merges multi-row samples to one BAM + one deconv_res per sample.
//

include { SPLIT_BAM_BY_DECONV_RES } from '../../modules/local/split_bam_by_deconv_res/main.nf'
include { BAM_TO_PILEUP } from '../../modules/local/bam_to_pileup/main.nf'
include { ESTIMATE_FF_HIGHER_PRECISION } from '../../modules/local/estimate_ff_higher_precision/main.nf'

workflow AIPT_FF {
    take:
    ch_samplesheet      // channel: [meta, clean_bam, deconv_res]  (collapsed by PREPARE_INPUTS)

    main:
    SPLIT_BAM_BY_DECONV_RES(
        ch_samplesheet,
        params.threshold
    )

    BAM_TO_PILEUP(
        SPLIT_BAM_BY_DECONV_RES.out.splitted_bam,
        file(params.snp_list)
    )

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
