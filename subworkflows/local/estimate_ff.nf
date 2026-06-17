//
// Estimate fetal fraction using bam
//

include { BAM_TO_PILEUP } from '../../modules/local/bam_to_pileup/main.nf'
include { SNP_TO_FF } from '../../modules/local/snp_to_ff/main.nf'
include { ESTIMATE_FF_HIGHER_PRECISION } from '../../modules/local/estimate_ff_higher_precision/main.nf'

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

    // When ff_precision is set, use iterative range-narrowing grid search
    // (estimate_ff_with_higher_precision.py); otherwise use the standard
    // fixed-step search (estimate_ff.py).
    if (params.ff_precision) {
        ESTIMATE_FF_HIGHER_PRECISION(
            ch_pileup,
            file(params.snp_list),
            params.snp_depth_threshold,
            params.snp_est_mode,
            params.ff_precision
        )
        ESTIMATE_FF_HIGHER_PRECISION.out.ff
            .set { ch_ff }
    } else {
        SNP_TO_FF(
            ch_pileup,
            params.snp_depth_threshold,
            params.snp_est_mode
        )
        SNP_TO_FF.out.ff
            .set { ch_ff }
    }

    emit:
    snp_pileup = ch_pileup
    snp_ff = ch_ff
}