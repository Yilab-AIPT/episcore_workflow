//
// AIPT ref-40 flexible per-chromosome read-count zscore (from deconv table).
//

include { CALC_ZSCORE_FLEXIBLE } from '../../modules/local/calc_zscore_flexible/main.nf'

workflow AIPT_ZSCORE {
    take:
    ch_samplesheet      // channel: [meta, clean_bam, deconv_res]

    main:
    ch_samplesheet
        .map { meta, clean_bam, deconv_res -> [meta, deconv_res] }
        .set { ch_deconv }

    def ref_matrix = file("${params.grid_search_result}/best_reference_matrix.tsv")
    def best_combo = file("${params.grid_search_result}/best_combo_zscore.csv")

    CALC_ZSCORE_FLEXIBLE(
        ch_deconv,
        best_combo,
        ref_matrix,
        file(params.cpg_recall_dir),
        params.zscore_mtcount
    )

    emit:
    zscore = CALC_ZSCORE_FLEXIBLE.out.zscore
}
