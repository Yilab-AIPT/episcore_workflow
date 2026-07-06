process CALC_ZSCORE_FLEXIBLE {
    tag "$meta.id"

    input:
    tuple val(meta), path(deconv_res)
    path(best_combo_zscore)
    path(reference_matrix)
    path(cpg_recall_dir)
    val(zscore_mtcount)

    output:
    tuple val(meta), path("*_zscore.tsv"), emit: zscore

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mtcount_arg = (zscore_mtcount != null && "${zscore_mtcount}" != 'null') ? "--mtcount ${zscore_mtcount}" : ""
    """
    calc_zscore_flexible.py \\
        --deconv-res ${deconv_res} \\
        --best-combo-zscore ${best_combo_zscore} \\
        --reference-matrix ${reference_matrix} \\
        --cpg-recall-dir ${cpg_recall_dir} \\
        --output-prefix ${prefix} \\
        ${mtcount_arg}
    """
}
