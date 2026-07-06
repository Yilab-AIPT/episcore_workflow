process CALC_EPISCORE_FLEXIBLE {
    tag "$meta.id"

    input:
    tuple val(meta), path(bedgraphs)
    path(best_combo_episcore)
    path(reference_matrix)
    path(cpg_recall_dir)
    val(beta_depth_threshold)

    output:
    tuple val(meta), path("*_episcore.tsv"), emit: episcore

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def depth_arg = (beta_depth_threshold != null && "${beta_depth_threshold}" != 'null') ? "--depth ${beta_depth_threshold}" : ""
    def bedgraph_args = bedgraphs.collect { "--bedgraph ${it}" }.join(' ')
    """
    calc_episcore_flexible.py \\
        ${bedgraph_args} \\
        --best-combo-episcore ${best_combo_episcore} \\
        --reference-matrix ${reference_matrix} \\
        --cpg-recall-dir ${cpg_recall_dir} \\
        --output-prefix ${prefix} \\
        ${depth_arg}
    """
}
