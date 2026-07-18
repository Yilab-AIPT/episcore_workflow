process CALC_EPISCORE_FLEXIBLE {
    tag "$meta.id"

    input:
    tuple val(meta), path(bedgraphs), path(depth_filtered_cpgs)
    path(best_combo_episcore)
    path(reference_matrix)
    path(cpg_recall_dir)

    output:
    tuple val(meta), path("*_episcore.tsv"), emit: episcore

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedgraph_args = bedgraphs.collect { "--bedgraph ${it}" }.join(' ')
    """
    calc_episcore_flexible.py \\
        ${bedgraph_args} \\
        --best-combo-episcore ${best_combo_episcore} \\
        --reference-matrix ${reference_matrix} \\
        --cpg-recall-dir ${cpg_recall_dir} \\
        --depth-filtered-cpgs ${depth_filtered_cpgs} \\
        --output-prefix ${prefix}
    """
}
