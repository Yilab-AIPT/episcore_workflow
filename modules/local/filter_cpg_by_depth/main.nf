process FILTER_CPG_BY_DEPTH {
    tag "$meta.id"

    input:
    tuple val(meta), path(target_bedgraph), path(background_bedgraph)
    val(depth)

    output:
    tuple val(meta), path("*_depth_filtered_cpgs.tsv"), emit: filtered_cpgs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_cpg_by_depth.py \\
        --target-bedgraph ${target_bedgraph} \\
        --background-bedgraph ${background_bedgraph} \\
        --depth ${depth} \\
        --output-prefix ${prefix}
    """
}
