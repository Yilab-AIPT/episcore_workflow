process EXTRACT_BETA_VALUE {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(target_bedgraph), path(background_bedgraph)
    path(cpg_list)

    output:
    tuple val(meta), path("*_beta_value.tsv.gz"), emit: beta_value

    script:
    """
    extract_beta_value.py \\
        --target-bedgraph ${target_bedgraph} \\
        --background-bedgraph ${background_bedgraph} \\
        --cpg-list ${cpg_list} \\
        --output-prefix ${meta.id} \\
        --ncpus ${task.cpus}
    """
}