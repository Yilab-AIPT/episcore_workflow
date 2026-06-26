process BETA_TO_EPISCORE {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(beta_value)
    path(reference_episcore_matrix)
    path(cpg_list)
    val(beta_depth_threshold)
    
    output:
    tuple val(meta), path("*_zscore.tsv"), emit: episcore
    
    script:
    def ref_matrix_arg = reference_episcore_matrix && reference_episcore_matrix.name != 'null' ? "--reference-episcore-matrix ${reference_episcore_matrix}" : ""
    def depth_arg = beta_depth_threshold != 'null' ? "--depth ${beta_depth_threshold}" : ""
    def args = task.ext.args ?: ''
    """
    beta_to_episcore.py \\
        --beta-value ${beta_value} \\
        --output-prefix ${meta.id} \\
        --ncpus ${task.cpus} \\
        --cpg-list ${cpg_list} \\
        ${ref_matrix_arg} \\
        ${depth_arg} \\
        ${args}
    """
}
