process BETA_TO_ZSCORE {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(beta_value)
    path(reference_beta_zscore_matrix)
    path(cpg_list)
    val(beta_depth_threshold)
    
    output:
    tuple val(meta), path("*_zscore.tsv"), emit: zscore
    
    script:
    def ref_matrix_arg = reference_beta_zscore_matrix && reference_beta_zscore_matrix.name != 'null' ? "--reference-beta-zscore-matrix ${reference_beta_zscore_matrix}" : ""
    def depth_arg = beta_depth_threshold != 'null' ? "--depth ${beta_depth_threshold}" : ""
    def args = task.ext.args ?: ''
    """
    beta_to_zscore.py \\
        --beta-value ${beta_value} \\
        --output-prefix ${meta.id} \\
        --ncpus ${task.cpus} \\
        --cpg-list ${cpg_list} \\
        ${ref_matrix_arg} \\
        ${depth_arg} \\
        ${args}
    """
}
