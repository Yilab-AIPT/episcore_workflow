process SNP_TO_FF {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(pileup_file)
    val(snp_depth_threshold)
    val(snp_est_mode)
    
    output:
    tuple val(meta), path("*_ff.tsv"), emit: ff
    
    script:
    def depth_arg = snp_depth_threshold != 'null' ? "--min-raw-depth ${snp_depth_threshold}" : ""
    def mode_arg = snp_est_mode != 'null' ? "--mode-list ${snp_est_mode}" : ""
    def args = task.ext.args ?: ''
    """
    estimate_ff.py \\
        --input-path ${pileup_file} \\
        --output-prefix ${meta.id} \\
        --ncpus ${task.cpus} \\
        ${depth_arg} \\
        ${mode_arg} \\
        ${args}
    """
}
