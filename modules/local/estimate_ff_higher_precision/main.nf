process ESTIMATE_FF_HIGHER_PRECISION {
    tag "$meta.id"

    input:
    tuple val(meta), path(pileup_file)
    path(known_sites_tsv)
    val(snp_depth_threshold)
    val(snp_est_mode)
    val(ff_precision)

    output:
    tuple val(meta), path("*_ff.tsv"), emit: ff

    script:
    def depth_arg = snp_depth_threshold != 'null' ? "--min-raw-depth ${snp_depth_threshold}" : ""
    def mode_arg = snp_est_mode != 'null' ? "--mode-list ${snp_est_mode}" : ""
    def known_sites_arg = known_sites_tsv?.name != 'null' ? "--known-sites ${known_sites_tsv}" : ""
    def ff_precision_arg = ff_precision != 'null' ? "--ff-precision ${ff_precision}" : ""
    def args = task.ext.args ?: ''
    """
    estimate_ff_with_higher_precision.py \\
        --input-path ${pileup_file} \\
        --output-prefix ${meta.id} \\
        --ncpus ${task.cpus} \\
        ${depth_arg} \\
        ${mode_arg} \\
        ${known_sites_arg} \\
        ${ff_precision_arg} \\
        ${args}
    """
}
