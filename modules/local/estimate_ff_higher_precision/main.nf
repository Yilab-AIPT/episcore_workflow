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
    // Groovy stringifies small floats as scientific notation (e.g. 1.0E-4),
    // but estimate_ff_with_higher_precision.py requires plain decimals (0.0001).
    def ff_precision_plain = ff_precision != 'null'
        ? new java.math.BigDecimal(ff_precision.toString()).stripTrailingZeros().toPlainString()
        : null
    def ff_precision_arg = ff_precision_plain ? "--ff-precision ${ff_precision_plain}" : ""
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
