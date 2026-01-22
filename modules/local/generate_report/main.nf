process GENERATE_REPORT {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(zscore), path(beta_value), path(snp_pileup), path(snp_ff)
    
    output:
    tuple val(meta), path("*_report.tsv"), emit: report
    
    script:
    def meta_arg = params.meta ? "--meta ${params.meta}" : ""
    """
    generate_report.py \\
        --sample-id ${meta.id} \\
        --zscore ${zscore} \\
        --beta-value ${beta_value} \\
        --snp-pileup ${snp_pileup} \\
        --snp-ff ${snp_ff} \\
        --output-prefix ${meta.id} \\
        ${meta_arg}
    """
}
