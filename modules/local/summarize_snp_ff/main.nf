process SUMMARIZE_SNP_FF {
    tag 'summary'

    publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true

    input:
    path ff_files, stageAs: 'ff/*'

    output:
    path 'ff_summary.tsv', emit: summary

    script:
    """
    summarize_snp_ff.py \\
        --input-dir ff \\
        --output ff_summary.tsv
    """
}
