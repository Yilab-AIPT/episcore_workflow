process MERGE_SCORES {
    tag "$meta.id"

    input:
    tuple val(meta), path(episcore_tsv), path(zscore_tsv), path(ff_tsv)
    val(skip_ezscore)
    path(ezscore_matrix)
    val(skip_ff)

    output:
    tuple val(meta), path("*_scores.tsv"), emit: scores

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ez_arg = skip_ezscore ? '--skip-ezscore' : "--ezscore-matrix ${ezscore_matrix}"
    def ff_arg = skip_ff ? '--skip-ff' : "--ff ${ff_tsv}"
    """
    merge_scores.py \\
        --sample ${meta.id} \\
        --episcore ${episcore_tsv} \\
        --zscore ${zscore_tsv} \\
        ${ff_arg} \\
        ${ez_arg} \\
        --output-prefix ${prefix}
    """
}
