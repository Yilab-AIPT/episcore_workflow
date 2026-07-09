process MERGE_SCORES {
    tag "$meta.id"

    input:
    tuple val(meta), path(episcore_tsv), path(zscore_tsv), path(ff_tsv)
    val(skip_ezscore)
    path(ezscore_matrix)

    output:
    tuple val(meta), path("*_scores.tsv"), emit: scores

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ez_arg = skip_ezscore ? '--skip-ezscore' : "--ezscore-matrix ${ezscore_matrix}"
    """
    merge_scores.py \\
        --sample ${meta.id} \\
        --episcore ${episcore_tsv} \\
        --zscore ${zscore_tsv} \\
        --ff ${ff_tsv} \\
        ${ez_arg} \\
        --output-prefix ${prefix}
    """
}
