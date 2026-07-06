process PLOT_SCORES {
    tag "$meta.id"

    input:
    tuple val(meta), path(scores_tsv)
    path(precomputed_scores)
    val(score_cutoff)

    output:
    tuple val(meta), path("*_scores.pdf"), emit: plot

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cutoff_arg = (score_cutoff != null && "${score_cutoff}" != 'null') ? "--threshold ${score_cutoff}" : ""
    """
    plot_scores.py \\
        --scores-tsv ${scores_tsv} \\
        --precomputed-tsv ${precomputed_scores} \\
        --output-prefix ${prefix} \\
        ${cutoff_arg}
    """
}
