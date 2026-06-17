process REPLACE_DECONV_PROB {
    tag "$meta.id"

    input:
    tuple val(meta), path(original_res), path(perturbed_res)

    output:
    tuple val(meta), path("*_perturbed_deconv_res.parquet"), emit: perturbed_deconv_res

    script:
    """
    replace_deconv_prob.py \\
        --original ${original_res} \\
        --perturbed ${perturbed_res} \\
        --output ${meta.id}_perturbed_deconv_res.parquet
    """
}
