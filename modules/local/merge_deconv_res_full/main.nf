process MERGE_DECONV_RES_FULL {
    tag "$meta.id"

    input:
    tuple val(meta), path(deconv_res_files, stageAs: "input_?")

    output:
    tuple val(meta), path("*_merged_deconv_res.parquet"), emit: merged_deconv_res

    script:
    """
    merge_deconv_res_full.py \\
        --inputs "\$(ls ${deconv_res_files instanceof List ? deconv_res_files.join(' ') : deconv_res_files} | tr '\\n' ' ')" \\
        --output ${meta.id}_merged_deconv_res.parquet
    """
}
