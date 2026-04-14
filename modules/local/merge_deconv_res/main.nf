process MERGE_DECONV_RES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(deconv_res_files, stageAs: "input_?")
    val(ncpgs)
    
    output:
    tuple val(meta), path("*_merged_deconv_res.parquet"), emit: merged_deconv_res
    
    script:
    """
    merge_deconv_res.py \\
        --inputs "\$(ls ${deconv_res_files.join(' ')} | tr '\\n' ' ')" \\
        --output ${meta.id}_merged_deconv_res.parquet \\
        --ncpgs ${ncpgs}
    """
}
