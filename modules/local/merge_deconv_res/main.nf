process MERGE_DECONV_RES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(deconv_res_files, stageAs: "input_?.txt")
    val(ncpgs)
    
    output:
    tuple val(meta), path("merged_deconv_res.txt"), emit: merged_deconv_res
    
    script:
    """
    merge_deconv_res.py \\
        --inputs "\$(ls ${deconv_res_files.join(' ')} | tr '\\n' ' ')" \\
        --output merged_deconv_res.txt \\
        --ncpgs ${ncpgs}
    """
}
