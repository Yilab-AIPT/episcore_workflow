process SPLIT_BAM_BY_DECONV_RES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bam_file), path(deconv_res_file)
    val(threshold)

    output:
    tuple val(meta), path("*_target.bam"), path("*_background.bam"), emit: splitted_bam
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    set -euo pipefail
    export LC_ALL=C

    # Extract read names based on probability threshold using Python script
    # Automatically detects file format (txt/tsv or parquet)
    # Memory-efficient: chunk-based processing with automatic deduplication
    filter_deconv_res.py \\
        --input ${deconv_res_file} \\
        --threshold ${threshold} \\
        --output-dir . \\
        ${args}
    
    # Use samtools view -N to extract reads (fastest method for BAM splitting)
    samtools view -@ ${task.cpus} -b -N target_reads.txt     -o ${prefix}_target.bam ${bam_file}
    samtools view -@ ${task.cpus} -b -N background_reads.txt -o ${prefix}_background.bam ${bam_file}
    samtools view -@ ${task.cpus} -b -N classified_reads.txt -U ${prefix}_unclassified.bam ${bam_file} > /dev/null
    """
}
