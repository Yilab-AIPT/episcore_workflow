process SPLIT_BAM_BY_THRESHOLDS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam_file), path(deconv_res_file)
    val(thresholds)

    output:
    tuple val(meta), path("*__thr_*_target.bam"), emit: target_bams

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C

    # Read the (large) deconv table once, emit one target read-name list per threshold.
    split_reads_by_thresholds.py \\
        --input ${deconv_res_file} \\
        --thresholds "${thresholds}" \\
        --output-dir .

    # One target BAM per threshold (target = reads with prob_class_1 >= threshold).
    for t in \$(echo "${thresholds}" | tr ',' ' '); do
        samtools view -@ ${task.cpus} -b -N target_thr_\${t}.txt \\
            -o ${prefix}__thr_\${t}_target.bam ${bam_file}
        rm -f target_thr_\${t}.txt
    done
    """
}
