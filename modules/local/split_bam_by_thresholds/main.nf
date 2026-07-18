process SPLIT_BAM_BY_THRESHOLDS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam_file), path(deconv_res_file)
    val(thresholds)
    val(background_threshold)

    output:
    tuple val(meta), path("*__thr_*_target.bam"), emit: target_bams
    tuple val(meta), path("*__thr_*_background.bam"), emit: background_bams, optional: true

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bg_arg = (background_threshold != null && "${background_threshold}" != 'null') \
        ? "--background-threshold ${background_threshold}" : ""
    def bg_thr_fmt = (background_threshold != null && "${background_threshold}" != 'null') \
        ? (background_threshold as BigDecimal).stripTrailingZeros().toPlainString() : ""
    """
    set -euo pipefail
    export LC_ALL=C

    # Read the (large) deconv table once, emit one target read-name list per threshold
    # (and optional background list for raw_total depth filtering).
    split_reads_by_thresholds.py \\
        --input ${deconv_res_file} \\
        --thresholds "${thresholds}" \\
        --output-dir . \\
        ${bg_arg}

    # One target BAM per threshold (target = reads with prob_class_1 >= threshold).
    for t in \$(echo "${thresholds}" | tr ',' ' '); do
        samtools view -@ ${task.cpus} -b -N target_thr_\${t}.txt \\
            -o ${prefix}__thr_\${t}_target.bam ${bam_file}
        rm -f target_thr_\${t}.txt
    done

    # Optional background BAM for depth-filter threshold.
    if [ -n "${bg_thr_fmt}" ] && [ -f background_thr_${bg_thr_fmt}.txt ]; then
        samtools view -@ ${task.cpus} -b -N background_thr_${bg_thr_fmt}.txt \\
            -o ${prefix}__thr_${bg_thr_fmt}_background.bam ${bam_file}
        rm -f background_thr_${bg_thr_fmt}.txt
    fi
    """
}
