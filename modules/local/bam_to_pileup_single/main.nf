process BAM_TO_PILEUP_SINGLE {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)
    path(known_sites_tsv)

    output:
    tuple val(meta), path("*_pileup.tsv.gz"), emit: pileup

    script:
    """
    samtools index ${bam}
    bam_to_pileup.py \\
        --input-bam ${bam} \\
        --known-sites ${known_sites_tsv} \\
        --output ${meta.id} \\
        --ncpus ${task.cpus}
    """
}
