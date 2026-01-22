process BAM_TO_PILEUP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(target_bam, stageAs: "target.bam"), path(background_bam, stageAs: "background.bam")
    path(known_sites_tsv)
    
    output:
    tuple val(meta), path("*_pileup.tsv.gz"), emit: pileup
    
    script:
    """
    samtools index target.bam
    bam_to_pileup.py \\
        --input-bam target.bam \\
        --known-sites ${known_sites_tsv} \\
        --output target \\
        --ncpus ${task.cpus}

    samtools index background.bam
    bam_to_pileup.py \\
        --input-bam background.bam \\
        --known-sites ${known_sites_tsv} \\
        --output background \\
        --ncpus ${task.cpus}

    # Merge target and background pileup
    merge_pileup.py \\
        --target_pileup target_pileup.tsv.gz \\
        --background_pileup background_pileup.tsv.gz \\
        --output-prefix ${meta.id}_pileup

    # Remove intermediate files
    rm -f target_pileup.tsv.gz background_pileup.tsv.gz
    """
}
