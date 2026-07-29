//
// Collapse a multi-row samplesheet to one clean BAM per sample.
//
// When a sample has multiple BAMs, merge them (SAMTOOLS_MERGE) and remove
// duplicates (PICARD_MARKDUPLICATES). Single-BAM samples pass through untouched.
//

include { SAMTOOLS_MERGE        } from '../../modules/nf-core/samtools/merge/main.nf'
include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main.nf'

workflow MERGE_CLEAN_BAM {
    take:
    ch_samplesheet // channel: [meta, clean_bam]

    main:
    ch_samplesheet
        .map { meta, clean_bam ->
            [meta.id.toString(), meta, clean_bam]
        }
        .groupTuple(by: 0)
        .map { groupKey, meta, clean_bam ->
            def new_meta = meta.first()
            def clean_bam_list = clean_bam.unique { it.toString() } as List
            [new_meta, clean_bam_list, clean_bam_list.size()]
        }
        .branch {
            multi:  it[2] > 1
            single: it[2] == 1
        }
        .set { ch_bam_branched }

    SAMTOOLS_MERGE(
        ch_bam_branched.multi.map { meta, clean_bam_list, count -> [meta, clean_bam_list] },
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)],
        [[:], []]
    )

    PICARD_MARKDUPLICATES(
        SAMTOOLS_MERGE.out.bam,
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)]
    )

    PICARD_MARKDUPLICATES.out.bam
        .mix(
            ch_bam_branched.single.map { meta, clean_bam_list, count ->
                [meta, clean_bam_list.first()]
            }
        )
        .set { ch_merged_bam }

    emit:
    merged_bam = ch_merged_bam // channel: [meta, clean_bam]
}
