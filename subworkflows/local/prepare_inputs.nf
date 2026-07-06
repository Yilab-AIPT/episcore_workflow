//
// Collapse the AIPT ref-40 samplesheet to a clean 1 sample : 1 BAM : 1 deconv_res
// relationship before any downstream scoring.
//
// A sample may appear on multiple rows (e.g. several lanes/runs, or several
// deconvolution results for the same clean BAM). Following the NIPT SPLIT_BAM
// pattern we group by sample, then:
//   * merge the distinct clean BAMs (SAMTOOLS_MERGE) and mark/remove duplicates
//     (PICARD_MARKDUPLICATES) when there is more than one BAM;
//   * merge the distinct deconv_res files (MERGE_DECONV_RES_FULL) when there is
//     more than one, preserving every column needed by the episcore split and
//     the read-count z-score.
// Single-BAM / single-deconv samples are passed through untouched.
//

include { SAMTOOLS_MERGE          } from '../../modules/nf-core/samtools/merge/main.nf'
include { PICARD_MARKDUPLICATES   } from '../../modules/nf-core/picard/markduplicates/main.nf'
include { MERGE_DECONV_RES_FULL   } from '../../modules/local/merge_deconv_res_full/main.nf'

workflow PREPARE_INPUTS {
    take:
    ch_samplesheet // channel: [meta, clean_bam, deconv_res]

    main:
    // --- Merge deconv_res per sample -----------------------------------------
    ch_samplesheet
        .map { meta, clean_bam, deconv_res ->
            [meta.id.toString(), meta, deconv_res]
        }
        .groupTuple(by: 0)
        .map { groupKey, meta, deconv_res ->
            def new_meta = meta.first()
            def deconv_res_list = deconv_res.unique { it.toString() } as List
            [new_meta, deconv_res_list, deconv_res_list.size()]
        }
        .branch {
            multi:  it[2] > 1
            single: it[2] == 1
        }
        .set { ch_deconv_branched }

    MERGE_DECONV_RES_FULL(
        ch_deconv_branched.multi.map { meta, deconv_res_list, count -> [meta, deconv_res_list] }
    )

    MERGE_DECONV_RES_FULL.out.merged_deconv_res
        .mix(
            ch_deconv_branched.single.map { meta, deconv_res_list, count ->
                [meta, deconv_res_list.first()]
            }
        )
        .map { meta, deconv_res -> [meta.id.toString(), meta, deconv_res] }
        .set { ch_final_deconv_res }

    // --- Merge + dedup clean BAMs per sample ---------------------------------
    ch_samplesheet
        .map { meta, clean_bam, deconv_res ->
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
        .map { meta, clean_bam -> [meta.id.toString(), meta, clean_bam] }
        .set { ch_final_bam }

    // --- Join into one 1:1:1 channel -----------------------------------------
    ch_final_deconv_res
        .join(ch_final_bam, by: 0)
        .map { groupKey, meta, deconv_res, meta2, clean_bam ->
            [meta, clean_bam, deconv_res]
        }
        .set { ch_prepared }

    emit:
    samplesheet = ch_prepared // channel: [meta, clean_bam, deconv_res]
}
