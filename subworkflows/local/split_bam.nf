//
// Hard filter the clean bam by deconvolution result (threshold = 0.5)
//

include { SPLIT_BAM_BY_DECONV_RES } from '../../modules/local/split_bam_by_deconv_res/main.nf'
include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main.nf'
include { MERGE_DECONV_RES } from '../../modules/local/merge_deconv_res/main.nf'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main.nf'

workflow SPLIT_BAM {
    take:
    ch_samplesheet // channel: deconv results with columns ['sample', 'clean_bam', 'deconv_res']

    main:
    // Group by sample and check if merge is needed for deconv results
    ch_samplesheet.map {
        meta, clean_bam, deconv_res ->
        def groupKey = meta.id.toString()
        return [groupKey, meta, deconv_res]
    }.groupTuple(by: 0)
    .map { groupKey, meta, deconv_res ->
        def new_meta = meta.first()
        def deconv_res_list = deconv_res.unique { it.toString() } as List
        return [new_meta, deconv_res_list, deconv_res_list.size()]
    }
    .branch {
        multi: it[2] > 1    // Multiple deconv_res files need merging
        single: it[2] == 1  // Single deconv_res file, no merge needed
    }
    .set { ch_deconv_res_branched }

    // Merge multiple deconv_res files
    MERGE_DECONV_RES(
        ch_deconv_res_branched.multi.map { meta, deconv_res_list, count -> [meta, deconv_res_list] },
        params.ncpgs
    )
    
    // Combine merged and single deconv_res into one channel
    MERGE_DECONV_RES.out.merged_deconv_res
        .mix(
            ch_deconv_res_branched.single.map { meta, deconv_res_list, count -> 
                [meta, deconv_res_list.first()]
            }
        )
        .map { meta, deconv_res ->
            def groupKey = meta.id.toString()
            return [groupKey, meta, deconv_res]
        }
        .set { ch_final_deconv_res }

    // Group by sample and check if merge is needed for clean bam files
    ch_samplesheet.map {
        meta, clean_bam, deconv_res ->
        def groupKey = meta.id.toString()
        return [groupKey, meta, clean_bam]
    }.groupTuple(by: 0)
    .map { groupKey, meta, clean_bam ->
        def new_meta = meta.first()
        def clean_bam_list = clean_bam.unique { it.toString() } as List
        return [new_meta, clean_bam_list, clean_bam_list.size()]
    }
    .branch {
        multi: it[2] > 1    // Multiple BAM files need merging
        single: it[2] == 1  // Single BAM file, no merge needed
    }
    .set { ch_clean_bam_branched }

    // Merge multiple BAM files
    SAMTOOLS_MERGE(
        ch_clean_bam_branched.multi.map { meta, clean_bam_list, count -> [meta, clean_bam_list] },
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)],
        [[:], []]
    )
    
    // Mark duplicates on merged BAMs
    PICARD_MARKDUPLICATES(
        SAMTOOLS_MERGE.out.bam,
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)],
    )
    
    // Combine merged+deduped BAMs and single BAMs into one channel
    PICARD_MARKDUPLICATES.out.bam
        .mix(
            ch_clean_bam_branched.single.map { meta, clean_bam_list, count -> 
                [meta, clean_bam_list.first()]
            }
        )
        .map { meta, clean_bam ->
            def groupKey = meta.id.toString()
            return [groupKey, meta, clean_bam]
        }
        .set { ch_final_bam }

    // Join deconv_res and BAM by sample ID
    ch_final_deconv_res
        .join(ch_final_bam, by: 0)
        .map { groupKey, meta, deconv_res, meta_2, clean_bam ->
            return [meta, clean_bam, deconv_res]
        }
        .set { ch_final_samplesheet }

    SPLIT_BAM_BY_DECONV_RES(
        ch_final_samplesheet,
        params.threshold
    )
    SPLIT_BAM_BY_DECONV_RES.out.splitted_bam
        .set { ch_splitted_bam }

    emit:
    splitted_bam = ch_splitted_bam
}