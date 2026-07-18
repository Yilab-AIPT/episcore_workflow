//
// AIPT ref-40 flexible per-chromosome episcore.
//
// Splits the clean BAM into one target BAM per distinct episcore probability
// threshold (plus background at params.threshold for raw_total depth filtering),
// extracts methylation with MethylDackel, filters CpGs by target+background
// depth, then computes episcore against the grid-search best combo + reference.
//

include { SPLIT_BAM_BY_THRESHOLDS } from '../../modules/local/split_bam_by_thresholds/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TARGET } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BACKGROUND } from '../../modules/nf-core/samtools/index/main.nf'
include { METHYLDACKEL_EXTRACT as METHYLDACKEL_TARGET } from '../../modules/nf-core/methyldackel/extract/main.nf'
include { METHYLDACKEL_EXTRACT as METHYLDACKEL_BACKGROUND } from '../../modules/nf-core/methyldackel/extract/main.nf'
include { FILTER_CPG_BY_DEPTH } from '../../modules/local/filter_cpg_by_depth/main.nf'
include { CALC_EPISCORE_FLEXIBLE } from '../../modules/local/calc_episcore_flexible/main.nf'

workflow AIPT_EPISCORE {
    take:
    ch_samplesheet      // channel: [meta, clean_bam, deconv_res]
    ep_thresholds       // val: comma-separated distinct episcore thresholds

    main:
    // Ensure depth-filter threshold (params.threshold, typically 0.5) is always
    // among the split targets so we have a target bedGraph to pair with background.
    def depth_thr = (params.threshold as BigDecimal).stripTrailingZeros().toPlainString()
    def split_thresholds = (ep_thresholds.split(',') as List)
    if (!split_thresholds.contains(depth_thr)) {
        split_thresholds = split_thresholds + [depth_thr]
    }
    def n_split_thresholds = split_thresholds.size()
    def split_thresholds_str = split_thresholds.join(',')

    SPLIT_BAM_BY_THRESHOLDS(ch_samplesheet, split_thresholds_str, params.threshold)

    // Fan out to one item per (sample, threshold) target BAM with a composite key.
    SPLIT_BAM_BY_THRESHOLDS.out.target_bams
        .flatMap { meta, bams ->
            def files = bams instanceof List ? bams : [bams]
            files.collect { bam ->
                def m = (bam.name =~ /__thr_([0-9.]+)_target/)
                def t = m ? m[0][1] : 'NA'
                def tmeta = [id: "${meta.id}__thr_${t}", sample: meta.id]
                return [tmeta, bam]
            }
        }
        .set { ch_target_bam }

    SAMTOOLS_INDEX_TARGET(ch_target_bam)
    ch_target_bam
        .join(SAMTOOLS_INDEX_TARGET.out.bai)
        .multiMap { tmeta, bam, bai ->
            bam_input: [tmeta, bam]
            bai_input: [tmeta, bai]
        }
        .set { ch_indexed_target }

    METHYLDACKEL_TARGET(
        ch_indexed_target.bam_input,
        ch_indexed_target.bai_input,
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)]
    )

    // Regroup per-threshold bedGraphs back to one item per sample.
    METHYLDACKEL_TARGET.out.bedgraph
        .map { tmeta, bedgraph -> [[id: tmeta.sample], bedgraph] }
        .groupTuple(size: n_split_thresholds)
        .set { ch_bedgraphs }

    // Background BAM at depth_thr → MethylDackel → join with matching target bedGraph.
    SPLIT_BAM_BY_THRESHOLDS.out.background_bams
        .flatMap { meta, bams ->
            def files = bams instanceof List ? bams : [bams]
            files.collect { bam ->
                def bmeta = [id: "${meta.id}__thr_${depth_thr}_bg", sample: meta.id]
                return [bmeta, bam]
            }
        }
        .set { ch_background_bam }

    SAMTOOLS_INDEX_BACKGROUND(ch_background_bam)
    ch_background_bam
        .join(SAMTOOLS_INDEX_BACKGROUND.out.bai)
        .multiMap { bmeta, bam, bai ->
            bam_input: [bmeta, bam]
            bai_input: [bmeta, bai]
        }
        .set { ch_indexed_background }

    METHYLDACKEL_BACKGROUND(
        ch_indexed_background.bam_input,
        ch_indexed_background.bai_input,
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)]
    )

    METHYLDACKEL_BACKGROUND.out.bedgraph
        .map { bmeta, bedgraph -> [[id: bmeta.sample], bedgraph] }
        .set { ch_background_bedgraph }

    // Target bedGraph at depth_thr only (for raw_total depth filter).
    METHYLDACKEL_TARGET.out.bedgraph
        .map { tmeta, bedgraph ->
            def m = (bedgraph.name =~ /__thr_([0-9.]+)_target/)
            def t = m ? m[0][1] : 'NA'
            return [[id: tmeta.sample], t, bedgraph]
        }
        .filter { _meta, t, _bedgraph -> t == depth_thr }
        .map { meta, _t, bedgraph -> [meta, bedgraph] }
        .join(ch_background_bedgraph)
        .set { ch_depth_pair }

    FILTER_CPG_BY_DEPTH(ch_depth_pair, params.beta_depth_threshold)

    def ref_matrix = file("${params.grid_search_result}/best_reference_matrix.tsv")
    def best_combo = file("${params.grid_search_result}/best_combo_episcore.csv")

    CALC_EPISCORE_FLEXIBLE(
        ch_bedgraphs.join(FILTER_CPG_BY_DEPTH.out.filtered_cpgs),
        best_combo,
        ref_matrix,
        file(params.cpg_recall_dir)
    )

    emit:
    episcore = CALC_EPISCORE_FLEXIBLE.out.episcore
}
