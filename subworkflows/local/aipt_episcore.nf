//
// AIPT ref-40 flexible per-chromosome episcore.
//
// Splits the clean BAM into one target BAM per distinct episcore probability
// threshold, extracts methylation with MethylDackel, then computes the
// per-chromosome episcore against the grid-search best combo + reference matrix.
//

include { SPLIT_BAM_BY_THRESHOLDS } from '../../modules/local/split_bam_by_thresholds/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'
include { METHYLDACKEL_EXTRACT } from '../../modules/nf-core/methyldackel/extract/main.nf'
include { CALC_EPISCORE_FLEXIBLE } from '../../modules/local/calc_episcore_flexible/main.nf'

workflow AIPT_EPISCORE {
    take:
    ch_samplesheet      // channel: [meta, clean_bam, deconv_res]
    ep_thresholds       // val: comma-separated distinct episcore thresholds

    main:
    def n_episcore_thresholds = ep_thresholds.split(',').size()

    SPLIT_BAM_BY_THRESHOLDS(ch_samplesheet, ep_thresholds)

    // Fan out to one item per (sample, threshold) target BAM with a composite key.
    SPLIT_BAM_BY_THRESHOLDS.out.target_bams
        .transpose()
        .map { meta, bam ->
            def m = (bam.name =~ /__thr_([0-9.]+)_target/)
            def t = m ? m[0][1] : 'NA'
            def tmeta = [id: "${meta.id}__thr_${t}", sample: meta.id]
            return [tmeta, bam]
        }
        .set { ch_target_bam }

    SAMTOOLS_INDEX(ch_target_bam)
    ch_target_bam
        .join(SAMTOOLS_INDEX.out.bai)
        .multiMap { tmeta, bam, bai ->
            bam_input: [tmeta, bam]
            bai_input: [tmeta, bai]
        }
        .set { ch_indexed }

    METHYLDACKEL_EXTRACT(
        ch_indexed.bam_input,
        ch_indexed.bai_input,
        [[:], file(params.fasta)],
        [[:], file(params.fasta_index)]
    )

    // Regroup per-threshold bedGraphs back to one item per sample.
    // size must match the number of distinct episcore thresholds so each sample
    // emits as soon as its bedGraphs are ready (default groupTuple blocks until
    // every upstream sample completes MethylDackel).
    METHYLDACKEL_EXTRACT.out.bedgraph
        .map { tmeta, bedgraph -> [[id: tmeta.sample], bedgraph] }
        .groupTuple(size: n_episcore_thresholds)
        .set { ch_bedgraphs }

    def ref_matrix = file("${params.grid_search_result}/best_reference_matrix.tsv")
    def best_combo = file("${params.grid_search_result}/best_combo_episcore.csv")

    CALC_EPISCORE_FLEXIBLE(
        ch_bedgraphs,
        best_combo,
        ref_matrix,
        file(params.cpg_recall_dir),
        params.beta_depth_threshold
    )

    emit:
    episcore = CALC_EPISCORE_FLEXIBLE.out.episcore
}
