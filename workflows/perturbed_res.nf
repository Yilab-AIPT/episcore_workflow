/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERGE_DECONV_RES } from '../modules/local/merge_deconv_res/main.nf'
include { REPLACE_DECONV_PROB } from '../modules/local/replace_deconv_prob/main.nf'
include { SPLIT_BAM_BY_DECONV_RES } from '../modules/local/split_bam_by_deconv_res/main.nf'
include { CALC_BETA_ZSCORE } from '../subworkflows/local/calc_beta_zscore.nf'
include { ESTIMATE_FF } from '../subworkflows/local/estimate_ff.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PERTURBED_RES {
    take:
    ch_samplesheet // channel: [meta(id=full_name, sample), clean_bam, perturbed_res, original_res]

    main:
    // 1. Merge the (shared) original deconv results per sample (e.g. 24 -> 12).
    //    Multiple perturbation conditions of a sample reference the same
    //    original files, so deduplicate before merging.
    ch_samplesheet
        .map { meta, clean_bam, perturbed_res, original_res ->
            [meta.sample, original_res]
        }
        .groupTuple(by: 0)
        .map { sample, original_res_list ->
            def uniq = original_res_list.unique { it.toString() } as List
            [[id: sample], uniq]
        }
        .set { ch_original_grouped }

    MERGE_DECONV_RES(
        ch_original_grouped,
        params.ncpgs
    )
    MERGE_DECONV_RES.out.merged_deconv_res
        .map { meta, merged -> [meta.id, merged] } // [sample, merged_original_res]
        .set { ch_merged_original }

    // 2. Collapse to one entry per perturbation condition (full_name).
    //    Each full_name appears once per original_res row; keep a single record
    //    carrying the sample key, clean_bam and perturbed_res.
    ch_samplesheet
        .map { meta, clean_bam, perturbed_res, original_res ->
            [meta.id, meta, clean_bam, perturbed_res]
        }
        .groupTuple(by: 0)
        .map { full_name, metas, clean_bams, perturbed_res_list ->
            def meta = metas.first()
            [meta.sample, meta, clean_bams.first(), perturbed_res_list.first()]
        }
        .set { ch_conditions } // [sample, meta, clean_bam, perturbed_res]

    // 3. Attach the per-sample merged original result to every condition and
    //    replace the original probabilities with the perturbed ones.
    ch_conditions
        .combine(ch_merged_original, by: 0)
        .map { sample, meta, clean_bam, perturbed_res, merged_original ->
            [meta, merged_original, perturbed_res]
        }
        .set { ch_replace_input } // [meta, original_res, perturbed_res]

    REPLACE_DECONV_PROB(ch_replace_input)
    REPLACE_DECONV_PROB.out.perturbed_deconv_res
        .set { ch_perturbed_deconv_res } // [meta, perturbed_deconv_res]

    // 4. Split the clean BAM by the perturbed deconv result. One clean BAM per
    //    sample, so no SAMTOOLS_MERGE / PICARD_MARKDUPLICATES is required.
    ch_conditions
        .map { sample, meta, clean_bam, perturbed_res -> [meta.id, clean_bam] }
        .set { ch_clean_bam_by_id }

    ch_perturbed_deconv_res
        .map { meta, deconv_res -> [meta.id, meta, deconv_res] }
        .join(ch_clean_bam_by_id, by: 0)
        .map { full_name, meta, deconv_res, clean_bam -> [meta, clean_bam, deconv_res] }
        .set { ch_split_input } // [meta, clean_bam, perturbed_deconv_res]

    SPLIT_BAM_BY_DECONV_RES(
        ch_split_input,
        params.threshold
    )
    SPLIT_BAM_BY_DECONV_RES.out.splitted_bam
        .set { ch_splitted_bam } // [meta, target_bam, background_bam]

    // Write a collated samplesheet of the splitted BAMs for downstream reuse.
    ch_splitted_bam
        .map { meta, target_bam, background_bam ->
            def target_bam_path = "${params.outdir}/split_bam_by_deconv_res/${target_bam.name}"
            def background_bam_path = "${params.outdir}/split_bam_by_deconv_res/${background_bam.name}"
            "${meta.id},${target_bam_path},${background_bam_path}\n"
        }
        .collectFile(
            name: 'splitted_bam.csv',
            storeDir: "${params.outdir}/samplesheet",
            seed: "sample,target_bam,background_bam\n"
        )

    // 5. Downstream: beta/z-score (episcore) and SNP-based fetal fraction.
    CALC_BETA_ZSCORE(ch_splitted_bam)
    CALC_BETA_ZSCORE.out.zscore
        .set { ch_beta_zscore }
    CALC_BETA_ZSCORE.out.beta_value
        .set { ch_beta_value }

    ESTIMATE_FF(ch_splitted_bam)
    ESTIMATE_FF.out.snp_pileup
        .set { ch_snp_pileup }
    ESTIMATE_FF.out.snp_ff
        .set { ch_snp_ff }

    emit:
    splitted_bam = ch_splitted_bam
    beta_value   = ch_beta_value
    zscore       = ch_beta_zscore
    snp_pileup   = ch_snp_pileup
    snp_ff       = ch_snp_ff
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
