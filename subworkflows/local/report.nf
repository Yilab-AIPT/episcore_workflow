//
// Generate final report from episcore and SNP-based fetal fraction results
//

include { GENERATE_REPORT } from '../../modules/local/generate_report/main.nf'
include { COLLECT_REPORTS } from '../../modules/local/collect_report/main.nf'

workflow REPORT {
    take:
    ch_episcore     // channel: [ meta, episcore_tsv ]
    ch_beta_value   // channel: [ meta, beta_value_tsv_gz ]
    ch_snp_pileup   // channel: [ meta, snp_pileup_tsv_gz ]
    ch_snp_ff       // channel: [ meta, snp_ff_tsv ]

    main:
    // Merge all channels by meta (sample)
    ch_episcore
        .join(ch_beta_value, by: 0)
        .join(ch_snp_pileup, by: 0)
        .join(ch_snp_ff, by: 0)
        .map { meta, episcore, beta_value, snp_pileup, snp_ff ->
            return [meta, episcore, beta_value, snp_pileup, snp_ff]
        }
        .set { ch_merged_input }

    // Generate individual reports per sample
    def meta_file = params.meta ? file(params.meta) : []
    GENERATE_REPORT(
        ch_merged_input,
        meta_file
    )
    
    // Collect all individual reports and merge into summary
    def ch_all_reports = GENERATE_REPORT.out.report
        .map { meta, report -> report }
        .collect()
    
    COLLECT_REPORTS(ch_all_reports)
    
    emit:
    // summary = COLLECT_REPORTS.out.summary
    report = GENERATE_REPORT.out.report
}
