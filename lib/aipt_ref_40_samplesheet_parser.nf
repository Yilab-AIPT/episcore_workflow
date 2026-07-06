// lib/aipt_ref_40_samplesheet_parser.nf

/*
Parses and validates the AIPT ref-40 samplesheet.
Columns: sample, clean_bam, deconv_res (same layout as the NIPT split_bam step).
@param samplesheet_path Path to the samplesheet CSV file.
@return                 Channel structure: [ meta, clean_bam, deconv_res ]
*/

def validateAndParseAiptRef40Samplesheet(samplesheet_path) {
    def required_columns = ['sample', 'clean_bam', 'deconv_res']

    Channel.fromPath(samplesheet_path)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def missing_cols = required_columns.findAll { !row.containsKey(it) }
            if (missing_cols) {
                error "Samplesheet missing required columns for step 'aipt_ref_40': ${missing_cols.join(', ')}"
            }

            def meta = [id: row.sample]
            def clean_bam = file(row.clean_bam, checkIfExists: true)
            def deconv_res = file(row.deconv_res, checkIfExists: true)
            return [ meta, clean_bam, deconv_res ]
        }
}
