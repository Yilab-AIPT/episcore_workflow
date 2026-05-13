// lib/grid_search_parameters_parser.nf

/*
Parses and validates the inputs to the grid-search workflow.

@param samplesheet_path Path to the input samplesheet CSV file, with columns: sample, clean_bam, deconv_res
@return                 A Channel of tuples: [ meta, clean_bam, deconv_res ]
                        where meta = [id: sample]
*/

def validateAndParseGridSearchParameters(samplesheet_path) {
    // 1. Required columns for the grid-search samplesheet
    def required_columns = ['sample', 'clean_bam', 'deconv_res']

    // 2. Parse samplesheet
    return Channel.fromPath(samplesheet_path)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def missing_cols = required_columns.findAll { !row.containsKey(it) }
            if (missing_cols) {
                error "Samplesheet missing required columns: ${missing_cols.join(', ')}"
            }

            def meta       = [id: row.sample]
            def clean_bam  = file(row.clean_bam,  checkIfExists: true)
            def deconv_res = file(row.deconv_res, checkIfExists: true)
            return [ meta, clean_bam, deconv_res ]
        }
}
