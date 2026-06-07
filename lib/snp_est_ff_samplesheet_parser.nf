// lib/snp_est_ff_samplesheet_parser.nf

/*
Parses and validates the input samplesheet for the SNP-based fetal fraction
estimation workflow (workflows/snp_est_ff.nf).

Two entry points are supported:
  - est_ff_from_bam    : start from clean BAM + deconvolution result and build
                         the SNP pileup internally before estimating FF.
  - est_ff_from_pileup : start directly from a pre-computed merged pileup.

@param samplesheet_path Path to the samplesheet CSV file.
@param step             The current entry point (est_ff_from_bam | est_ff_from_pileup).
@return                 A Channel structure:
                          est_ff_from_bam    -> [ meta, clean_bam, deconv_res ]
                          est_ff_from_pileup -> [ meta, pileup ]
*/

def validateAndParseSnpFFSamplesheet(samplesheet_path, step) {

    // 1. Define required columns for different entry steps
    def required_columns = [
        'est_ff_from_bam'    : ['sample', 'clean_bam', 'deconv_res'],
        'est_ff_from_pileup' : ['sample', 'pileup']
    ]

    // Validate if the provided step is a valid entry point
    if (!required_columns.containsKey(step)) {
        error "Unknown step '${step}'. Valid options are: ${required_columns.keySet().join(', ')}"
    }

    // 2. Create Channel and parse CSV
    Channel.fromPath(samplesheet_path)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // 2.1 Basic Validation: Check if required columns exist
            def missing_cols = required_columns[step].findAll { !row.containsKey(it) }
            if (missing_cols) {
                error "Samplesheet missing required columns for step '${step}': ${missing_cols.join(', ')}"
            }

            // 2.2 Construct Meta Map
            def meta = [id: row.sample]

            // 2.3 Dispatch parsing logic based on the step
            if (step == 'est_ff_from_bam') {
                return parseEstFfFromBamRow(row, meta)
            } else if (step == 'est_ff_from_pileup') {
                return parseEstFfFromPileupRow(row, meta)
            } else {
                return null
            }
        }
}

// --- Internal Helper Functions ---

// Parse row for est_ff_from_bam step
// Returns: [ meta, clean_bam, deconv_res ]
def parseEstFfFromBamRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    def deconv_res = file(row.deconv_res, checkIfExists: true)
    return [ meta, clean_bam, deconv_res ]
}

// Parse row for est_ff_from_pileup step
// Returns: [ meta, pileup ]
def parseEstFfFromPileupRow(row, meta) {
    def pileup = file(row.pileup, checkIfExists: true)
    return [ meta, pileup ]
}
