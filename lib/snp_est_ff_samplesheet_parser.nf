// lib/snp_est_ff_samplesheet_parser.nf

/*
Parses and validates the input samplesheet for the SNP-based fetal fraction
estimation workflow (workflows/snp_est_ff.nf).

Three entry points are supported:
  - est_ff_from_bam                      : clean BAM + deconvolution result;
                                           builds the SNP pileup internally.
  - est_ff_from_bam_without_deconv_res   : clean BAM only; merges multi-BAM
                                           samples before pileup + FF estimation.
  - est_ff_from_pileup                   : pre-computed merged pileup.

@param samplesheet_path Path to the samplesheet CSV file.
@param step             The current entry point.
@return                 A Channel structure:
                          est_ff_from_bam                      -> [ meta, clean_bam, deconv_res ]
                          est_ff_from_bam_without_deconv_res   -> [ meta, clean_bam ]
                          est_ff_from_pileup                   -> [ meta, pileup ]
*/

def validateAndParseSnpFFSamplesheet(samplesheet_path, step) {

    // 1. Define required columns for different entry steps
    def required_columns = [
        'est_ff_from_bam'                      : ['sample', 'clean_bam', 'deconv_res'],
        'est_ff_from_bam_without_deconv_res'   : ['sample', 'clean_bam'],
        'est_ff_from_pileup'                   : ['sample', 'pileup']
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
            } else if (step == 'est_ff_from_bam_without_deconv_res') {
                return parseEstFfFromBamWithoutDeconvResRow(row, meta)
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

// Parse row for est_ff_from_bam_without_deconv_res step
// Returns: [ meta, clean_bam ]
def parseEstFfFromBamWithoutDeconvResRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    return [ meta, clean_bam ]
}

// Parse row for est_ff_from_pileup step
// Returns: [ meta, pileup ]
def parseEstFfFromPileupRow(row, meta) {
    def pileup = file(row.pileup, checkIfExists: true)
    return [ meta, pileup ]
}
