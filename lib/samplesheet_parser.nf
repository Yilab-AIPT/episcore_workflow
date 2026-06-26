// lib/samplesheet_parser.nf

/*
Parses and validates the input samplesheet based on the requested entry step.
@param samplesheet_path Path to the samplesheet CSV file.
@param step             The current entry point (split_bam | episcore).
@return                 A Channel structure: [ meta, [files...] ]
*/

def validateAndParseSamplesheet(samplesheet_path, step) {
    
    // 1. Define required columns for different entry steps
    def required_columns = [
        'split_bam' : ['sample', 'clean_bam', 'deconv_res'],
        'episcore' : ['sample', 'target_bam', 'background_bam']
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
            if (step == 'split_bam') {
                return parseSplitBamRow(row, meta)
            } else if (step == 'episcore') {
                return parseEpiscoreRow(row, meta)
            } else {
                return null
            }
        }
}

// --- Internal Helper Functions ---

// Parse row for Split BAM Step
// Returns: [ meta, clean_bam, deconv_res ]
def parseSplitBamRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    def deconv_res = file(row.deconv_res, checkIfExists: true)
    return [ meta, clean_bam, deconv_res ]
}

// Parse row for Episcore Step
// Returns: [ meta, target_bam, background_bam ]
def parseEpiscoreRow(row, meta) {
    def target_bam = file(row.target_bam, checkIfExists: true)
    def background_bam = file(row.background_bam, checkIfExists: true)
    return [ meta, target_bam, background_bam ]
}
