// lib/samplesheet_parser.nf

/*
Parses and validates the input samplesheet based on the requested entry step.
@param samplesheet_path Path to the samplesheet CSV file.
@param step             The current entry point (preprocess | deconv | beta_zscore | rc_zscore).
@return                 A Channel structure: [ meta, [files...] ]
*/

def validateAndParseSamplesheet(samplesheet_path, step) {
    
    // 1. Define required columns for different entry steps
    def required_columns = [
        'preprocess': ['sample', 'fastq_1', 'fastq_2'],
        'deconv'  : ['sample', 'clean_bam'],
        'beta_zscore' : ['sample', 'target_bam', 'background_bam'],
        'split_bam' : ['sample', 'clean_bam', 'deconv_res'],
        'rc_zscore' : ['sample', 'clean_bam', 'deconv_res']
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
            // nf-core standard: includes id, single_end status, etc.
            def meta = [id: row.sample]

            // 2.3 Dispatch parsing logic based on the step
            if (step == 'preprocess') {
                return parsePreprocessRow(row, meta)
            } else if (step == 'deconv') {
                return parseDeconvRow(row, meta)
            } else if (step == 'split_bam') {
                return parseSplitBamRow(row, meta)
            } else if (step == 'beta_zscore') {
                return parseBetaZscoreRow(row, meta)
            } else if (step == 'rc_zscore') {
                return parseRCZScoreRow(row, meta)
            } else {
                return null
            }
        }
}

// --- Internal Helper Functions ---

// Parse row for Mapping Step
// Returns: [ meta, [fastq1, fastq2] ]
def parsePreprocessRow(row, meta) {
    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def fastq_2 = file(row.fastq_2, checkIfExists: true)
    
    def reads = [fastq_1, fastq_2]
    
    return [ meta, reads ]
}

// Parse row for Model Inference Step
// Returns: [ meta, clean_bam ]
def parseDeconvRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    // You can handle bam indices here if necessary
    return [ meta, clean_bam ]
}

// Parse row for Split BAM Step
// Returns: [ meta, target_bam, background_bam ]
def parseSplitBamRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    def deconv_res = file(row.deconv_res, checkIfExists: true)
    return [ meta, clean_bam, deconv_res ]
}

// Parse row for Beta Z-score Step
// Returns: [ meta, target_bam, background_bam ]
def parseBetaZscoreRow(row, meta) {
    def target_bam = file(row.target_bam, checkIfExists: true)
    def background_bam = file(row.background_bam, checkIfExists: true)
    return [ meta, target_bam, background_bam ]
}

// Parse row for RC Z-score Step
// Returns: [ meta, clean_bam, deconv_res ]
def parseRCZScoreRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    def deconv_res = file(row.deconv_res, checkIfExists: true)
    
    return [ meta, clean_bam, deconv_res ]
}