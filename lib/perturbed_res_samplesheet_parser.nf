// lib/perturbed_res_samplesheet_parser.nf

/*
Parses and validates the input samplesheet for the methylation-perturbation
workflow (workflows/perturbed_res.nf, params.step = 'perturbed_res').

The samplesheet records, for each sample, a single clean BAM, the original
deconvolution read-probability file(s), and a per-perturbation-condition file
whose reads have had their methylation status randomly perturbed and were then
re-inferred. Downstream the workflow merges the per-sample original results,
overwrites the original read probabilities with the perturbed ones, and runs
split_bam + beta/z-score + SNP FF for every perturbation condition.

Required columns:
  - sample        : sample identifier (one clean_bam per sample)
  - full_name     : unique perturbation-condition identifier ({sample}_*)
  - clean_bam     : clean BAM for the sample
  - perturbed_res : read-probability file with perturbed methylation status (txt/parquet)
  - original_res  : original read-probability file (txt/parquet)

@param samplesheet_path Path to the samplesheet CSV file.
@param step             The current entry point (perturbed_res).
@return                 A Channel structure: [ meta, clean_bam, perturbed_res, original_res ]
                        where meta = [id: full_name, sample: sample].
*/

def validateAndParsePerturbedResSamplesheet(samplesheet_path, step) {

    // 1. Define required columns for the supported entry step
    def required_columns = [
        'perturbed_res' : ['sample', 'full_name', 'clean_bam', 'perturbed_res', 'original_res']
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

            // 2.2 Construct Meta Map: id tracks the perturbation condition,
            //     sample groups the (shared) clean_bam and original_res files.
            def meta = [id: row.full_name, sample: row.sample]

            // 2.3 Parse row
            return parsePerturbedResRow(row, meta)
        }
}

// --- Internal Helper Functions ---

// Parse row for the perturbed_res step
// Returns: [ meta, clean_bam, perturbed_res, original_res ]
def parsePerturbedResRow(row, meta) {
    def clean_bam = file(row.clean_bam, checkIfExists: true)
    def perturbed_res = file(row.perturbed_res, checkIfExists: true)
    def original_res = file(row.original_res, checkIfExists: true)
    return [ meta, clean_bam, perturbed_res, original_res ]
}
