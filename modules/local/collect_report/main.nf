process COLLECT_REPORTS {
    tag "collect_report"
    
    input:
    path(reports, stageAs: "report_?.tsv")
    
    output:
    path("summary_report.tsv"), emit: summary
    
    script:
    """
    # Extract header from first file
    head -n 1 ${reports[0]} > summary_report.tsv
    
    # Append all data rows (skip header) from all files
    for file in ${reports}; do
        tail -n +2 "\$file" >> summary_report.tsv
    done
    """
}