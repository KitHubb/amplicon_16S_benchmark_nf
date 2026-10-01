process TRIMMING_SUMMARY {
    label 'process_low'

    input:
    path summaries
    path summary_script

    output:
    path 'trimming_summary.tsv', emit: tsv
    path 'trimming_summary_mqc.json', emit: multiqc

    script:
    """
    python3 ${summary_script} *.trimming.summary.json
    """
}
