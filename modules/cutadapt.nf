process CUTADAPT {
    tag meta.id
    label 'process_medium'

    input:
    tuple val(meta), path(r1), path(r2)
    val trimming_options
    val primers
    path trimming_script

    output:
    tuple val(meta),
          path("${meta.id}.R1.trimmed.fastq.gz"),
          path("${meta.id}.R2.trimmed.fastq.gz"),
          emit: reads
    path "${meta.id}*.cutadapt.json", emit: json
    path "${meta.id}.trimming.summary.json", emit: summary
    path "${meta.id}.cutadapt.log",  emit: cutadapt_log
    path 'versions.yml',              emit: versions

    script:
    def config = trimming_options + [primers: primers]
    def encoded = groovy.json.JsonOutput.toJson(config).bytes.encodeBase64().toString()
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 ${quote.call(trimming_script)} --config '${encoded}' --sample ${quote.call(meta.id)} --cores ${task.cpus} ${quote.call(r1)} ${quote.call(r2)}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
