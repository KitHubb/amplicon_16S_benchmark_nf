process QIIME_IMPORT {
    tag run_label
    label 'process_medium'

    input:
    val run_label
    tuple val(sample_ids), path(forward_reads, stageAs: 'forward/read??/*'), path(reverse_reads, stageAs: 'reverse/read??/*')

    output:
    path 'manifest.tsv',      emit: manifest
    path 'demux-paired.qza',  emit: demux
    path 'demux-summary.qzv', emit: summary
    path 'versions.yml',      emit: versions

    script:
    def forward = forward_reads instanceof List ? forward_reads : [forward_reads]
    def reverse = reverse_reads instanceof List ? reverse_reads : [reverse_reads]
    def records = sample_ids.withIndex().collect { id, i -> [id, forward[i].toString(), reverse[i].toString()] }
    def encoded = groovy.json.JsonOutput.toJson(records).bytes.encodeBase64().toString()
    """
    export TMPDIR="\$PWD/tmp" XDG_CACHE_HOME="\$PWD/cache"
    export MPLCONFIGDIR="\$PWD/mpl"
    mkdir -p "\$TMPDIR" "\$XDG_CACHE_HOME" "\$MPLCONFIGDIR"

    python - <<'PYMANIFEST'
    import base64, csv, json
    from pathlib import Path
    rows = json.loads(base64.b64decode('${encoded}'))
    if not rows or len({r[0] for r in rows}) != len(rows):
        raise SystemExit('Import requires nonempty, unique sample IDs')
    with open('manifest.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter=chr(9))
        writer.writerow(['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath'])
        for sample, forward, reverse in rows:
            writer.writerow([sample, str(Path(forward).absolute()), str(Path(reverse).absolute())])
    PYMANIFEST

    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path manifest.tsv \
      --input-format PairedEndFastqManifestPhred33V2 \
      --output-path demux-paired.qza

    qiime demux summarize \
      --i-data demux-paired.qza \
      --o-visualization demux-summary.qzv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      qiime2: \$(qiime --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
