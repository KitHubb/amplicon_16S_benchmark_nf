process QIIME_DADA2_BENCHMARK {
    tag meta.config_id
    label 'process_high'

    input:
    tuple val(meta), val(options)
    path demux

    output:
    tuple val(meta),
          path('table.qza'),
          path('rep-seqs.qza'),
          path('denoising-stats.qza'),
          path('denoising-stats.qzv'),
          path('table-summary.qzv'),
          path('rep-seqs.qzv'),
          path('denoising-stats.tsv'),
          path('feature-table.biom'),
          path('rep-seqs.fasta'),
          path('benchmark-parameters.tsv'),
          path('versions.yml'),
          emit: results

    script:
    """
    export TMPDIR="\$PWD/tmp" NUMBA_CACHE_DIR="\$PWD/numba" MPLCONFIGDIR="\$PWD/mpl" XDG_CACHE_HOME="\$PWD/cache"
    mkdir -p "\$TMPDIR" "\$NUMBA_CACHE_DIR" "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    qiime dada2 denoise-paired \
      --p-pooling-method ${options.dada2_pooling_method} \
      --i-demultiplexed-seqs ${demux} \
      --p-trim-left-f 0 --p-trim-left-r 0 \
      --p-trunc-len-f ${options.dada2_trunc_len_f} \
      --p-trunc-len-r ${options.dada2_trunc_len_r} \
      --p-max-ee-f ${options.dada2_max_ee_f} \
      --p-max-ee-r ${options.dada2_max_ee_r} \
      --p-trunc-q 2 \
      --p-min-overlap ${options.min_overlap} \
      --p-chimera-method consensus \
      --p-n-threads ${task.cpus} \
      --o-table table.qza \
      --o-representative-sequences rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza

    qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
    qiime feature-table summarize --i-table table.qza --o-visualization table-summary.qzv
    qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

    qiime tools export --input-path denoising-stats.qza --output-path exported-stats
    qiime tools export --input-path table.qza --output-path exported-table
    qiime tools export --input-path rep-seqs.qza --output-path exported-repseqs
    mv exported-stats/stats.tsv denoising-stats.tsv
    mv exported-table/feature-table.biom feature-table.biom
    mv exported-repseqs/dna-sequences.fasta rep-seqs.fasta

    printf 'config_id\ttrunc_len_f\ttrunc_len_r\tmax_ee_f\tmax_ee_r\tmerge_id\tmerge_mode\tmin_overlap\ttrunc_q\tpooling\tchimera\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t2\tindependent\tconsensus\n' \
      '${meta.config_id}' '${options.dada2_trunc_len_f}' '${options.dada2_trunc_len_r}' \
      '${options.dada2_max_ee_f}' '${options.dada2_max_ee_r}' '${options.merge_id}' '${options.merge_mode}' '${options.min_overlap}' \
      > benchmark-parameters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      qiime2: \$(qiime --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
