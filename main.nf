nextflow.enable.dsl = 2

params.reads = null
params.outdir = 'results'

include { QC as RAW_QC; QC as CLEAN_QC } from './subworkflows/qc'
include { TRIMMING_SUMMARY } from './modules/trimming_summary'
include { MULTIQC as TRIMMING_MULTIQC } from './modules/multiqc'
include { CUTADAPT } from './modules/cutadapt'
include { QIIME_IMPORT } from './modules/qiime_import'
include { QIIME_DADA2 } from './modules/qiime_dada2'
include { QIIME_TAXONOMY } from './modules/qiime_taxonomy'
include { QIIME_PHYLOGENY } from './modules/qiime_phylogeny'
include { QIIME_DIVERSITY } from './modules/qiime_diversity'
include { TRIMM_OPTIMAL_WORKFLOW } from './subworkflows/trimm_optimal'

include { PRIMER_TEST } from './tests/primers'

workflow {
    ['skip_cutadapt', 'skip_qc', 'qc_only'].each { key ->
        if (!(params[key].toString() in ['true', 'false'])) {
            error "${key} must be true or false"
        }
    }
    def trimming_options = TrimmingConfig.resolve(params)
    def dada2_options = [dada2_pooling_method: params.dada2_pooling_method, dada2_trunc_len_f: params.dada2_trunc_len_f, dada2_trunc_len_r: params.dada2_trunc_len_r, dada2_max_ee_f: params.dada2_max_ee_f, dada2_max_ee_r: params.dada2_max_ee_r]
    def taxonomy_options = [taxonomy_confidence: params.taxonomy_confidence, optimal_min_sample_reads: params.optimal_min_sample_reads]
    def diversity_options = [sampling_depth: params.sampling_depth, alpha_max_depth: params.alpha_max_depth]
    if (!(params.stop_after in ['none', 'raw_qc', 'cutadapt', 'clean_qc'])) {
        error 'stop_after must be one of: none, raw_qc, cutadapt, clean_qc'
    }
    if (params.qc_only.toString().toBoolean() && !(params.stop_after in ['none', 'clean_qc'])) {
        error 'qc_only is an alias for stop_after clean_qc; do not combine it with a different stop_after'
    }
    def stopAfter = params.qc_only.toString().toBoolean() ? 'clean_qc' : params.stop_after
    if (params.skip_qc.toString().toBoolean() && stopAfter in ['raw_qc', 'clean_qc']) {
        error 'raw_qc/clean_qc requires skip_qc=false'
    }
    if (params.skip_cutadapt.toString().toBoolean() && stopAfter in ['cutadapt', 'clean_qc']) {
        error 'cutadapt/clean_qc requires skip_cutadapt=false; use stop_after raw_qc for original reads'
    }
    if (!(params.dada2_pooling_method instanceof CharSequence) || !(params.dada2_pooling_method in ['independent', 'pseudo'])) {
        error 'dada2_pooling_method must be one of: independent, pseudo'
    }
    log.info "DADA2 pooling method: ${params.dada2_pooling_method}"
    if (params.test) {
        PRIMER_TEST()
    } else {
        primers = null
        if (!params.skip_cutadapt.toString().toBoolean() && stopAfter != 'raw_qc') {
            primers = PrimerConfig.resolve(params)
            log.info "Region: ${primers.region}; primer_f: ${primers.primer_f}; primer_r: ${primers.primer_r}; adapter_f: ${primers.adapter_f}; adapter_r: ${primers.adapter_r}"
        }
        if (!params.reads) {
            error 'Provide --reads, e.g. /data/run/*_{1,2}.fastq.gz'
        }
        if (stopAfter == 'none' && !params.metadata) {
            error 'Provide --metadata with a QIIME 2-compatible TSV file'
        }

        if (stopAfter == 'none' && !params.classifier) {
            error 'Provide --classifier with a QIIME 2 classifier artifact'
        }

        reads_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true, flat: true)
            .map { id, r1, r2 -> tuple([id: id.toString(), single_end: false], r1, r2) }

        if (!params.skip_qc.toString().toBoolean()) {
            RAW_QC(reads_ch, '01_raw_qc')
        }
        if (stopAfter != 'raw_qc') {
            if (params.skip_cutadapt.toString().toBoolean()) {
                log.info 'Skipping Cutadapt; importing original paired reads'
                import_reads = reads_ch
            } else {
                CUTADAPT(reads_ch, trimming_options, Channel.value(primers), file("${projectDir}/bin/trim_reads.py"))
                TRIMMING_SUMMARY(CUTADAPT.out.summary.collect(), file("${projectDir}/bin/summarize_trimming.py"))
                if (!params.skip_qc.toString().toBoolean()) {
                    TRIMMING_MULTIQC(TRIMMING_SUMMARY.out.multiqc, '02_cutadapt_q20')
                }
                import_reads = CUTADAPT.out.reads
                if (!params.skip_qc.toString().toBoolean() && stopAfter != 'cutadapt') {
                    CLEAN_QC(import_reads, '03_clean_qc')
                }
            }
            if (stopAfter == 'none') {
                QIIME_IMPORT(params.run_label, import_reads.collect(flat: false).map { rows ->
                    tuple(rows.collect { it[0].id }, rows.collect { it[1] }, rows.collect { it[2] })
                })

                metadata_ch = Channel.value(file(params.metadata, checkIfExists: true))
                classifier_file = file(params.classifier, checkIfExists: true)
                classifier_ch = Channel.value(classifier_file)

                taxonomy_label = params.taxonomy_label?.toString()?.trim()
                if (!taxonomy_label || taxonomy_label.equalsIgnoreCase('AUTO')) {
                    classifier_path_lower = classifier_file.toString().toLowerCase()
                    taxonomy_label = classifier_path_lower.contains('/gg2/')   ? 'GG2'   :
                                     classifier_path_lower.contains('/silva/') ? 'SILVA' :
                                     classifier_path_lower.contains('/gtdb')   ? 'GTDB'  :
                                     classifier_file.baseName
                }
                taxonomy_label = taxonomy_label.replaceAll(/[^A-Za-z0-9._-]/, '_')
                taxonomy_label_ch = Channel.value(taxonomy_label)

                // 5. Standard DADA2, or optional parameter sweep followed by rule-based selection.
                if (params.trimm_optimal.toString().toBoolean()) {
                    TRIMM_OPTIMAL_WORKFLOW(QIIME_IMPORT.out.demux, classifier_ch, taxonomy_label_ch,
                        file(params.trimm_combinations ?: "${projectDir}/params/trimming/${params.region.toString().trim().toLowerCase()}_10bp.tsv", checkIfExists: true),
                        dada2_options, taxonomy_options)
                    selected_table_ch = TRIMM_OPTIMAL_WORKFLOW.out.table
                    selected_repseq_ch = TRIMM_OPTIMAL_WORKFLOW.out.repseq
                } else {
                    QIIME_DADA2(dada2_options, params.run_label, QIIME_IMPORT.out.demux)
                    selected_table_ch = QIIME_DADA2.out.table
                    selected_repseq_ch = QIIME_DADA2.out.repseq
                }

                // 6. Final SILVA taxonomy assignment and taxa bar plot using the selected ASVs.
                QIIME_TAXONOMY(
                    taxonomy_options,
                    params.run_label,
                    selected_repseq_ch,
                    selected_table_ch,
                    metadata_ch,
                    classifier_ch,
                    taxonomy_label_ch
                )

                // 7. MAFFT alignment, masking, FastTree and midpoint-rooted phylogeny.
                QIIME_PHYLOGENY(params.run_label, selected_repseq_ch)

                // 8. Optional phylogenetic alpha/beta diversity, rarefaction, PCoA and Emperor plots.
                if (params.diversity_enabled.toString().toBoolean()) {
                    QIIME_DIVERSITY(
                        diversity_options,
                        params.run_label,
                        selected_table_ch,
                        QIIME_PHYLOGENY.out.rooted_tree,
                        metadata_ch
                    )
                }
            }
        }
    }
}
