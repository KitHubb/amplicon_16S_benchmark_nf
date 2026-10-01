// Standalone consumer: run with -C to exclude the pipeline configuration.
include { QIIME_DADA2 as INDEPENDENT; QIIME_DADA2 as PSEUDO } from '../modules/qiime_dada2'
include { QIIME_TAXONOMY } from '../modules/qiime_taxonomy'
include { QIIME_PHYLOGENY } from '../modules/qiime_phylogeny'
include { QIIME_DIVERSITY } from '../modules/qiime_diversity'
include { TRIMM_OPTIMAL_WORKFLOW } from '../subworkflows/trimm_optimal'

workflow {
    def fixture = file(params.fixture, checkIfExists: true)
    def metadata = file(params.metadata, checkIfExists: true)
    def candidates = file(params.candidates, checkIfExists: true)
    def dada = [dada2_pooling_method: 'independent', dada2_trunc_len_f: 250, dada2_trunc_len_r: 200, dada2_max_ee_f: 3, dada2_max_ee_r: 5]
    def tax = [taxonomy_confidence: 0.83, optimal_min_sample_reads: 71]
    INDEPENDENT(dada, 'standalone', fixture)
    PSEUDO(dada + [dada2_pooling_method: 'pseudo'], 'standalone', fixture)
    QIIME_TAXONOMY(tax, 'standalone', PSEUDO.out.repseq, PSEUDO.out.table, metadata, fixture, 'TEST')
    QIIME_PHYLOGENY('standalone', PSEUDO.out.repseq)
    QIIME_DIVERSITY([sampling_depth: 37, alpha_max_depth: 83], 'standalone', PSEUDO.out.table, QIIME_PHYLOGENY.out.rooted_tree, metadata)
    TRIMM_OPTIMAL_WORKFLOW(Channel.value(fixture), Channel.value(fixture), Channel.value('TEST'), candidates,
        dada + [dada2_pooling_method: 'pseudo'], tax)
}
