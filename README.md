# amplicon_16S_benchmark_nf

Minimal Nextflow DSL2 pipeline for Stage 1 of the [2026 16S V1-V5 benchmark design](https://github.com/KitHubb/2026_16S_V1V5_test).

## Scope

The benchmark workflow performs:

1. FastQC/MultiQC on paired FASTQ files.
2. One shared Cutadapt preprocessing pass.
3. One shared QIIME 2 paired-end import.
4. One independent q2-dada2 task for every supported row in the benchmark parameter TSV.
5. Flat export of QIIME artifacts and plain TSV/BIOM/FASTA files for later scoring.

Taxonomy is not run in benchmark mode. This repository does not run tree construction or alpha/beta-diversity analysis.

## Stage 1 parameters

The parameter TSV must contain these columns:

```text
config_id  trunc_len_f  trunc_len_r  max_ee_f  max_ee_r  merge_id  merge_mode  min_overlap
```

The supplied design contains 165 rows:

- 55 `M12` rows: q2-dada2 with `minOverlap=12`.
- 55 `M10` rows: q2-dada2 with `minOverlap=10`.
- 55 `CAT` rows: reserved for standalone R/DADA2 `justConcatenate`.

CAT rows are validated but skipped until a pinned standalone R/DADA2 image and module are available. They are never replaced by an overlap setting.

The workflow fixes `truncQ=2`, `maxMismatch=0` (q2-dada2 default), pooling=`independent`, and chimera=`consensus`. M10/M12 rows must use `merge_mode=overlap`; CAT rows must use `merge_mode=just_concatenate`.

## Run

Keep run-specific paths in an external YAML file:

```yaml
benchmark_parameters: "<benchmark-parameters.tsv>"
reads: "<input-dir>/*_{R1,R2}_001.fastq.gz"
outdir: "<outdir>"
run_label: "benchmark"
region: "V1V3"
primer_f: "AGRGTTYGATYMTGGCTCAG"
primer_r: "CCGTCAATTCMTTTGAGTTT"
quality: 20
min_length: 50
discard_untrimmed: true
skip_cutadapt: false
skip_qc: false
stop_after: "none"
trimm_optimal: false
dada2_pooling_method: "independent"
qiime_sif: "file:///data/software/singularity/qiime2_amplicon_2025.7.sif"
```

Run in the foreground so Nextflow progress remains visible:

```bash
nextflow -log <outdir>/nextflow.log run <pipeline-dir>/main.nf \
  -profile singularity \
  -params-file <base.yaml> \
  -work-dir <work-dir> \
  -with-report <outdir>/nextflow-report.html \
  -with-timeline <outdir>/nextflow-timeline.html \
  -with-trace <outdir>/nextflow-trace.tsv \
  -with-dag <outdir>/nextflow-dag.html \
  -resume
```

Do not put results, work files, logs, or temporary files under `/data/software/nextflow/`.

## Outputs

Shared preprocessing outputs use the existing numbered stage directories. Each DADA2 configuration is published directly under `outdir` using an option-labelled directory such as:

```text
T01_E22_M12_F450_R450_EE2-2_overlap12/
```

There is no additional stage directory inside a configuration directory. Each successful configuration publishes:

- `table.qza` and `feature-table.biom`
- `rep-seqs.qza` and `rep-seqs.fasta`
- `denoising-stats.qza`, `denoising-stats.qzv`, and `denoising-stats.tsv`
- `table-summary.qzv` and `rep-seqs.qzv`
- `benchmark-parameters.tsv`
- `versions.yml`

Failed parameter tasks do not terminate the entire sweep. Use the Nextflow trace and report to identify failed or zero-output configurations.

## Containers and resources

Tool versions are pinned by the original pipeline configuration. A local image can be supplied with `qiime_sif`; otherwise the QIIME 2 amplicon 2025.7 image is used.

Default resource labels are defined in `conf/base.config`. Override executor, queue, CPU, memory, time, and concurrency in an external site config rather than committing project-specific settings.

## Validation

Non-executing workflow validation:

```bash
nextflow run . -preview -profile singularity \
  -params-file <base.yaml> \
  -work-dir <preview-work-dir>
```

Repository checks:

```bash
python3 -m json.tool nextflow_schema.json >/dev/null
nextflow config . >/dev/null
python3 tests/test_reuse.py
```

The preview and mocked reuse tests do not establish biological validity. Run a small representative end-to-end dataset before launching the complete benchmark.
