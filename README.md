# amplicon_16S_benchmark_nf

This repository keeps the original `amplicon_16S_qiime_nf` workflow and adds a minimal Stage 1 benchmark mode for the [2026 16S V1-V5 benchmark design](https://github.com/KitHubb/2026_16S_V1V5_test).

When `benchmark_parameters` is supplied, FastQC, Cutadapt and QIIME 2 import run once. The imported reads are then fanned out to one q2-dada2 task per M10/M12 row. Taxonomy, phylogeny and diversity are intentionally skipped in this mode. CAT (`justConcatenate`) rows are also skipped until the standalone R/DADA2 image and module described by the study design are available; they are never substituted with an overlap setting.

```bash
nextflow run <pipeline-dir>/main.nf \
  -profile singularity \
  -params-file <base.yaml> \
  --benchmark_parameters <benchmark-parameters.tsv> \
  --outdir <outdir> \
  -work-dir <work-dir> \
  -with-report <outdir>/nextflow-report.html \
  -with-timeline <outdir>/nextflow-timeline.html \
  -with-trace <outdir>/nextflow-trace.tsv \
  -with-dag <outdir>/nextflow-dag.html \
  -resume
```

Benchmark outputs use one flat option-labelled directory per configuration, for example `OUTDIR/T01_E22_M12_F450_R450_EE2-2_overlap12/`. No additional stage directory is created inside a configuration directory. Each directory contains QIIME 2 artifacts and visualizations plus `denoising-stats.tsv`, `feature-table.biom`, `rep-seqs.fasta`, `benchmark-parameters.tsv`, and tool versions. The TSV must use unique safe `config_id` values and the columns `trunc_len_f`, `trunc_len_r`, `max_ee_f`, `max_ee_r`, `merge_id`, `merge_mode`, and `min_overlap`. Stage 1 fixes `truncQ=2`, pooling=`independent`, chimera=`consensus`, and the q2-dada2 default `maxMismatch=0`.

[English](#english) · [한국어](#korean)

<a id="english"></a>

## English

### Contents

1. [Quick start](#en-start)
2. [FASTQ and metadata](#en-inputs)
3. [Region and primers](#en-primers)
4. [Truncation optimization](#en-optimization)
5. [Containers and cache](#en-containers)
6. [Server resources](#en-resources)
7. [Analysis settings and customization](#en-options)
8. [Results](#en-outputs)
9. [Tests and releases](#en-tests)
10. [Code structure](#en-structure)
11. [Additional trimming and QC-only runs](#en-trimming)

<a id="en-start"></a>

### 1. Quick start

This pipeline analyzes Illumina paired-end 16S rRNA reads from the V1V3 or V3V4 region. It can compare DADA2 forward/reverse truncation lengths, select the highest-ranked completed candidate, and use its ASVs for downstream analysis.

**Workflow:** FASTQ → FastQC/MultiQC → Cutadapt → FastQC/MultiQC → QIIME 2 import → DADA2 (fixed lengths or candidate comparison) → taxonomy → MAFFT/FastTree → diversity.

Use Linux, Java 17+, Nextflow 26.04.2 or later (tested with 26.04.2), and one container runtime: Docker, Singularity, or Apptainer. Prepare your FASTQ files, metadata TSV, and a classifier QZA compatible with QIIME 2 amplicon 2025.7. SILVA, GG2, or other compatible classifiers can be used.

Create an external `analysis.yml` (the filename `nextflow_params.yml` also works). Replace the example paths with your own:

```yaml
region: V3V4
reads: "/path/to/reads/*_{1,2}.fastq.gz"
metadata: "/path/to/metadata.tsv"
classifier: "/path/to/classifier.qza"
outdir: "/path/to/results"
run_label: "experiment01"
trimm_optimal: true
taxonomy_label: SILVA
taxonomy_confidence: 0.7
diversity_enabled: true
sampling_depth: 1000
```

Run the released pipeline directly from GitHub; cloning the repository is optional:

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile singularity \
  -params-file /path/to/analysis.yml \
  -work-dir /path/to/work
```

Choose `-profile docker`, `-profile singularity` (the default runtime), or `-profile apptainer`. Images are downloaded on first use and cached. Analysis tools do not need separate host installations.

For a local checkout:

```bash
nextflow run /path/to/amplicon_16S_qiime_nf/main.nf \
  -profile singularity -params-file /path/to/analysis.yml \
  -work-dir /path/to/work
```

Use `region: V1V3` for V1V3 data. Leave `trimm_combinations` unset to use the region's default candidates. This example enables optimization; the code default for `trimm_optimal` is `false`.

Set `diversity_enabled: false` to skip diversity analysis. Adjust `sampling_depth` to your sample read counts. To resume an analysis, keep the same launch directory and work directory and add `-resume`. Changed inputs, options, or module structure can cause tasks to run again. For the latest development code, replace `-r v0.1.0` with `-r main -latest`.

<a id="en-inputs"></a>

### 2. FASTQ and metadata

The workflow pairs reads using the filename pattern in `reads` and takes sample annotations from `metadata`. It does **not** accept a CSV samplesheet of FASTQ paths through `--input` or `--samplesheet`. The QIIME 2 import manifest is generated automatically after Cutadapt.

Prepare one R1/R2 pair per sample:

```text
reads/
  S01_1.fastq.gz
  S01_2.fastq.gz
  S02_1.fastq.gz
  S02_2.fastq.gz
```

For these files, use `reads: "/path/to/reads/*_{1,2}.fastq.gz"`. The extracted sample IDs are `S01` and `S02`. For filenames ending in `_R1_001.fastq.gz` and `_R2_001.fastq.gz`, use `*_{R1,R2}_001.fastq.gz` instead.

Prepare a tab-separated metadata file with `sample-id` as the first column:

```tsv
sample-id	group
S01	control
S02	treatment
```

To create this example with actual tab characters:

```bash
printf 'sample-id\tgroup\nS01\tcontrol\nS02\ttreatment\n' > metadata.tsv
```

Set `metadata` in your analysis YAML to the file's absolute path. Sample IDs must be unique and exactly match the IDs extracted from the filenames. Replace or extend the example `group` column with your experiment's conditions and batches. Use tabs, not commas or spaces. Run V1V3 and V3V4 separately with their corresponding reads and metadata.

<a id="en-primers"></a>

### 3. Region and primers

The default region is V1V3.

| Region | Forward primer (5′→3′) | Reverse primer (5′→3′) |
|---|---|---|
| V1V3 | `AGAGTTTGATCCTGGCTCAG` | `ATTACCGCGGCTGCTGG` |
| V3V4 | `CCTACGGGNGGCWGCAG` | `GACTACHVGGGTATCTAATCC` |

Choose `--region V3V4`, use a YAML preset from [params/regions](params/regions), or add explicit primer options to your run command:

```bash
--region V3V4 --primer_f 'CCTACGGGNGGCWGCAG' --primer_r 'GACTACHVGGGTATCTAATCC'
```

**Priority: command line > external parameter YAML > region defaults.** Overriding one primer leaves the other at its region default. See [custom_primers.yml](params/examples/custom_primers.yml) for a YAML example.

Primers accept IUPAC DNA letters, including lowercase. Read-through adapters are the reverse complements of the opposite primers unless `adapter_f` or `adapter_r` is set explicitly. Default V1V3 primers retain the adapters `CCAGCAGCCGCGGTAAT` and `CTGAGCCAGGATCAAACTCT`. The final sequences are printed in the run log.

<a id="en-optimization"></a>

### 4. Truncation optimization

`dada2_pooling_method` selects `independent` (default) or `pseudo` for both standard DADA2 and every trimming candidate. Other values fail at workflow startup. Set it in your parameter YAML or pass `--dada2_pooling_method pseudo`. QIIME 2 pseudo-pooling is not R DADA2 full `pool=TRUE`. Existing output paths are unchanged; choose different `outdir` values when comparing runs.

Set `trimm_optimal: true` to enable candidate comparison. Truncation lengths apply to the reads imported into QIIME 2 after Cutadapt.

1. Run DADA2 independently for each forward/reverse pair in the candidate TSV.
2. Classify candidate ASVs and calculate read-retention and species-assignment metrics.
3. Rank completed candidates by the criteria below, in descending order. Later criteria break ties.
4. Pass the highest-ranked feature table and representative sequences to taxonomy, phylogeny, and diversity analysis.

| Priority | Metric | Calculation |
|---|---|---|
| 1 | `species_input_yield` | Species-assigned reads / DADA2 input reads × 100 |
| 2 | `species_read_pct` | Species-assigned reads / final feature-table reads × 100 |
| 3 | `nonchimeric_pct` | Non-chimeric reads / DADA2 input reads × 100 |

“Optimal” means best among the supplied candidates for the chosen classifier and ranking criteria. It does not establish a universal biological optimum. Metrics are aggregated across the dataset, rather than ranked per sample. Species detection uses `s__` or `s:` taxonomy markers and excludes names such as unknown and uncultured.

The region selects [v1v3_10bp.tsv](params/trimming/v1v3_10bp.tsv) (6 candidates) or [v3v4_10bp.tsv](params/trimming/v3v4_10bp.tsv) (10 candidates). Override it with `trimm_combinations: "/path/to/candidates.tsv"` or `--trimm_combinations /path/to/candidates.tsv`:

```tsv
name	trunc_len_f	trunc_len_r
F280_R270	280	270
F280_R260	280	260
F270_R250	270	250
```

Use actual tabs and unique candidate names. Current validation requires **F > R** and both lengths in **10-nt increments**. Default candidate design uses F + R − 20 ≥ 500 bp for V1V3 and ≥ 460 bp for V3V4. Both standard and optimization DADA2 use a minimum overlap of 12 bp. These design bounds are not verified biological maxima; adjust candidates to your read lengths and quality. See the [candidate design notes](params/trimming/README.md) for sources, length assumptions, and overlap calculations.

Candidates with failed DADA2 or taxonomy tasks are excluded. Selection fails if no candidates succeed. `optimal_min_sample_reads` defaults to 10000: low-read and zero-merge sample counts are reported but do not filter or rank candidates.

For fixed-length analysis, set `trimm_optimal: false` and use `dada2_trunc_len_f` and `dada2_trunc_len_r`. Their default value, `0`, disables fixed-length truncation. Optimization uses the TSV lengths instead of these two parameters.

<a id="en-containers"></a>

### 5. Containers and cache

Each tool uses a pinned public image unless you supply a local override.

| Tool | Public image | Local SIF/IMG parameter |
|---|---|---|
| FastQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | `fastqc_sif` |
| MultiQC | `quay.io/biocontainers/multiqc:1.27--pyhdfd78af_0` | `multiqc_sif` |
| Cutadapt | `quay.io/biocontainers/cutadapt:5.2--py311hc303176_2` | `cutadapt_sif` |
| QIIME 2 | `quay.io/qiime2/amplicon:2025.7` | `qiime_sif` |

FastQC and MultiQC use separate images and separate settings. The combined QC image and shared path setting are no longer used. QIIME 2 analysis steps share one QIIME 2 image.

To use existing images, add any of these settings to `analysis.yml` or `nextflow_params.yml`:

```yaml
fastqc_sif: "/absolute/path/to/fastqc.sif"
multiqc_sif: "/absolute/path/to/multiqc.img"
cutadapt_sif: "/absolute/path/to/cutadapt.sif"
qiime_sif: "/absolute/path/to/qiime2.sif"
```

Use `-profile singularity` or `-profile apptainer` for local SIF/IMG files. The image must be readable by that runtime; Docker cannot use these files. Omit a tool's setting to download its public image automatically. A local override takes priority for that tool only.

First use requires registry access and enough disk space, especially for QIIME 2. Singularity/Apptainer image references use `docker://`; Docker references omit that prefix. See the [Nextflow container documentation](https://docs.seqera.io/nextflow/container/singularity).

Set a persistent Singularity image cache before running:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/container-cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
```

For Apptainer, use `NXF_APPTAINER_CACHEDIR`. Docker manages its own image storage. On HPC, use a writable cache shared by all compute nodes. Image cache reuse avoids downloading images again; `-resume` separately controls reuse of completed task results.

**Existing lab server configuration.** Add these QC paths to your analysis YAML while retaining the reads, metadata, classifier, and other analysis settings:

```yaml
fastqc_sif: "/data/software/singularity/quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"
multiqc_sif: "/data/software/singularity/quay.io-biocontainers-multiqc-1.27--pyhdfd78af_0.img"
```

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile singularity -params-file /path/to/nextflow_params.yml
```

This server also has `/data/software/singularity/amplicon_16S_qc.config` with the same QC paths. To test it:

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile test,singularity \
  -c /data/software/singularity/amplicon_16S_qc.config
```

These absolute paths and the config file are specific to this server. Other users should set their own local paths or omit overrides to use automatic downloads.

<a id="en-resources"></a>

### 6. Server resources

Keep data paths and analysis parameters in an external YAML. Use an external `site.config` for the executor, queue, CPU, memory, and time limits. The default executor in [nextflow.config](nextflow.config) is `local`.

[conf/base.config](conf/base.config) defines these per-task defaults:

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_low` | 8 | 16 GB | 12 h |
| `process_medium` | 8 | 24 GB | 24 h |
| `process_high` | 16 | 64 GB | 72 h |

Optimization candidates can run concurrently. Adjust resources for your server and dataset. Example `site.config`:

```groovy
process {
    executor = 'local'
    maxForks = 1
    withLabel: process_low {
        cpus = 2
        memory = '4 GB'
    }
    withLabel: process_medium {
        cpus = 4
        memory = '8 GB'
    }
    withLabel: process_high {
        cpus = 8
        memory = '32 GB'
    }
}
```

Add `-c /path/to/site.config` to the run command. The values above illustrate the syntax; they are not suitable for every dataset. `maxForks` limits concurrency per process, not across all processes combined. Set the scheduler and queue to match your institution's HPC configuration. Set the work directory with `-work-dir`.

<a id="en-options"></a>

### 7. Analysis settings and customization

Optional stage skipping (both default to `false`):

```yaml
skip_cutadapt: true
skip_qc: true
```

`skip_cutadapt` bypasses all Cutadapt primer/adapter removal, quality trimming, and minimum-length filtering; original paired reads go directly to QIIME import. `discard_untrimmed` controls untrimmed-pair filtering (default `true`); see section 11 for additional trimming and QC-only runs. `skip_qc` disables all FastQC/MultiQC tasks, including the trimming summary report. If only Cutadapt is skipped, raw QC runs once and post-Cutadapt QC is omitted. Sample IDs still come from the `reads` pairing pattern; no trimmed filename suffix is required. Skipped stages do not produce new outputs; use a fresh output directory to avoid confusing old reports with the current run.

Use parameters first. Edit tool modules when you need to change commands or logic.

| Purpose | Parameters or defaults | Implementation |
|---|---|---|
| FastQC / MultiQC | Image overrides and CPU settings | [fastqc.nf](modules/fastqc.nf), [multiqc.nf](modules/multiqc.nf) |
| Region, primers, adapters | `region`, `primer_f`, `primer_r`, `adapter_f`, `adapter_r` | [PrimerConfig.groovy](lib/PrimerConfig.groovy) |
| Cutadapt filtering | `quality: 20`, `min_length: 50` | [cutadapt.nf](modules/cutadapt.nf); primer pass uses `-n 2`; `discard_untrimmed: true` by default |
| QIIME import | `reads` filename pattern | [qiime_import.nf](modules/qiime_import.nf) |
| Fixed DADA2 | `dada2_trunc_len_f: 0`, `dada2_trunc_len_r: 0`, `dada2_max_ee_f: 2`, `dada2_max_ee_r: 4` | [qiime_dada2.nf](modules/qiime_dada2.nf) |
| Candidate validation | `trimm_optimal`, `trimm_combinations` | [subworkflows/trimm_optimal.nf](subworkflows/trimm_optimal.nf), [params/trimming](params/trimming) |
| Candidate execution | The same `dada2_max_ee_*` settings | `DADA2_TRIMM_SWEEP` in [modules/trimm_optimal.nf](modules/trimm_optimal.nf) |
| Ranking and species rules | `optimal_min_sample_reads: 10000` changes reporting only | `TAXONOMY_TRIMM_SWEEP`, `SELECT_TRIMM_OPTIMAL` in the same module |
| Taxonomy | `classifier`, `taxonomy_confidence: 0.7`, `taxonomy_label: AUTO` | [qiime_taxonomy.nf](modules/qiime_taxonomy.nf) and candidate taxonomy in `modules/trimm_optimal.nf` |
| MAFFT / FastTree | No separate analysis parameters | [qiime_phylogeny.nf](modules/qiime_phylogeny.nf) |
| Diversity / rarefaction | `diversity_enabled: true`, `sampling_depth: 1000`, `alpha_max_depth: 10000` | [qiime_diversity.nf](modules/qiime_diversity.nf) |

For hard-coded DADA2 options such as `trunc-q` or chimera handling, check both standard and optimization modules. Keep candidate and final taxonomy commands consistent. When adding a parameter, update [nextflow.config](nextflow.config), [nextflow_schema.json](nextflow_schema.json), the relevant module, and both README language sections.

When changing image versions, check command options and classifier compatibility. Record the image versions and classifier used for the final analysis.

<a id="en-outputs"></a>

### 8. Results

All results are written under `outdir`.

| Directory | Contents |
|---|---|
| `01_raw_qc/`, `03_clean_qc/` | FastQC/MultiQC before and after filtering |
| `02_cutadapt_q20/` | Trimmed FASTQ, pass logs/JSON, summary TSV and trimming MultiQC; the directory name is fixed even if `quality` changes |
| `04_qiime2_import/` | Manifest, paired-end QZA, quality summary |
| `05_dada2/` | Standard DADA2 results and summaries when optimization is disabled |
| `05_trimm_optimal_dada2/` | DADA2 outputs for each candidate |
| `05_trimm_optimal_<LABEL>/` | Candidate taxonomy and metrics |
| `05_trimm_optimal_<LABEL>/selected/` | Full ranking, selected lengths, selected QZA files |
| `06_taxonomy/` | Outputs such as `taxonomy_<LABEL>.qza` and `taxa-bar-plots_<LABEL>.qzv` |
| `07_phylogeny/`, `08_diversity/` | Phylogeny and, if enabled, diversity results |

After optimization, review:

- `all_parameter_results.tsv`: rankings and metrics for successful candidates.
- `optimal_selection.tsv` and `optimal_truncation.txt`: selected F/R lengths and scores.
- `selected-table.qza`, `selected-rep-seqs.qza`, and `selected-denoising-stats.qza`: artifacts passed to downstream analysis.

`taxonomy_label: AUTO` infers labels such as SILVA, GG2, or GTDB from the classifier path. Set an explicit label to distinguish analyses. Optimization does not generate the standard-mode `05_dada2/table-summary.qzv`; inspect the selected table and stats when choosing diversity depth.

Keep raw reads, QIIME artifacts, container images, work directories, and participant metadata out of Git.

<a id="en-tests"></a>

### 9. Tests and releases

Run the bundled synthetic test without biological input data, metadata, or a classifier:

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile test,docker
```

Use `test,singularity` or `test,apptainer` for those runtimes. In a local checkout, also run the output verifier:

```bash
nextflow run . -profile test,singularity
python3 tests/check_containers.py results/primer-test
```

A successful run completes five tasks across four images:

| Task | Check |
|---|---|
| Cutadapt | Primer assertions and actual trimming; each mate must contain one 60-base T sequence |
| FastQC | Generate reports from the trimmed reads |
| MultiQC | Combine FastQC results into a report |
| QIIME import | Import paired reads and generate a demux summary |
| QIIME plugin check | Load required DADA2, feature-classifier, phylogeny, diversity, and taxa actions |

Reports are in `results/primer-test/container_qc/`, QIIME artifacts in `results/primer-test/04_qiime2_import/`, and the plugin log in `results/primer-test/container_checks/qiime-container.txt`.

Omit local image overrides to check public image retrieval. Set a new writable image-cache directory to test a cold cache. Rerun without `-resume` to execute tools again using cached images. Full DADA2, classification, phylogeny, diversity, and biological validity of the selected truncation pair are outside this test's coverage. See [tests/README.md](tests/README.md) for YAML/CLI precedence checks.

[GitHub Actions](.github/workflows/ci.yml) runs the Docker test and output verifier on pushes and pull requests. A `v*` tag publishes a GitHub release only after that tag's test passes. See [CHANGELOG.md](CHANGELOG.md) for release changes.

<a id="en-structure"></a>

### 10. Code structure

| Path | Purpose |
|---|---|
| [main.nf](main.nf) | Stage connections, branching, and input checks |
| [nextflow.config](nextflow.config) | Default parameters, containers, and runtime profiles |
| [nextflow_schema.json](nextflow_schema.json) | Parameter schema |
| [conf/base.config](conf/base.config) | Default task resources |
| [conf/modules.config](conf/modules.config) | Containers, publication paths and candidate-failure policy |
| [lib/PrimerConfig.groovy](lib/PrimerConfig.groovy) | Region primers and adapter resolution |
| [params/regions/](params/regions/) | Region preset YAML files |
| [params/examples/](params/examples/) | Custom primer and QC-only YAML examples |
| [lib/TrimmingConfig.groovy](lib/TrimmingConfig.groovy) | Trimming option validation |
| [bin/](bin/) | Cutadapt pass execution and summary generation |
| [params/trimming/](params/trimming/) | Default F/R candidate TSV files and design notes |
| [modules/](modules/) | Tool process commands |
| [subworkflows/qc.nf](subworkflows/qc.nf) | Shared raw/clean FastQC → MultiQC workflow |
| [subworkflows/trimm_optimal.nf](subworkflows/trimm_optimal.nf) | Candidate input, execution, ranking, and selection |
| [tests/](tests/) | Synthetic FASTQ, primer/container tests, output verifiers |

Module options are passed explicitly from `main.nf`; production modules and subworkflows do not read global `params`. Containers and publication paths are configured in `conf/modules.config`. See [module interfaces](modules/README.md) for reuse in other workflows. Existing user YAML options are unchanged.

Individual tools are called from `main.nf`. QC and truncation optimization remain subworkflows because they connect multiple steps.

<a id="en-trimming"></a>

### 11. Additional trimming and QC-only runs

To use the existing primer/read-through and quality filters without extra adapter or poly-tail removal:

```yaml
skip_cutadapt: false
discard_untrimmed: true
trim_poly_a: false
trim_poly_g: false
extra_adapters_f: []
extra_adapters_r: []
```

Set `discard_untrimmed: false` to retain pairs without primer/read-through matches. Primer removal, quality trimming and minimum-length filtering still run with either value. This block can be added to the existing analysis YAML; it does not change the end point of the run.

These options are available in the current checkout/development branch, not in the existing `v0.1.0` tag. Use [qc_only.yml](params/examples/qc_only.yml) as a starting point. QC-only runs require paired reads, but no metadata or classifier.

```yaml
reads: "/path/to/reads/*_{1,2}.fastq.gz"
outdir: "/path/to/qc_results"
region: V3V4
qc_only: true
# Adapters must be written in the orientation present in each read (5′→3′).
extra_adapters_f:
  - CTGTCTCTTATACACATCT
  - AGATGTGTATAAGAGACAG
extra_adapters_r:
  - CTGTCTCTTATACACATCT
  - AGATGTGTATAAGAGACAG
trim_poly_a: true
trim_poly_g: true
poly_min_length: 10
discard_untrimmed: false
```

The two example sequences are a Nextera transposase adapter and its reverse complement. Include only orientations relevant to your library; extra adapters are not reverse-complemented automatically. Each direction accepts multiple IUPAC DNA sequences. The default lists are empty, both poly options default to `false`, and `discard_untrimmed` defaults to `true`.

`discard_untrimmed: false` keeps pairs without primer/read-through matches; **it does not disable primer removal**. With `true`, the pair is discarded if either mate has no primer/read-through match, as before. A read-through match also counts as trimmed, so this is not a strict “5′ primer required” check. Additional adapters or tails do not rescue a pair discarded at this stage. `skip_cutadapt: true` skips all trimming.

Cutadapt runs in separate passes:

1. Existing Q20 (configurable `quality`) and primer/read-through removal, `-n 2`, minimum length and optional untrimmed-pair filtering.
2. If enabled, exact 3′ poly-A/G tail removal, including adjacent qualifying A/G runs.
3. If supplied, extra 3′ adapters, including Nextera. This pass can remove multiple matches.
4. If both features are enabled, another poly-tail pass removes tails exposed by extra-adapter removal.
5. A final minimum-length filter removes pairs made too short by the added trimming.

The first length filter preserves the previous default behavior and avoids processing pairs already too short; lengths cannot increase in later passes. With no extra adapters or poly trimming, only the original pass runs. Poly trimming uses non-internal Cutadapt 3′ adapters with zero errors and no indels, rather than `--nextseq-trim`. Only an exact consecutive terminal run of at least `poly_min_length` bases is removed; shorter runs and internal homopolymers remain. Biological terminal homopolymers can also meet this rule. See the [Cutadapt adapter guide](https://cutadapt.readthedocs.io/en/v5.2/guide.html#non-internal-5-and-3-adapters).

| End point | Option | Completed stages |
|---|---|---|
| Raw QC | `--stop_after raw_qc` | Raw FastQC and MultiQC |
| Trimming | `--stop_after cutadapt` | Raw QC, Cutadapt, trimming summary and its MultiQC report |
| Clean QC | `--stop_after clean_qc` or `--qc_only` | All of the above plus clean FastQC/MultiQC |
| Full analysis (default) | `--stop_after none` | Continue through QIIME 2 and downstream analysis |

`skip_qc: true` disables all FastQC/MultiQC tasks but still produces trimming JSON/TSV. It cannot be combined with `raw_qc`/`clean_qc`. `skip_cutadapt: true` cannot be combined with `cutadapt`/`clean_qc`; use `raw_qc` to inspect original reads only. `qc_only` cannot be combined with a different explicit end point.

To inspect raw QC and then continue, keep the same launch directory, work directory and trimming parameters:

```bash
nextflow run . -profile singularity -params-file analysis.yml \
  --qc_only false --stop_after raw_qc -work-dir /path/to/work
nextflow run . -profile singularity -params-file analysis.yml \
  --qc_only false --stop_after clean_qc -work-dir /path/to/work -resume
```

To continue into QIIME 2, provide `metadata` and `classifier`, set `--qc_only false --stop_after none`, and keep `-resume`. Changing trimming options reruns affected tasks. Reusing a result directory does not delete old outputs from stages omitted in the new run; use the current trace/log to identify stages actually executed.

`02_cutadapt_q20/` contains the final paired FASTQ, individual pass JSON, a combined command/log file, per-sample `*.trimming.summary.json` with parameters and commands, and `trimming_summary.tsv` for all samples. `trimming_summary_mqc.json` adds a table to `02_cutadapt_q20/multiqc/multiqc_report.html` (unless QC is skipped).

| TSV field | Unit and meaning |
|---|---|
| `input_pairs`, `retained_pairs` | Initial and final read pairs |
| `primer_matches`, `readthrough_matches` | Match events across both mates in the primer pass |
| `nextera_matches`, `other_extra_adapter_matches` | Match events for the two listed Nextera sequences versus other configured adapters |
| `poly_a_matches`, `poly_g_matches` | Exact tail-removal events across both mates and poly passes |
| `untrimmed_r1_reads`, `untrimmed_r2_reads` | Reads with no primer/read-through match before filtering, even when retained |
| `discarded_untrimmed_pairs`, `too_short_pairs` | Pairs removed by each filter, summed across passes |

Poly-A removal targets terminal runs, not internal homopolymers. A trailing non-A base prevents an exact A-tail match. Quality or read-through trimming can remove an A-rich end before the poly pass, leaving `poly_a_matches` at zero.

Match counts are not unique read/pair counts. If one adapter removes a downstream adapter or tail with it, the downstream sequence is not counted as another detected match. Pair counts satisfy `input_pairs = retained_pairs + discarded_untrimmed_pairs + too_short_pairs`.

Run the synthetic regression locally (CI uses Docker):

```bash
python3 tests/test_trimming.py --profile singularity --outdir results/trimming-test
```

Use a new test output directory. The test checks exact sequences and qualities for two samples, Nextera orientations, other adapters, read-through, adjacent/long/internal/short poly tails, untrimmed pairs, minimum length, MultiQC, early stops, invalid parameters and cached tasks across `-resume`. It does not use research data or run DADA2.

---

<a id="korean"></a>

## 한국어

### 목차

1. [빠른 실행](#ko-start)
2. [FASTQ와 샘플 정보](#ko-inputs)
3. [영역과 프라이머](#ko-primers)
4. [절단 길이 최적화](#ko-optimization)
5. [컨테이너와 캐시](#ko-containers)
6. [서버 자원 설정](#ko-resources)
7. [분석 설정과 코드 수정](#ko-options)
8. [결과 확인](#ko-outputs)
9. [테스트와 릴리스](#ko-tests)
10. [코드 구조](#ko-structure)
11. [추가 trimming과 QC 전용 실행](#ko-trimming)

<a id="ko-start"></a>

### 1. 빠른 실행

Illumina paired-end 16S rRNA의 V1V3 또는 V3V4 영역을 분석하는 파이프라인입니다. DADA2의 forward/reverse 절단 길이를 비교해 완료된 후보 중 가장 높은 순위의 조합을 선택하고, 해당 ASV로 후속 분석을 진행할 수 있습니다.

**분석 흐름:** FASTQ → FastQC/MultiQC → Cutadapt → FastQC/MultiQC → QIIME 2 가져오기 → DADA2(고정 길이 또는 후보 비교) → 분류 → MAFFT/FastTree → 다양성 분석.

Linux, Java 17+, Nextflow 26.04.2 이상(검증 버전: 26.04.2)과 Docker·Singularity·Apptainer 중 하나를 준비합니다. 입력 데이터는 FASTQ, metadata TSV, QIIME 2 amplicon 2025.7과 호환되는 classifier QZA입니다. SILVA, GG2 등 호환되는 분류기를 사용할 수 있습니다.

저장소 밖에 `analysis.yml`을 작성합니다. 파일명은 `nextflow_params.yml`로 정해도 됩니다. 아래 경로를 실제 경로로 바꾸세요.

```yaml
region: V3V4
reads: "/path/to/reads/*_{1,2}.fastq.gz"
metadata: "/path/to/metadata.tsv"
classifier: "/path/to/classifier.qza"
outdir: "/path/to/results"
run_label: "experiment01"
trimm_optimal: true
taxonomy_label: SILVA
taxonomy_confidence: 0.7
diversity_enabled: true
sampling_depth: 1000
```

저장소를 복제하지 않고 GitHub에서 릴리스 버전을 바로 실행할 수 있습니다.

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile singularity \
  -params-file /path/to/analysis.yml \
  -work-dir /path/to/work
```

실행 환경에 따라 `-profile docker`, `-profile singularity`(기본 실행 환경), `-profile apptainer`를 선택합니다. 이미지는 첫 실행 때 내려받고 이후 캐시를 재사용합니다. 분석 도구를 서버에 각각 설치할 필요는 없습니다.

저장소를 로컬에 복제했다면 다음과 같이 실행합니다.

```bash
nextflow run /path/to/amplicon_16S_qiime_nf/main.nf \
  -profile singularity -params-file /path/to/analysis.yml \
  -work-dir /path/to/work
```

V1V3 데이터는 `region: V1V3`로 설정합니다. 영역별 기본 후보를 사용하려면 `trimm_combinations`를 생략합니다. 위 예시는 최적화를 켠 설정이며, 코드의 `trimm_optimal` 기본값은 `false`입니다.

다양성 분석을 생략하려면 `diversity_enabled: false`로 설정합니다. `sampling_depth`는 샘플별 read 수에 맞게 조정하세요. 분석을 재개할 때는 같은 실행 폴더와 작업 폴더를 유지하고 `-resume`을 추가합니다. 입력·옵션·모듈 구조가 바뀌면 해당 작업이 다시 실행될 수 있습니다. 최신 개발 코드는 `-r v0.1.0` 대신 `-r main -latest`로 실행합니다.

<a id="ko-inputs"></a>

### 2. FASTQ와 샘플 정보

`reads`의 파일명 패턴으로 read 쌍을 묶고, `metadata`로 샘플 정보를 읽습니다. **FASTQ 경로를 나열한 CSV samplesheet(`--input`, `--samplesheet`)는 지원하지 않습니다.** QIIME 2 가져오기에 쓰는 manifest는 Cutadapt 처리 후 자동 생성합니다.

샘플마다 R1/R2 한 쌍을 준비합니다.

```text
reads/
  S01_1.fastq.gz
  S01_2.fastq.gz
  S02_1.fastq.gz
  S02_2.fastq.gz
```

위 파일은 `reads: "/path/to/reads/*_{1,2}.fastq.gz"`로 지정합니다. 추출되는 샘플 ID는 `S01`, `S02`입니다. 파일명이 `_R1_001.fastq.gz`, `_R2_001.fastq.gz`로 끝나면 `*_{R1,R2}_001.fastq.gz` 패턴을 사용하세요.

첫 열이 `sample-id`인 탭 구분 샘플 정보 파일을 준비합니다.

```tsv
sample-id	group
S01	control
S02	treatment
```

터미널에서 실제 탭 문자가 들어간 예제 파일을 만들려면 다음 명령을 사용합니다.

```bash
printf 'sample-id\tgroup\nS01\tcontrol\nS02\ttreatment\n' > metadata.tsv
```

분석 YAML의 `metadata`에 파일의 절대 경로를 지정합니다. 샘플 ID는 중복 없이 파일명에서 추출한 ID와 정확히 일치해야 합니다. `group`은 예시이므로 실험 조건·배치 등의 열로 바꾸거나 확장하세요. 열 구분자는 쉼표나 공백이 아닌 탭입니다. V1V3와 V3V4는 각 영역에 해당하는 reads와 metadata로 나누어 실행합니다.

<a id="ko-primers"></a>

### 3. 영역과 프라이머

기본 영역은 V1V3입니다.

| 영역 | Forward 프라이머 (5′→3′) | Reverse 프라이머 (5′→3′) |
|---|---|---|
| V1V3 | `AGAGTTTGATCCTGGCTCAG` | `ATTACCGCGGCTGCTGG` |
| V3V4 | `CCTACGGGNGGCWGCAG` | `GACTACHVGGGTATCTAATCC` |

`--region V3V4`로 선택하거나 [params/regions](params/regions)의 YAML을 사용합니다. 직접 지정하려면 실행 명령에 다음 옵션을 추가하세요.

```bash
--region V3V4 --primer_f 'CCTACGGGNGGCWGCAG' --primer_r 'GACTACHVGGGTATCTAATCC'
```

**적용 우선순위: 명령줄 > 외부 파라미터 YAML > 영역 기본값.** 한쪽 프라이머만 지정하면 나머지는 해당 영역의 기본값을 사용합니다. YAML 예시는 [custom_primers.yml](params/examples/custom_primers.yml)을 참고하세요.

프라이머는 소문자를 포함한 IUPAC DNA 문자를 허용합니다. Read-through adapter는 반대쪽 프라이머의 역상보 서열로 계산하며, `adapter_f` 또는 `adapter_r`를 직접 지정하면 그 값을 사용합니다. V1V3 기본 프라이머의 adapter는 `CCAGCAGCCGCGGTAAT`, `CTGAGCCAGGATCAAACTCT`입니다. 최종 서열은 실행 로그에 표시합니다.

<a id="ko-optimization"></a>

### 4. 절단 길이 최적화

`dada2_pooling_method`는 일반 DADA2와 모든 trimming 후보에 적용되며 `independent`(기본값) 또는 `pseudo`를 선택합니다. 다른 값은 workflow 시작 시 거부합니다. 파라미터 YAML에 지정하거나 `--dada2_pooling_method pseudo`로 전달하세요. QIIME 2 pseudo-pooling은 R DADA2의 완전한 `pool=TRUE`와 다릅니다. 기존 결과 경로는 유지하므로 실행을 비교할 때는 `outdir`를 다르게 지정하세요.

`trimm_optimal: true`로 후보 비교를 활성화합니다. 절단 길이는 Cutadapt 처리 후 QIIME 2로 가져온 reads에 적용합니다.

1. 후보 TSV의 forward/reverse 조합마다 DADA2를 독립 실행합니다.
2. 후보별 ASV를 분류하고 read 유지율과 species 분류 지표를 계산합니다.
3. 완료된 후보를 아래 기준으로 내림차순 정렬합니다. 앞 기준이 같으면 다음 기준으로 비교합니다.
4. 최상위 후보의 feature table과 대표 서열로 분류·계통수·다양성 분석을 진행합니다.

| 우선순위 | 지표 | 계산 |
|---|---|---|
| 1 | `species_input_yield` | Species로 분류된 reads / DADA2 입력 reads × 100 |
| 2 | `species_read_pct` | Species로 분류된 reads / 최종 feature table reads × 100 |
| 3 | `nonchimeric_pct` | Non-chimeric reads / DADA2 입력 reads × 100 |

여기서 “최적”은 지정한 분류기와 순위 기준에 따라 제공된 후보 중 가장 높은 순위라는 뜻입니다. 보편적인 생물학적 최적값을 보장하지 않습니다. 샘플별 순위가 아닌 전체 데이터의 합산 지표를 사용합니다. Species 판정에는 taxonomy의 `s__` 또는 `s:` 표기를 사용하며 unknown, uncultured 등의 이름은 제외합니다.

영역에 따라 [v1v3_10bp.tsv](params/trimming/v1v3_10bp.tsv)의 6개 후보 또는 [v3v4_10bp.tsv](params/trimming/v3v4_10bp.tsv)의 10개 후보를 사용합니다. 다른 후보를 쓰려면 `trimm_combinations: "/path/to/candidates.tsv"` 또는 `--trimm_combinations /path/to/candidates.tsv`로 지정합니다.

```tsv
name	trunc_len_f	trunc_len_r
F280_R270	280	270
F280_R260	280	260
F270_R250	270	250
```

열은 실제 탭으로 구분하고 후보 이름은 고유하게 지정합니다. 현재 코드는 **F > R**, 두 길이 모두 **10 nt 단위**인지 확인합니다. 기본 후보의 설계 기준은 V1V3에서 F + R − 20 ≥ 500 bp, V3V4에서 ≥ 460 bp입니다. 일반·최적화 DADA2 모두 최소 overlap은 12 bp입니다. 이 설계 기준이 검증된 생물학적 최대 길이를 뜻하지는 않으므로 실제 read 길이와 품질에 맞춰 조정하세요. 출처·길이 가정·overlap 계산은 [후보 설계 설명](params/trimming/README.md)에 정리되어 있습니다.

DADA2 또는 분류 작업이 실패한 후보는 비교에서 제외합니다. 성공한 후보가 없으면 선택 단계가 실패합니다. `optimal_min_sample_reads` 기본값은 10000이며, 기준보다 read가 적은 샘플 수와 merge가 0인 샘플 수는 보고용입니다. 후보 제외나 순위 계산에는 사용하지 않습니다.

고정 길이로 분석하려면 `trimm_optimal: false`와 `dada2_trunc_len_f`, `dada2_trunc_len_r`를 사용합니다. 두 길이의 기본값 `0`은 고정 길이 절단을 하지 않는다는 뜻입니다. 최적화 모드에서는 이 두 파라미터 대신 TSV의 길이를 사용합니다.

<a id="ko-containers"></a>

### 5. 컨테이너와 캐시

도구별 로컬 경로를 지정하지 않으면 아래 버전의 공개 이미지를 사용합니다.

| 도구 | 공개 이미지 | 로컬 SIF/IMG 파라미터 |
|---|---|---|
| FastQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | `fastqc_sif` |
| MultiQC | `quay.io/biocontainers/multiqc:1.27--pyhdfd78af_0` | `multiqc_sif` |
| Cutadapt | `quay.io/biocontainers/cutadapt:5.2--py311hc303176_2` | `cutadapt_sif` |
| QIIME 2 | `quay.io/qiime2/amplicon:2025.7` | `qiime_sif` |

FastQC와 MultiQC는 이미지와 설정을 각각 따로 사용합니다. 통합 QC 이미지와 공통 경로 설정은 더 이상 사용하지 않습니다. QIIME 2 분석 단계는 하나의 QIIME 2 이미지를 공유합니다.

보유한 이미지를 사용하려면 필요한 항목을 `analysis.yml` 또는 `nextflow_params.yml`에 추가합니다.

```yaml
fastqc_sif: "/absolute/path/to/fastqc.sif"
multiqc_sif: "/absolute/path/to/multiqc.img"
cutadapt_sif: "/absolute/path/to/cutadapt.sif"
qiime_sif: "/absolute/path/to/qiime2.sif"
```

로컬 SIF/IMG 파일은 `-profile singularity` 또는 `-profile apptainer`로 실행합니다. 해당 실행 환경에서 읽을 수 있는 이미지여야 하며, Docker는 이 파일을 사용할 수 없습니다. 경로를 생략한 도구는 공개 이미지를 자동으로 내려받습니다. 로컬 경로를 지정하면 해당 도구에만 우선 적용합니다.

첫 실행에는 이미지 저장소 접속과 충분한 디스크 공간이 필요합니다. 특히 QIIME 2 이미지의 용량을 고려하세요. Singularity/Apptainer 이미지 주소에는 `docker://`를 붙이고 Docker에서는 생략합니다. 자세한 내용은 [Nextflow 컨테이너 문서](https://docs.seqera.io/nextflow/container/singularity)를 참고하세요.

실행 전에 Singularity 이미지 캐시 경로를 지정하면 여러 분석에서 재사용할 수 있습니다.

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/container-cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
```

Apptainer는 `NXF_APPTAINER_CACHEDIR`를 사용합니다. Docker의 이미지 저장 공간은 Docker가 관리합니다. HPC에서는 모든 계산 노드가 접근할 수 있고 쓰기 가능한 캐시를 사용하세요. 이미지 캐시는 재다운로드를 줄이는 기능이며, 완료된 분석 결과를 재사용하는 `-resume`과는 별개입니다.

**현재 연구실 서버 설정.** 기존 분석 YAML에 아래 QC 경로를 넣습니다. reads·metadata·classifier 등 나머지 분석 설정은 함께 유지합니다.

```yaml
fastqc_sif: "/data/software/singularity/quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"
multiqc_sif: "/data/software/singularity/quay.io-biocontainers-multiqc-1.27--pyhdfd78af_0.img"
```

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile singularity -params-file /path/to/nextflow_params.yml
```

이 서버에는 동일한 QC 경로를 담은 `/data/software/singularity/amplicon_16S_qc.config`도 있습니다. 다음 명령으로 테스트할 수 있습니다.

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile test,singularity \
  -c /data/software/singularity/amplicon_16S_qc.config
```

위 절대 경로와 config 파일은 해당 서버 전용입니다. 외부 사용자는 자신의 로컬 경로를 지정하거나 경로 설정을 생략해 자동 다운로드를 사용하세요.

<a id="ko-resources"></a>

### 6. 서버 자원 설정

데이터 경로와 분석 파라미터는 외부 YAML에 둡니다. 실행 방식·queue·CPU·메모리·시간 제한은 외부 `site.config`로 설정합니다. [nextflow.config](nextflow.config)의 기본 실행 방식은 `local`입니다.

[conf/base.config](conf/base.config)에 정의된 작업 하나당 기본 자원은 다음과 같습니다.

| 라벨 | CPU 수 | 메모리 | 시간 |
|---|---|---|---|
| `process_low` | 8 | 16 GB | 12 h |
| `process_medium` | 8 | 24 GB | 24 h |
| `process_high` | 16 | 64 GB | 72 h |

최적화 후보는 동시에 실행될 수 있으므로 서버와 데이터 크기에 맞춰 자원을 조정하세요. 다음은 `site.config` 예시입니다.

```groovy
process {
    executor = 'local'
    maxForks = 1
    withLabel: process_low {
        cpus = 2
        memory = '4 GB'
    }
    withLabel: process_medium {
        cpus = 4
        memory = '8 GB'
    }
    withLabel: process_high {
        cpus = 8
        memory = '32 GB'
    }
}
```

실행 명령에 `-c /path/to/site.config`를 추가합니다. 위 수치는 문법 예시이며 모든 데이터에 적합한 값은 아닙니다. `maxForks`는 프로세스별 동시 작업 수를 제한하며, 전체 프로세스를 합친 동시 작업 수 제한은 아닙니다. HPC의 스케줄러와 queue는 기관 환경에 맞춰 지정합니다. 작업 폴더는 `-work-dir`로 설정합니다.

<a id="ko-options"></a>

### 7. 분석 설정과 코드 수정

단계 생략 옵션은 둘 다 기본값이 `false`입니다.

```yaml
skip_cutadapt: true
skip_qc: true
```

`skip_cutadapt`는 Cutadapt의 프라이머·adapter 제거, 품질 절단, 최소 길이 필터를 모두 생략하고 원본 read 쌍을 QIIME으로 가져옵니다. `discard_untrimmed`로 untrimmed-pair 제거를 선택합니다(기본값 `true`). 추가 trimming과 QC 전용 실행은 11절을 참고하세요. `skip_qc`는 trimming 집계 보고서를 포함한 모든 FastQC/MultiQC를 생략합니다. Cutadapt만 생략하면 원본 QC만 한 번 수행하고 처리 후 QC는 생략합니다. 샘플 ID는 기존 `reads` 패턴으로 추출하며 trimmed 파일명 접미사는 필요하지 않습니다. 생략한 단계의 새 결과는 생성하지 않으므로 이전 보고서와 혼동하지 않도록 새 결과 폴더를 사용하세요.

먼저 파라미터로 조정하고, 명령이나 로직을 바꿔야 할 때 도구 모듈을 수정합니다.

| 변경 목적 | 파라미터 또는 기본값 | 구현 위치 |
|---|---|---|
| FastQC / MultiQC | 개별 이미지 경로와 CPU 설정 | [fastqc.nf](modules/fastqc.nf), [multiqc.nf](modules/multiqc.nf) |
| 영역·프라이머·adapter | `region`, `primer_f`, `primer_r`, `adapter_f`, `adapter_r` | [PrimerConfig.groovy](lib/PrimerConfig.groovy) |
| Cutadapt 필터 | `quality: 20`, `min_length: 50` | [cutadapt.nf](modules/cutadapt.nf); primer pass는 `-n 2`; 기본값 `discard_untrimmed: true` |
| QIIME 가져오기 | `reads` 파일명 패턴 | [qiime_import.nf](modules/qiime_import.nf) |
| 고정 길이 DADA2 | `dada2_trunc_len_f: 0`, `dada2_trunc_len_r: 0`, `dada2_max_ee_f: 2`, `dada2_max_ee_r: 4` | [qiime_dada2.nf](modules/qiime_dada2.nf) |
| 후보 조건 검사 | `trimm_optimal`, `trimm_combinations` | [subworkflows/trimm_optimal.nf](subworkflows/trimm_optimal.nf), [params/trimming](params/trimming) |
| 후보 실행 | 동일한 `dada2_max_ee_*` 설정 | [modules/trimm_optimal.nf](modules/trimm_optimal.nf)의 `DADA2_TRIMM_SWEEP` |
| 순위·species 판정 | `optimal_min_sample_reads: 10000`은 보고 기준만 변경 | 같은 모듈의 `TAXONOMY_TRIMM_SWEEP`, `SELECT_TRIMM_OPTIMAL` |
| 분류 | `classifier`, `taxonomy_confidence: 0.7`, `taxonomy_label: AUTO` | [qiime_taxonomy.nf](modules/qiime_taxonomy.nf)와 `modules/trimm_optimal.nf`의 후보 분류 |
| MAFFT / FastTree | 별도 분석 파라미터 없음 | [qiime_phylogeny.nf](modules/qiime_phylogeny.nf) |
| 다양성·rarefaction | `diversity_enabled: true`, `sampling_depth: 1000`, `alpha_max_depth: 10000` | [qiime_diversity.nf](modules/qiime_diversity.nf) |

DADA2의 `trunc-q`나 chimera 처리처럼 명령에 고정된 옵션은 일반 모듈과 최적화 모듈을 모두 확인합니다. 후보 분류와 최종 분류 명령도 일관되게 유지하세요. 파라미터를 추가하면 [nextflow.config](nextflow.config), [nextflow_schema.json](nextflow_schema.json), 해당 모듈, README의 두 언어 설명을 함께 갱신합니다.

이미지 버전을 바꿀 때는 명령 옵션과 분류기 호환성을 확인합니다. 최종 분석에 사용한 이미지 버전과 분류기는 함께 기록하세요.

<a id="ko-outputs"></a>

### 8. 결과 확인

모든 결과는 `outdir` 아래에 저장합니다.

| 폴더 | 내용 |
|---|---|
| `01_raw_qc/`, `03_clean_qc/` | 필터링 전·후 FastQC/MultiQC |
| `02_cutadapt_q20/` | 처리된 FASTQ, pass별 로그·JSON, 집계 TSV 및 trimming MultiQC. `quality`를 바꿔도 폴더명은 고정 |
| `04_qiime2_import/` | Manifest, paired-end QZA, 품질 요약 |
| `05_dada2/` | 최적화를 끈 일반 DADA2 결과와 요약 |
| `05_trimm_optimal_dada2/` | 후보별 DADA2 결과 |
| `05_trimm_optimal_<LABEL>/` | 후보별 분류 결과와 지표 |
| `05_trimm_optimal_<LABEL>/selected/` | 전체 순위, 선택 길이, 선택된 QZA |
| `06_taxonomy/` | `taxonomy_<LABEL>.qza`, `taxa-bar-plots_<LABEL>.qzv` 등 |
| `07_phylogeny/`, `08_diversity/` | 계통수와 활성화한 경우 다양성 결과 |

최적화 후에는 다음 파일을 먼저 확인합니다.

- `all_parameter_results.tsv`: 성공한 후보의 순위와 평가 지표.
- `optimal_selection.tsv`, `optimal_truncation.txt`: 선택된 F/R 길이와 점수.
- `selected-table.qza`, `selected-rep-seqs.qza`, `selected-denoising-stats.qza`: 후속 분석에 전달한 산출물.

`taxonomy_label: AUTO`는 classifier 경로에서 SILVA·GG2·GTDB 등의 라벨을 추론합니다. 분석을 구분하려면 라벨을 직접 지정하세요. 최적화 모드에서는 일반 모드의 `05_dada2/table-summary.qzv`를 생성하지 않습니다. 다양성 분석 깊이는 선택된 table과 stats를 확인해 정합니다.

원시 reads, QIIME 산출물, 컨테이너 이미지, 작업 폴더, 참여자 정보는 Git에 올리지 않습니다.

<a id="ko-tests"></a>

### 9. 테스트와 릴리스

실제 연구 데이터·샘플 정보·분류기 없이 제공된 합성 데이터로 테스트할 수 있습니다.

```bash
nextflow run KitHubb/amplicon_16S_qiime_nf -r v0.1.0 \
  -profile test,docker
```

실행 환경에 따라 `test,singularity` 또는 `test,apptainer`로 바꿉니다. 로컬에 저장소가 있다면 결과 검증도 실행하세요.

```bash
nextflow run . -profile test,singularity
python3 tests/check_containers.py results/primer-test
```

정상 실행 시 이미지 4종으로 다음 작업 5개를 완료합니다.

| 작업 | 확인 내용 |
|---|---|
| Cutadapt | 프라이머 조건 검사와 실제 절단. 각 mate에 60개의 T로 된 서열 한 개가 남는지 확인 |
| FastQC | 처리된 reads의 보고서 생성 |
| MultiQC | FastQC 결과를 통합한 보고서 생성 |
| QIIME 가져오기 | Paired reads를 가져오고 demux 요약 생성 |
| QIIME 플러그인 검사 | 필요한 DADA2·feature-classifier·phylogeny·diversity·taxa 명령 로딩 |

보고서는 `results/primer-test/container_qc/`, QIIME 산출물은 `results/primer-test/04_qiime2_import/`, 플러그인 로그는 `results/primer-test/container_checks/qiime-container.txt`에 저장합니다.

공개 이미지 다운로드를 확인하려면 로컬 이미지 경로를 생략합니다. 빈 캐시에서 확인하려면 새로 만든 쓰기 가능한 이미지 캐시 폴더를 지정하세요. 캐시된 이미지로 도구를 다시 실행하려면 `-resume` 없이 재실행합니다. 전체 DADA2·분류·계통수·다양성 분석과 선택 길이의 생물학적 타당성은 이 테스트의 검증 범위가 아닙니다. YAML과 명령줄의 우선순위 검사는 [tests/README.md](tests/README.md)를 참고하세요.

[GitHub Actions](.github/workflows/ci.yml)는 push와 pull request마다 Docker 테스트와 결과 검증을 실행합니다. `v*` 태그의 테스트가 통과하면 GitHub 릴리스를 게시합니다. 버전별 변경은 [CHANGELOG.md](CHANGELOG.md)에서 확인할 수 있습니다.

<a id="ko-structure"></a>

### 10. 코드 구조

| 경로 | 역할 |
|---|---|
| [main.nf](main.nf) | 단계 연결·분기·입력 확인 |
| [nextflow.config](nextflow.config) | 기본 파라미터·컨테이너·실행 프로필 |
| [nextflow_schema.json](nextflow_schema.json) | 파라미터 명세 |
| [conf/base.config](conf/base.config) | 작업 자원 기본값 |
| [lib/PrimerConfig.groovy](lib/PrimerConfig.groovy) | 영역별 프라이머와 adapter 해석 |
| [params/regions/](params/regions/) | 영역 선택용 YAML |
| [params/examples/](params/examples/) | 사용자 프라이머 및 QC 전용 YAML 예시 |
| [lib/TrimmingConfig.groovy](lib/TrimmingConfig.groovy) | Trimming 옵션 검증 |
| [bin/](bin/) | Cutadapt pass 실행과 통계 집계 |
| [params/trimming/](params/trimming/) | 기본 F/R 후보 TSV와 설계 설명 |
| [modules/](modules/) | 도구별 실행 명령 |
| [subworkflows/qc.nf](subworkflows/qc.nf) | 처리 전·후 공통 FastQC → MultiQC 흐름 |
| [subworkflows/trimm_optimal.nf](subworkflows/trimm_optimal.nf) | 후보 입력·실행·순위·선택 연결 |
| [tests/](tests/) | 합성 FASTQ·프라이머/컨테이너 테스트·결과 검증 |

개별 도구는 `main.nf`에서 호출합니다. 여러 단계를 연결하는 QC와 절단 길이 최적화는 subworkflow로 구성합니다.

---

분석 설정은 `main.nf`에서 모듈 입력으로 명시적으로 전달하며, 분석 모듈과 서브워크플로는 전역 `params`를 직접 읽지 않습니다. 컨테이너·결과 경로·후보 실패 처리 정책은 `conf/modules.config`에서 관리합니다. 다른 workflow에서의 재사용 방법은 [모듈 입력 설명](modules/README.md)을 참고하세요. 기존 사용자 YAML 옵션은 유지합니다.

<a id="ko-trimming"></a>

### 11. 추가 trimming과 QC 전용 실행

추가 adapter와 poly-tail 제거 없이 기존 primer/read-through 제거와 품질 필터만 사용하려면 기존 분석 YAML에 다음을 지정합니다.

```yaml
skip_cutadapt: false
discard_untrimmed: true
trim_poly_a: false
trim_poly_g: false
extra_adapters_f: []
extra_adapters_r: []
```

Primer/read-through가 검출되지 않은 pair도 남기려면 `discard_untrimmed: false`로 바꿉니다. 어느 값을 쓰든 프라이머 제거, 품질 절단, 최소 길이 필터는 실행됩니다. 이 설정은 실행 종료 지점을 바꾸지 않습니다.

아래 옵션은 현재 체크아웃·개발 브랜치에 추가된 기능이며 기존 `v0.1.0` 태그에는 없습니다. [qc_only.yml](params/examples/qc_only.yml)을 복사해 시작하세요. QC 전용 실행에는 paired reads가 필요하며 metadata와 classifier는 필요하지 않습니다.

```yaml
reads: "/path/to/reads/*_{1,2}.fastq.gz"
outdir: "/path/to/qc_results"
region: V3V4
qc_only: true
# 각 read에 실제로 나타나는 방향(5′→3′)의 adapter 서열을 지정합니다.
extra_adapters_f:
  - CTGTCTCTTATACACATCT
  - AGATGTGTATAAGAGACAG
extra_adapters_r:
  - CTGTCTCTTATACACATCT
  - AGATGTGTATAAGAGACAG
trim_poly_a: true
trim_poly_g: true
poly_min_length: 10
discard_untrimmed: false
```

예제의 두 서열은 Nextera transposase adapter와 그 역상보 서열입니다. 라이브러리에 필요한 방향만 지정하세요. 추가 adapter의 역상보 서열은 자동 생성하지 않습니다. 방향별로 여러 IUPAC DNA 서열을 입력할 수 있습니다. 기본값은 빈 adapter 목록, `trim_poly_a: false`, `trim_poly_g: false`, `discard_untrimmed: true`입니다.

`discard_untrimmed: false`는 primer/read-through가 검출되지 않은 pair도 보존하지만 **프라이머 제거 자체를 끄지는 않습니다**. `true`이면 기존처럼 어느 한 mate라도 primer/read-through가 검출되지 않은 pair를 제거합니다. Read-through 검출도 trimmed로 취급하므로 “5′ 프라이머 필수” 검사와는 다릅니다. 이 단계에서 제거된 pair는 추가 adapter나 tail이 있어도 복구되지 않습니다. `skip_cutadapt: true`는 모든 trimming을 생략합니다.

Cutadapt를 다음과 같이 나누어 실행합니다.

1. 기존 Q20(`quality`로 변경 가능), primer/read-through 제거, `-n 2`, 최소 길이 필터 및 선택적 untrimmed-pair 제거.
2. 활성화한 경우 정확히 일치하는 3′ poly-A/G tail 제거. 연속된 A/G 구간도 각각 기준을 만족하면 제거합니다.
3. 지정한 경우 Nextera를 포함한 추가 3′ adapter 제거. 여러 번의 일치를 처리할 수 있습니다.
4. Poly-tail과 추가 adapter를 모두 사용하면 adapter 제거 후 드러난 tail을 한 번 더 처리합니다.
5. 추가 trimming으로 짧아진 pair에 최종 최소 길이 필터를 적용합니다.

첫 길이 필터는 기존 기본 동작을 유지하면서 이미 너무 짧은 pair의 후속 처리를 줄입니다. 이후 단계에서 서열 길이가 늘어나지는 않습니다. 추가 adapter와 poly 제거를 모두 끄면 기존 첫 pass만 실행합니다. Poly 제거에는 `--nextseq-trim` 대신 Cutadapt의 말단 전용 3′ adapter를 사용하며 오류와 indel을 허용하지 않습니다. 말단에서 같은 염기가 `poly_min_length` 이상 정확히 연속될 때만 제거하므로 더 짧은 tail과 내부 homopolymer는 유지합니다. 생물학적 서열의 말단 homopolymer도 이 조건을 만족하면 제거될 수 있습니다. [Cutadapt adapter 설명](https://cutadapt.readthedocs.io/en/v5.2/guide.html#non-internal-5-and-3-adapters)을 참고하세요.

| 종료 지점 | 옵션 | 수행 단계 |
|---|---|---|
| 원본 QC | `--stop_after raw_qc` | 원본 FastQC·MultiQC |
| Trimming | `--stop_after cutadapt` | 원본 QC, Cutadapt, trimming 집계 및 해당 MultiQC 보고서 |
| 처리 후 QC | `--stop_after clean_qc` 또는 `--qc_only` | 위 단계와 처리 후 FastQC·MultiQC |
| 전체 분석(기본값) | `--stop_after none` | QIIME 2와 후속 분석까지 계속 실행 |

`skip_qc: true`는 모든 FastQC·MultiQC를 생략하지만 trimming JSON·TSV는 생성합니다. `raw_qc`·`clean_qc`와는 함께 사용할 수 없습니다. `skip_cutadapt: true`는 `cutadapt`·`clean_qc`와 함께 사용할 수 없으며 원본만 확인하려면 `raw_qc`를 사용하세요. `qc_only`와 다른 명시적 종료 지점도 함께 지정할 수 없습니다.

원본 QC를 확인한 뒤 계속하려면 실행 위치, work directory, trimming 파라미터를 유지합니다.

```bash
nextflow run . -profile singularity -params-file analysis.yml \
  --qc_only false --stop_after raw_qc -work-dir /path/to/work
nextflow run . -profile singularity -params-file analysis.yml \
  --qc_only false --stop_after clean_qc -work-dir /path/to/work -resume
```

QIIME 2까지 계속하려면 `metadata`와 `classifier`를 제공하고 `--qc_only false --stop_after none`으로 변경한 뒤 `-resume`을 유지하세요. Trimming 옵션을 바꾸면 관련 task가 다시 실행됩니다. 같은 결과 폴더를 재사용해도 이번에 생략한 단계의 이전 결과가 삭제되지는 않습니다. 실제 수행 단계는 현재 trace·로그로 확인하세요.

`02_cutadapt_q20/`에는 최종 paired FASTQ, pass별 JSON, 통합 명령·로그, 파라미터와 실행 명령이 담긴 샘플별 `*.trimming.summary.json`, 전체 샘플의 `trimming_summary.tsv`가 생성됩니다. `trimming_summary_mqc.json`은 `02_cutadapt_q20/multiqc/multiqc_report.html`에 통계 표를 표시합니다(QC 생략 시 HTML 제외).

| TSV 항목 | 단위 및 의미 |
|---|---|
| `input_pairs`, `retained_pairs` | 최초 입력 및 최종 보존 read pair 수 |
| `primer_matches`, `readthrough_matches` | 첫 pass에서 양쪽 mate의 해당 서열 검출·제거 횟수 |
| `nextera_matches`, `other_extra_adapter_matches` | 예제의 두 Nextera 서열과 그 외 추가 adapter의 검출·제거 횟수 |
| `poly_a_matches`, `poly_g_matches` | 양쪽 mate와 모든 poly pass의 정확한 tail 제거 횟수 |
| `untrimmed_r1_reads`, `untrimmed_r2_reads` | 필터 전 primer/read-through가 검출되지 않은 read 수. 보존한 read도 포함 |
| `discarded_untrimmed_pairs`, `too_short_pairs` | 각 필터에서 제거된 pair 수를 pass 전체에 걸쳐 합산 |

Poly-A 제거는 말단의 연속 A에 적용하며 내부 homopolymer는 유지합니다. 끝에 A가 아닌 염기가 붙으면 정확한 A-tail 조건에 맞지 않습니다. 앞선 품질·read-through 처리에서 A 구간이 잘리면 `poly_a_matches`는 0일 수 있습니다.

검출 횟수는 고유 read·pair 수와 다릅니다. 한 adapter를 자르면서 그 뒤의 다른 adapter나 tail까지 함께 제거되면, 뒤쪽 서열을 별도의 검출로 중복 집계하지 않습니다. Pair 수는 `input_pairs = retained_pairs + discarded_untrimmed_pairs + too_short_pairs`를 만족합니다.

합성 데이터 회귀 테스트는 다음과 같이 실행합니다(CI에서는 Docker 사용).

```bash
python3 tests/test_trimming.py --profile singularity --outdir results/trimming-test
```

새 테스트 결과 폴더를 지정하세요. 두 샘플의 정확한 서열·품질값, Nextera 방향, 기타 adapter, read-through, 인접·긴·내부·짧은 poly 구간, untrimmed pair, 최소 길이 필터, MultiQC, 조기 종료, 잘못된 파라미터, `-resume`의 task 캐시 재사용을 확인합니다. 연구 데이터나 DADA2 실행은 포함하지 않습니다.
