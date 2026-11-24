# Module interfaces / 모듈 입력

## English

Production modules and subworkflows do not read global `params`. The calling workflow validates user parameters and passes values explicitly. `conf/modules.config` defines this pipeline's containers, publication paths and candidate-failure policy; `conf/base.config` defines resources. A separate consumer must supply its own runtime/container configuration and publication rules.

| Process | Inputs, in order |
|---|---|
| `FASTQC` | `(meta, r1, r2)`, `stage` |
| `MULTIQC` | report files, `stage` |
| `CUTADAPT` | `(meta, r1, r2)`, trimming options map, primers map, path to `trim_reads.py` |
| `TRIMMING_SUMMARY` | per-sample summary files, path to `summarize_trimming.py` |
| `QIIME_IMPORT` | run label, `(sample_ids, forward_reads, reverse_reads)` |
| `QIIME_DADA2` | DADA2 options map, run label, demux artifact |
| `QIIME_TAXONOMY` | taxonomy options map, run label, representative sequences, table, metadata, classifier, taxonomy label |
| `DADA2_TRIMM_SWEEP` | DADA2 options map, `(candidate_id, trunc_f, trunc_r, demux)` |
| `TAXONOMY_TRIMM_SWEEP` | taxonomy options map, `(candidate_id, trunc_f, trunc_r, table, repseq, stats)`, classifier, taxonomy label |
| `SELECT_TRIMM_OPTIMAL` | taxonomy label, candidate artifact/metrics files |

Option map keys:

- Trimming: `extra_adapters_f`, `extra_adapters_r`, `trim_poly_a`, `trim_poly_g`, `poly_min_length`, `discard_untrimmed`, `quality`, `min_length`.
- Primers: `primer_f`, `primer_r`, `adapter_f`, `adapter_r`; the pipeline also records `region`.
- DADA2: `dada2_pooling_method`, `dada2_max_ee_f`, `dada2_max_ee_r`; standard DADA2 additionally requires `dada2_trunc_len_f`, `dada2_trunc_len_r`. Sweep truncation lengths come from its tuple input.
- Taxonomy: `taxonomy_confidence`; candidate evaluation additionally requires `optimal_min_sample_reads`.

`QC` takes a reads channel and stage. `TRIMM_OPTIMAL_WORKFLOW` takes demux, classifier, taxonomy label, a candidate TSV **Path**, DADA2 options and taxonomy options. Use value channels for shared demux/classifier/label inputs to the sweep. Output channel names and filenames are unchanged.

`tests/reuse.nf` is an independent consumer. `python3 tests/test_reuse.py` runs it without the pipeline configuration, checks the real generated commands with test doubles, then verifies the pipeline's publication paths. It does not perform QIIME biological analysis. Tests using real containers remain separate. Changes to process interfaces can invalidate caches from older pipeline revisions; resume between runs of the new revision is covered by the trimming regression.

## 한국어

분석 모듈과 서브워크플로는 전역 `params`를 직접 읽지 않습니다. 호출하는 workflow가 사용자 설정을 검증하고 입력값으로 전달합니다. 이 파이프라인의 컨테이너·결과 저장 경로·후보 실패 처리 정책은 `conf/modules.config`, 자원 설정은 `conf/base.config`에서 관리합니다. 다른 파이프라인에서 가져다 쓰려면 실행환경·컨테이너와 결과 저장 설정을 직접 지정해야 합니다.

| Process | 순서대로 전달할 입력 |
|---|---|
| `FASTQC` | `(meta, r1, r2)`, `stage` |
| `MULTIQC` | 보고서 파일들, `stage` |
| `CUTADAPT` | `(meta, r1, r2)`, trimming 설정 map, 프라이머 map, `trim_reads.py` 경로 |
| `TRIMMING_SUMMARY` | 샘플별 통계 파일들, `summarize_trimming.py` 경로 |
| `QIIME_IMPORT` | 실행 이름, `(sample_ids, forward_reads, reverse_reads)` |
| `QIIME_DADA2` | DADA2 설정 map, 실행 이름, demux artifact |
| `QIIME_TAXONOMY` | 분류 설정 map, 실행 이름, 대표 서열, table, metadata, classifier, 분류 이름 |
| `DADA2_TRIMM_SWEEP` | DADA2 설정 map, `(candidate_id, trunc_f, trunc_r, demux)` |
| `TAXONOMY_TRIMM_SWEEP` | 분류 설정 map, `(candidate_id, trunc_f, trunc_r, table, repseq, stats)`, classifier, 분류 이름 |
| `SELECT_TRIMM_OPTIMAL` | 분류 이름, 후보 artifact·지표 파일들 |

설정 map에 필요한 항목:

- Trimming: `extra_adapters_f`, `extra_adapters_r`, `trim_poly_a`, `trim_poly_g`, `poly_min_length`, `discard_untrimmed`, `quality`, `min_length`.
- 프라이머: `primer_f`, `primer_r`, `adapter_f`, `adapter_r`. 현재 파이프라인은 `region`도 기록합니다.
- DADA2: `dada2_pooling_method`, `dada2_max_ee_f`, `dada2_max_ee_r`. 일반 DADA2는 `dada2_trunc_len_f`, `dada2_trunc_len_r`도 필요합니다. 후보 비교에서는 tuple의 절단 길이를 사용합니다.
- 분류: `taxonomy_confidence`. 후보 평가에는 `optimal_min_sample_reads`도 필요합니다.

`QC`는 reads channel과 stage를 받습니다. `TRIMM_OPTIMAL_WORKFLOW`는 demux, classifier, 분류 이름, 후보 TSV의 **Path**, DADA2 설정, 분류 설정을 받습니다. Sweep에 공통으로 전달할 demux·classifier·분류 이름은 value channel을 사용하세요. 출력 channel 이름과 파일명은 기존과 같습니다.

`tests/reuse.nf`는 별도의 호출 workflow 예제입니다. `python3 tests/test_reuse.py`는 파이프라인 config 없이 이 workflow를 실행하여 실제 생성 명령을 테스트 대역으로 검사하고, 이어서 기존 결과 저장 경로도 검증합니다. QIIME 생물학적 분석을 수행하는 테스트는 아닙니다. 실제 컨테이너를 사용하는 테스트는 별도로 유지합니다. Process 입력 변경으로 이전 코드 버전의 캐시가 무효화될 수 있으며, 새 버전 내 실행 간 resume은 trimming 회귀 테스트에서 확인합니다.
