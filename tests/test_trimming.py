#!/usr/bin/env python3
"""Synthetic-only end-to-end regression, early stops and Nextflow resume."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
F = 'ACGTACGTACGT'
R = 'CGATCGATCGAT'
N = 'CTGTCTCTTATACACATCT'
RC = 'AGATGTGTATAAGAGACAG'
OTHER = 'GATCTAGTACGATCGTAC'
INSERT = 'TTCTTC' * 10


def read_fastq(path):
    with gzip.open(str(path), 'rt') as handle:
        lines = handle.read().splitlines()
    assert len(lines) % 4 == 0
    records = {}
    for index in range(0, len(lines), 4):
        name, sequence, plus, qualities = lines[index:index + 4]
        assert name.startswith('@') and plus == '+' and len(sequence) == len(qualities)
        assert name[1:] not in records
        assert qualities == 'I' * len(sequence)
        records[name[1:]] = sequence
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', default='docker', choices=['docker', 'singularity', 'apptainer'])
    parser.add_argument('--outdir', default='results/trimming-test')
    args = parser.parse_args()
    base = Path(args.outdir).resolve()
    base.mkdir(parents=True, exist_ok=False)
    fixtures = base / 'fixtures'
    fixtures.mkdir()
    # name, mate1 insert/tail, mate2 insert/tail, primer present in R1/R2, expected R1/R2
    cases = [
        ('nextera_poly', INSERT + N + 'A'*14 + 'G'*12, INSERT + RC + 'G'*13, True, True, INSERT, INSERT),
        ('no_primer', INSERT + N + 'A'*10, INSERT + RC + 'G'*10, False, False, INSERT, INSERT),
        ('one_primer', INSERT, INSERT, True, False, INSERT, INSERT),
        ('short_final', 'TTC' + N + 'A'*10, INSERT + RC, True, True, None, None),
        ('short_primer', 'TTC', 'TTC', True, True, None, None),
        ('below_minimum', INSERT + 'A'*9, INSERT + 'G'*9, True, True, INSERT + 'A'*9, INSERT + 'G'*9),
        ('at_minimum', INSERT + 'A'*10, INSERT + 'G'*10, True, True, INSERT, INSERT),
        ('long_tail', INSERT + 'A'*150, INSERT + 'G'*150, True, True, INSERT, INSERT),
        ('internal_tail', INSERT + 'A'*12 + 'TTC', INSERT + 'G'*12 + 'TTC', True, True, INSERT + 'A'*12 + 'TTC', INSERT + 'G'*12 + 'TTC'),
        ('revealed_tail', INSERT + 'A'*10 + OTHER, INSERT + 'G'*10 + OTHER, True, True, INSERT, INSERT),
        ('readthrough', INSERT + 'ATCGATCGATCG', INSERT + 'ACGTACGTACGT', True, True, INSERT, INSERT),
        ('unchanged', INSERT, INSERT, True, True, INSERT, INSERT),
    ]
    for sample in ('sample1', 'sample2'):
        for mate in (1, 2):
            with gzip.open(str(fixtures / (sample + '_' + str(mate) + '.fastq.gz')), 'wt') as handle:
                for name, r1, r2, p1, p2, _, _ in cases:
                    seq = ((F if p1 else '') + r1) if mate == 1 else ((R if p2 else '') + r2)
                    handle.write('@' + name + '\n' + seq + '\n+\n' + 'I'*len(seq) + '\n')
    config = dict(reads=str(fixtures / '*_{1,2}.fastq.gz'), primer_f=F, primer_r=R,
                  extra_adapters_f=[N, OTHER], extra_adapters_r=[RC, OTHER], trim_poly_a=True,
                  trim_poly_g=True, poly_min_length=10, discard_untrimmed=False, min_length=30)
    resource = base / 'resources.config'
    resource.write_text("process { cpus = 2; memory = '4 GB'; withLabel: process_low { cpus = 2; memory = '4 GB' }; withLabel: process_medium { cpus = 2; memory = '4 GB' } }\n")
    def run(name, updates, resume=False, fail=None):
        params = dict(config, outdir=str(base / 'output'))
        params.update(updates)
        path = base / (name + '.json')
        path.write_text(json.dumps(params))
        trace = base / (name + '.trace.tsv')
        command = ['nextflow', 'run', str(ROOT), '-profile', args.profile, '-params-file', str(path),
                   '-c', str(resource), '-work-dir', str(base / 'work'), '-ansi-log', 'false', '-with-trace', str(trace)]
        if resume:
            command.append('-resume')
        with open(str(base / (name + '.log')), 'w') as log:
            result = subprocess.run(command, cwd=str(ROOT), stdout=log, stderr=subprocess.STDOUT,
                                    env=dict(os.environ, NXF_DISABLE_CHECK_LATEST='true'))
        output = (base / (name + '.log')).read_text(encoding='utf-8')
        if fail:
            assert result.returncode != 0 and fail in output, output
            print('PASS validation: ' + name, flush=True)
            return []
        assert result.returncode == 0, output
        with open(str(trace)) as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert rows and all(row['status'] in ('COMPLETED', 'CACHED') for row in rows), rows
        assert all('QIIME' not in row['name'] for row in rows), rows
        print('PASS Nextflow: ' + name, flush=True)
        return rows

    raw = run('raw', {'stop_after': 'raw_qc'})
    assert all('RAW_QC:' in row['name'] for row in raw)
    assert not (base / 'output/02_cutadapt_q20').exists()
    trimmed = run('cutadapt', {'stop_after': 'cutadapt'}, resume=True)
    assert all(row['status'] == 'CACHED' for row in trimmed if 'RAW_QC:' in row['name'])
    assert not (base / 'output/03_clean_qc').exists()
    clean = run('clean', {'qc_only': True}, resume=True)
    assert all(row['status'] == 'CACHED' for row in clean if 'CLEAN_QC:' not in row['name'])
    assert any('CLEAN_QC:' in row['name'] for row in clean)
    repeated = run('repeat', {'stop_after': 'clean_qc'}, resume=True)
    assert all(row['status'] == 'CACHED' for row in repeated), repeated

    def check_results(directory, discard):
        for sample in ('sample1', 'sample2'):
            for mate in (1, 2):
                expected = {row[0]: row[5 if mate == 1 else 6] for row in cases
                            if row[5] is not None and (not discard or (row[3] and row[4]))}
                actual = read_fastq(directory / (sample + '.R' + str(mate) + '.trimmed.fastq.gz'))
                assert actual == expected, (sample, mate, expected, actual)
            report = json.loads((directory / (sample + '.trimming.summary.json')).read_text(encoding='utf-8'))
            counts = report['metrics']
            assert counts['input_pairs'] == len(cases)
            assert counts['retained_pairs'] == (8 if discard else 10), counts
            assert counts['too_short_pairs'] == 2, counts
            assert counts['discarded_untrimmed_pairs'] == (2 if discard else 0), counts
            assert counts['untrimmed_r1_reads'] == 1 and counts['untrimmed_r2_reads'] == 2, counts
            expected_matches = dict(primer_matches=21, readthrough_matches=2, nextera_matches=4 if discard else 6,
                                    other_extra_adapter_matches=2, poly_a_matches=5 if discard else 6, poly_g_matches=5 if discard else 6)
            for key, value in expected_matches.items():
                assert counts[key] == value, (key, counts)
            assert counts['input_pairs'] == counts['retained_pairs'] + counts['too_short_pairs'] + counts['discarded_untrimmed_pairs']
            assert report['passes'] == ['primer', 'poly_before_extra', 'extra', 'poly_after_extra', 'final']
        with open(str(directory / 'trimming_summary.tsv')) as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert [row['sample'] for row in rows] == ['sample1', 'sample2']
        custom = json.loads((directory / 'trimming_summary_mqc.json').read_text(encoding='utf-8'))
        assert set(custom['data']) == {'sample1', 'sample2'}

    check_results(base / 'output/02_cutadapt_q20', False)
    multiqc = base / 'output/02_cutadapt_q20/multiqc/multiqc_report_data'
    assert (base / 'output/02_cutadapt_q20/multiqc/multiqc_report.html').is_file()
    multiqc_data = json.loads((multiqc / 'multiqc_data.json').read_text(encoding='utf-8'))
    saved = multiqc_data['report_saved_raw_data']['multiqc_trimming_summary_table']
    assert set(saved) == {'sample1', 'sample2'} and saved['sample1']['retained_pairs'] == 10, saved
    run('discard', {'stop_after': 'cutadapt', 'skip_qc': True, 'discard_untrimmed': True, 'outdir': str(base / 'discard')})
    check_results(base / 'discard/02_cutadapt_q20', True)
    assert not (base / 'discard/02_cutadapt_q20/multiqc').exists()
    edge = base / 'edge_fixtures'
    edge.mkdir()
    def edge_case(name, reads, updates, expected, passes):
        for mate, seq in enumerate(reads, 1):
            with gzip.open(str(edge / ('edge_' + str(mate) + '.fastq.gz')), 'wt') as handle:
                handle.write('@edge\n' + seq + '\n+\n' + 'I'*len(seq) + '\n')
        options = dict(reads=str(edge / '*_{1,2}.fastq.gz'), stop_after='cutadapt', skip_qc=True,
                       outdir=str(base / name), extra_adapters_f=[], extra_adapters_r=[], trim_poly_a=False, trim_poly_g=False)
        options.update(updates)
        run(name, options)
        directory = base / name / '02_cutadapt_q20'
        for mate, seq in enumerate(expected, 1):
            assert read_fastq(directory / ('edge.R' + str(mate) + '.trimmed.fastq.gz')) == ({} if seq is None else {'edge': seq})
        report = json.loads((directory / 'edge.trimming.summary.json').read_text(encoding='utf-8'))
        assert report['passes'] == passes, report
    edge_case('poly_a_only', [F + INSERT + 'A'*10, R + INSERT + 'G'*10], {'trim_poly_a': True},
              [INSERT, INSERT + 'G'*10], ['primer', 'poly_before_extra', 'final'])
    edge_case('poly_g_only', [F + INSERT + 'A'*10, R + INSERT + 'G'*10], {'trim_poly_g': True},
              [INSERT + 'A'*10, INSERT], ['primer', 'poly_before_extra', 'final'])
    edge_case('extra_f_only', [F + INSERT + N, R + INSERT + RC], {'extra_adapters_f': [N]},
              [INSERT, INSERT + RC], ['primer', 'extra', 'final'])
    edge_case('empty_output', [INSERT, INSERT], {'discard_untrimmed': True, 'trim_poly_a': True},
              [None, None], ['primer', 'poly_before_extra', 'final'])
    for name, updates, message in [
        ('bad_stop', {'stop_after': 'invalid'}, 'stop_after must be one of'),
        ('bad_poly_length', {'poly_min_length': 0}, 'poly_min_length must be an integer'),
        ('bad_adapters', {'extra_adapters_f': N}, 'extra_adapters_f must be a YAML list'),
        ('bad_boolean', {'discard_untrimmed': 'maybe'}, 'discard_untrimmed must be true or false'),
        ('conflict_qc', {'stop_after': 'raw_qc', 'skip_qc': True}, 'requires skip_qc=false'),
        ('conflict_cutadapt', {'qc_only': True, 'skip_cutadapt': True}, 'requires skip_cutadapt=false'),
    ]:
        run(name, updates, fail=message)
    print('PASS exact sequences, metrics, MultiQC, early stops, discard and resume (synthetic data only)', flush=True)


if __name__ == '__main__':
    main()
