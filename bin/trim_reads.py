#!/usr/bin/env python3
"""Run separate Cutadapt passes and record explicit match/read/pair statistics."""
import argparse
import base64
import gzip
import json
import os
import shlex
import subprocess

NEXTERA = {'CTGTCTCTTATACACATCT', 'AGATGTGTATAAGAGACAG'}


def max_length(paths):
    longest = 0
    for path in paths:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for index, line in enumerate(handle):
                if index % 4 == 1:
                    longest = max(longest, len(line.strip()))
    return longest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--cores', required=True)
    parser.add_argument('reads', nargs=2)
    args = parser.parse_args()
    config = json.loads(base64.b64decode(args.config))
    sample = args.sample
    current = args.reads
    reports = []
    temporary = []
    commands = []

    def run(stage, options, final=False):
        nonlocal current
        outputs = ([sample + '.R1.trimmed.fastq.gz', sample + '.R2.trimmed.fastq.gz']
                   if final else ['intermediate.' + stage + '.R1.fastq.gz', 'intermediate.' + stage + '.R2.fastq.gz'])
        report = sample + ('.cutadapt.json' if stage == 'primer' else '.' + stage + '.cutadapt.json')
        command = ['cutadapt', '--cores', args.cores] + options + ['--json', report, '-o', outputs[0], '-p', outputs[1]] + current
        commands.append(command)
        with open(sample + '.cutadapt.log', 'a') as log:
            log.write('\n### ' + stage + '\n' + ' '.join(shlex.quote(x) for x in command) + '\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        with open(report) as handle:
            reports.append((stage, json.load(handle)))
        current = outputs
        if not final:
            temporary.extend(outputs)

    extra = bool(config['extra_adapters_f'] or config['extra_adapters_r'])
    poly = config['trim_poly_a'] or config['trim_poly_g']
    advanced = extra or poly
    p = config['primers']
    options = ['-g', 'primer_f=' + p['primer_f'], '-a', 'readthrough_f=' + p['adapter_f'],
               '-G', 'primer_r=' + p['primer_r'], '-A', 'readthrough_r=' + p['adapter_r'],
               '-n', '2', '-q', str(config['quality']), '-Q', str(config['quality']),
               '--minimum-length', str(config['min_length'])]
    if config['discard_untrimmed']:
        options.append('--discard-untrimmed')
    run('primer', options, final=not advanced)
    if advanced:
        longest = max_length(current)
        def poly_pass(stage):
            options = ['-e', '0', '--no-indels', '-O', str(config['poly_min_length']),
                       '-n', str(max(1, longest // config['poly_min_length'] + 1))]
            for base in ('A', 'G'):
                if config['trim_poly_' + base.lower()]:
                    # X forbids internal matches; a longer-than-read adapter allows
                    # exact suffixes of any length >= poly_min_length, without a cap.
                    adapter = 'poly_' + base.lower() + '=' + base + '{' + str(max(longest + 1, config['poly_min_length'])) + '}X'
                    options.extend(['-a', adapter, '-A', adapter])
            run(stage, options)
        if poly:
            poly_pass('poly_before_extra')
        if extra:
            options = ['-n', str(max(1, longest // 3 + 1))]
            for direction, flag in (('f', '-a'), ('r', '-A')):
                for index, sequence in enumerate(config['extra_adapters_' + direction]):
                    name = ('nextera' if sequence in NEXTERA else 'extra') + '_' + direction + '_' + str(index + 1)
                    options.extend([flag, name + '=' + sequence])
            run('extra', options)
        if poly and extra:
            poly_pass('poly_after_extra')
        run('final', ['--minimum-length', str(config['min_length'])], final=True)

    first = reports[0][1]['read_counts']
    metrics = {'input_pairs': first['input'], 'primer_matches': 0, 'readthrough_matches': 0,
               'nextera_matches': 0, 'other_extra_adapter_matches': 0, 'poly_a_matches': 0, 'poly_g_matches': 0,
               'untrimmed_r1_reads': first['input'] - first['read1_with_adapter'],
               'untrimmed_r2_reads': first['input'] - first['read2_with_adapter'],
               'discarded_untrimmed_pairs': 0, 'too_short_pairs': 0,
               'retained_pairs': reports[-1][1]['read_counts']['output']}
    for stage, report in reports:
        counts = report['read_counts']['filtered']
        metrics['discarded_untrimmed_pairs'] += counts.get('discard_untrimmed') or 0
        metrics['too_short_pairs'] += counts.get('too_short') or 0
        for key in ('adapters_read1', 'adapters_read2'):
            for adapter in report.get(key) or []:
                name = adapter['name']
                metric = next((key for prefix, key in [('primer_', 'primer_matches'), ('readthrough_', 'readthrough_matches'),
                              ('nextera_', 'nextera_matches'), ('extra_', 'other_extra_adapter_matches'),
                              ('poly_a', 'poly_a_matches'), ('poly_g', 'poly_g_matches')] if name.startswith(prefix)), None)
                if metric:
                    metrics[metric] += adapter['total_matches']
    assert metrics['input_pairs'] == metrics['retained_pairs'] + metrics['too_short_pairs'] + metrics['discarded_untrimmed_pairs']
    with open(sample + '.trimming.summary.json', 'w') as handle:
        json.dump({'sample': sample, 'metrics': metrics, 'parameters': config,
                   'commands': commands, 'passes': [name for name, _ in reports]}, handle, indent=2)
    for path in temporary:
        os.remove(path)


if __name__ == '__main__':
    main()
