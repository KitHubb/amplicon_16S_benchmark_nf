"""Check actual outputs from all four container images after the test workflow."""
import argparse
import gzip
import csv
import runpy
import sys
from pathlib import Path
from zipfile import ZipFile

parser = argparse.ArgumentParser()
parser.add_argument('outdir')
parser.add_argument('--skip-cutadapt', action='store_true')
parser.add_argument('--skip-qc', action='store_true')
args = parser.parse_args()
root = Path(args.outdir)
if args.skip_cutadapt:
    assert not (root / '02_cutadapt_q20').exists(), 'Skipped Cutadapt produced outputs (use a fresh outdir)'
else:
    runpy.run_path(str(Path(__file__).with_name('check_trimmed.py')), run_name='__main__')
if args.skip_qc:
    assert not (root / 'container_qc').exists(), 'Skipped QC produced outputs (use a fresh outdir)'
else:
    for mate in (1, 2):
        name = f'r{mate}_fastqc' if args.skip_cutadapt else f'synthetic.R{mate}.trimmed_fastqc'
        with ZipFile(root / 'container_qc' / 'fastqc' / f'{name}.zip') as archive:
            assert any(p.endswith('/fastqc_data.txt') for p in archive.namelist())
    report = root / 'container_qc' / 'multiqc' / 'multiqc_report.html'
    assert report.stat().st_size > 0, f'Missing or empty report: {report}'
with (root / '04_qiime2_import' / 'manifest.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert [r['sample-id'] for r in rows] == ['synthetic']
if args.skip_cutadapt:
    with ZipFile(root / '04_qiime2_import' / 'demux-paired.qza') as archive:
        for mate in (1, 2):
            member = next(n for n in archive.namelist() if '/data/' in n and f'_R{mate}_' in n and n.endswith('.fastq.gz'))
            actual = gzip.decompress(archive.read(member)).decode().splitlines()
            expected = (Path(__file__).parent / 'fixtures' / f'r{mate}.fastq').read_text().splitlines()
            assert actual == expected, 'Raw reads changed while Cutadapt was skipped'
for name in ('demux-paired.qza', 'demux-summary.qzv'):
    with ZipFile(root / '04_qiime2_import' / name) as archive:
        assert archive.testzip() is None, f'Corrupt QIIME artifact: {name}'
        assert any(p.endswith('/metadata.yaml') for p in archive.namelist())
assert 'PASS: required QIIME actions are available' in (root / 'container_checks' / 'qiime-container.txt').read_text()
print('PASS: enabled stage outputs, import sample IDs, and skip behavior verified')
