#!/usr/bin/env python3
"""Check real module commands with test doubles, without pipeline params or QIIME analysis."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
base = Path(sys.argv[1] if len(sys.argv)>1 else 'results/module-reuse-test').resolve()
base.mkdir(parents=True, exist_ok=False)
fixture=base/'fixture.txt'
fixture.write_text('TEST DOUBLE - NOT RESEARCH DATA\n')
metadata=base/'metadata.tsv'
metadata.write_text('sample-id\tgroup\nS1\ttest\n')
candidates=base/'candidates.tsv'
candidates.write_text('name\ttrunc_len_f\ttrunc_len_r\nF250_R200\t250\t200\nF240_R190\t240\t190\n')
config=base/'standalone.config'
config.write_text("process { executor = 'local'; cpus = 1; shell = ['/bin/bash', '-euo', 'pipefail'] }\n")
command=['nextflow','-C',str(config),'run',str(ROOT/'tests/reuse.nf'),'--fixture',str(fixture),'--metadata',str(metadata),'--candidates',str(candidates),'-work-dir',str(base/'work'),'-with-trace',str(base/'trace.tsv'),'-ansi-log','false']
runtime=base/'runtime'
runtime.mkdir()
(runtime/'python').symlink_to(sys.executable)
env=dict(os.environ,NXF_DISABLE_CHECK_LATEST='true',PATH=str(runtime)+os.pathsep+str(ROOT/'tests/mocks')+os.pathsep+os.environ['PATH'])
with (base/'run.log').open('w') as log:
    result=subprocess.run(command,cwd=str(base),env=env,stdout=log,stderr=subprocess.STDOUT)
assert result.returncode==0,(base/'run.log').read_text(encoding='utf-8')
with (base/'trace.tsv').open() as handle:
    rows=list(csv.DictReader(handle,delimiter='\t'))
assert len(rows)==10 and all(row['status']=='COMPLETED' for row in rows),rows
commands=[]
for path in (base/'work').glob('*/*/mock_commands.jsonl'):
    commands.extend(json.loads(line) for line in path.read_text().splitlines())
def option(args,key): return args[args.index(key)+1]
dada=[a for a in commands if a[:2]==['dada2','denoise-paired']]
assert len(dada)==4
assert sorted(option(a,'--p-pooling-method') for a in dada)==['independent','pseudo','pseudo','pseudo']
assert all(option(a,'--p-max-ee-f')=='3' and option(a,'--p-max-ee-r')=='5' for a in dada)
assert sorted(option(a,'--p-trunc-len-f') for a in dada)==['240','250','250','250']
tax=[a for a in commands if a[:2]==['feature-classifier','classify-sklearn']]
assert len(tax)==3 and all(option(a,'--p-confidence')=='0.83' for a in tax)
core=next(a for a in commands if a[:2]==['diversity','core-metrics-phylogenetic'])
rare=next(a for a in commands if a[:2]==['diversity','alpha-rarefaction'])
assert option(core,'--p-sampling-depth')=='37' and option(rare,'--p-max-depth')=='83'
selection=list((base/'work').glob('*/*/optimal_selection.tsv'))
assert len(selection)==1
with selection[0].open() as handle:
    winner=next(csv.DictReader(handle,delimiter='\t'))
assert winner['samples_below_min_reads']=='1'
for directory in ['modules','subworkflows']:
    for path in (ROOT/directory).glob('*.nf'):
        assert 'params.' not in path.read_text(),path
print('PASS: standalone module inputs, pooling/truncation/EE/confidence/depth, sweep selection; 10 tasks. QIIME commands mocked, no biological validation.')

# Exercise the actual pipeline publish rules with mocked tool commands.
config.write_text(config.read_text() + "\nparams.outdir = '" + str(base/'published') + "'\nincludeConfig '" + str(ROOT/'conf/modules.config') + "'\n")
command[command.index('-work-dir')+1]=str(base/'publish-work')
command[command.index('-with-trace')+1]=str(base/'publish-trace.tsv')
with (base/'publish.log').open('w') as log:
    result=subprocess.run(command,cwd=str(base),env=env,stdout=log,stderr=subprocess.STDOUT)
assert result.returncode==0,(base/'publish.log').read_text(encoding='utf-8')
for relative in ['05_dada2/table.qza', '06_taxonomy/taxonomy_TEST.qza', '07_phylogeny/rooted-tree.qza',
                 '08_diversity/alpha-rarefaction.qzv', '05_trimm_optimal_dada2/F250_R200.table.qza',
                 '05_trimm_optimal_TEST/F250_R200.metrics.tsv', '05_trimm_optimal_TEST/selected/optimal_selection.tsv']:
    assert (base/'published'/relative).is_file(),relative
print('PASS: pipeline publishing paths preserved (mock outputs only)')
