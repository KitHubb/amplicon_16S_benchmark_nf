#!/usr/bin/env python3
"""Combine per-sample trimming metrics into TSV and MultiQC custom content."""
import csv
import json
import sys

rows = []
for path in sys.argv[1:]:
    with open(path) as handle:
        rows.append(json.load(handle))
rows.sort(key=lambda row: row['sample'])
assert rows and len({row['sample'] for row in rows}) == len(rows), 'Empty or duplicate sample summaries'
columns = list(rows[0]['metrics'])
with open('trimming_summary.tsv', 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['sample'] + columns)
    for row in rows:
        writer.writerow([row['sample']] + [row['metrics'][key] for key in columns])
report = {'id': 'trimming_summary', 'section_name': 'Multi-pass trimming', 'plot_type': 'table',
          'description': 'Input/retained/filtered counts are read pairs. Untrimmed R1/R2 counts are reads without a primer/read-through match before filtering. Adapter/poly counts are match events across both mates and passes, not unique pairs; downstream sequence removed with an adapter is not counted again.',
          'pconfig': {'id': 'trimming_summary_table', 'title': 'Trimming per sample'},
          'headers': {key: {'title': key.replace('_', ' '), 'format': '{:,.0f}'} for key in columns},
          'data': {row['sample']: row['metrics'] for row in rows}}
with open('trimming_summary_mqc.json', 'w') as handle:
    json.dump(report, handle, indent=2)
