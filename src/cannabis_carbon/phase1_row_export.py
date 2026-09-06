"""Lossless collection-row export for large Phase 1 evidence reports."""
import hashlib
import json
from pathlib import Path


def encode(report, sha):
    collections = {k: v for k, v in report.items() if isinstance(v, list) and v and all(isinstance(x, dict) for x in v)}
    metadata = {'report_fields': {k: v for k, v in report.items() if k not in collections},
                'collections': {k: len(v) for k, v in collections.items()}}
    entries = [('metadata', 'report', metadata)]
    for k, values in collections.items():
        for i, value in enumerate(values):
            identity = value.get('cannabisdb_id') or value.get('id') or value.get('compound_id') or 'record'
            entries.append((k, f'{i:08d}:{identity}', value))
    return [{'record_kind': kind, 'record_id': identity, 'record_json': json.dumps(value, separators=(',', ':')),
             'report_sha256': sha} for kind, identity, value in entries]


def decode(rows):
    if len({(r['record_kind'], r['record_id']) for r in rows}) != len(rows):
        raise ValueError('Duplicate exported record')
    if len({r['report_sha256'] for r in rows}) != 1:
        raise ValueError('Mixed report digests')
    metadata = [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'metadata' and r['record_id'] == 'report']
    if len(metadata) != 1:
        raise ValueError('Exactly one metadata record required')
    meta = metadata[0]; result = dict(meta['report_fields'])
    if {r['record_kind'] for r in rows} != {'metadata', *meta['collections']}:
        raise ValueError('Unexpected collection')
    for k, count in meta['collections'].items():
        records = sorted((r for r in rows if r['record_kind'] == k), key=lambda r: r['record_id'])
        if len(records) != count or any(not r['record_id'].startswith(f'{i:08d}:') for i, r in enumerate(records)):
            raise ValueError('Missing or invalid collection ordinal')
        result[k] = [json.loads(r['record_json']) for r in records]
    return result


def run(name):
    path = Path('data/reports/phase1-' + name + '.json')
    payload = path.read_bytes(); report = json.loads(payload)
    rows = encode(report, hashlib.sha256(payload).hexdigest())
    if decode(list(reversed(rows))) != report:
        raise ValueError('Lossless unordered export replay failed')
    Path('data/derived/phase1-' + name + '.ndjson').write_text(''.join(json.dumps(r, separators=(',', ':')) + '\n' for r in rows))
    Path('docs/data/' + name + '.json').write_bytes(payload)
    print(json.dumps({'records': len(rows), 'report_sha256': rows[0]['report_sha256']}))
