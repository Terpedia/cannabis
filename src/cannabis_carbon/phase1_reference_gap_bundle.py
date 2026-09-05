"""Lossless Terpedia export of reference-gap audits and experimental leads."""
import hashlib
import json
from pathlib import Path

SOURCES = (
    'data/reports/phase1-missing-reference-review.json',
    'data/reports/phase1-nonplant-reference-review.json',
    'data/curation/glucuronolactone-literature-review.json',
    'data/curation/glucuronolactone-fulltext-review.json',
    'data/curation/biopterin-reference-review.json',
    'data/reports/phase1-biopterin-lead-references.json',
    'data/reports/phase1-biopterin-lead-search.json',
    'data/reports/phase1-pteridine-family-references.json',
    'data/reports/phase1-pteridine-family-search.json',
)


def build(sources=SOURCES):
    documents = []
    for name in sources:
        payload = Path(name).read_bytes()
        document = json.loads(payload)
        for path, digest in document.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                raise ValueError(f'Stale document lineage: {name}: {path}')
        documents.append({'source_path': name, 'source_sha256': hashlib.sha256(payload).hexdigest(),
                          'document': document})
    return {'schema': 'cannabis-reference-gap-bundle-v1', 'documents': documents,
            'summary': {'source_documents': len(documents), 'new_exact_enzyme_assignments': 0,
                        'candidate_model_changed': False},
            'claim_boundary': 'Lossless reference-gap reviews, generic-substrate searches and experimental hypotheses. '
            'Negative searches do not establish enzyme absence. Generic substrate annotations and in-vitro chemical '
            'interconversion do not close exact Cannabis reaction gaps. Historical balanced network and atom tracing unchanged.'}


def run(sources=SOURCES, name='reference-gap-bundle'):
    report = build(sources)
    path = Path('data/reports/phase1-' + name + '.json')
    path.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    Path('docs/data/' + name + '.json').write_bytes(path.read_bytes())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows = [{'record_kind': 'source_document', 'record_id': d['source_path'],
             'record_json': json.dumps(d, separators=(',', ':')), 'report_sha256': digest}
            for d in report['documents']]
    rows.append({'record_kind': 'metadata', 'record_id': 'report',
                 'record_json': json.dumps({k: v for k, v in report.items() if k != 'documents'}, separators=(',', ':')),
                 'report_sha256': digest})
    Path('data/derived/phase1-' + name + '.ndjson').write_text(
        ''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows))
    print(json.dumps({'records': len(rows), 'sha256': digest}))


if __name__ == '__main__':
    run()
