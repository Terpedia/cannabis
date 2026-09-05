"""Validate and export the curated occurrence review without upgrading claims."""
import hashlib
import json
from pathlib import Path

from .phase1_scope import write_rows


def validate(review, table, branches):
    rows = review['rows']
    expected = {d['cannabisdb_id'] for d in branches['decisions']}
    if len(rows) != len(expected) or {r['cannabisdb_id'] for r in rows} != expected:
        raise ValueError('Incomplete or duplicate priority accession inventory')
    originals = {r['accession']: r for r in table['rows']}
    for row in rows:
        original = originals[row['cannabisdb_id']]
        if row['name'] != original['name'] or row['source_reference_pmids'] != [r['pubmed_id'] for r in original['references']]:
            raise ValueError('Source reference inventory changed')
        if not row.get('occurrence_status') or not row.get('next_action'):
            raise ValueError('Missing evidence status or next action')
        for observation in row.get('observations', []):
            if observation['snapshot'] not in review['source_sha256']:
                raise ValueError('Untracked observation source')
            if observation['pdf_page_one_based'] < 1 or not observation.get('url'):
                raise ValueError('Missing observation locator')
    return [('metadata', 'review', {k: v for k, v in review.items() if k != 'rows'}),
            *(('occurrence', r['cannabisdb_id'], r) for r in rows)]


def run():
    path = Path('data/curation/priority-occurrence-review.json')
    review = json.loads(path.read_text())
    for source, sha in review['source_sha256'].items():
        if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
            raise ValueError('Source checksum mismatch')
    records = validate(review, json.loads(Path('data/terpedia/cannabisdb-compounds.json').read_text()),
                       json.loads(Path('data/reports/phase1-identity-branches.json').read_text()))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    count = write_rows(records, sha, Path('data/derived/priority-occurrence-review.ndjson'))
    print(json.dumps({'rows': count, 'report_sha256': sha}))


if __name__ == '__main__':
    run()
