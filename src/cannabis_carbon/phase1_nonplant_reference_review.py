"""Broaden the two missing-reference families without claiming enzyme absence."""
import copy
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from .phase1_new_references import attach


def run():
    source = Path('data/reports/phase1-missing-reference-review.json')
    parent = json.loads(source.read_text())
    rows = copy.deepcopy(parent['rows'])
    masters = sorted({f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()})
    query = '(' + ' OR '.join('cc_catalytic_activity:"' + m.lower() + '"' for m in masters) + ') AND NOT taxonomy_id:33090 AND reviewed:false AND fragment:false'
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': query,
        'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
    raw = Path('data/raw/phase1-nonplant-reference-review'); raw.mkdir(parents=True, exist_ok=True)
    snapshot, cache = raw / 'nonplant-unreviewed.tsv', raw / 'lookup.json'
    if not cache.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        snapshot.write_bytes(payload)
        lookup = {'url': url, 'query': query, 'requested_master_ids': masters, 'status': 'retrieved',
            'snapshot': str(snapshot), 'sha256': hashlib.sha256(payload).hexdigest(),
            'retrieved_at': datetime.now(timezone.utc).isoformat(), 'review_status': 'unreviewed',
            'taxonomy_scope': 'Excludes Viridiplantae (33090)'}
        cache.write_text(json.dumps(lookup) + '\n')
    lookup = json.loads(cache.read_text())
    if lookup['url'] != url:
        raise ValueError('Cached request mismatch')
    proteins = attach(rows, [lookup])
    for row in rows:
        row['lookup_status'] = 'unreviewed-nonplant-references-found' if row['reference_matches'] else 'no-unreviewed-nonplant-reference-returned'
        row['next_step'] = 'Screen exact-family sequences with unreviewed-reference uncertainty.' if row['reference_matches'] else 'Review primary literature and exact enzyme-independent chemistry; do not infer enzyme absence.'
        for match in row['reference_matches']:
            match['review_status'] = 'unreviewed'
    report = {'schema': 'cannabis-carbon.phase1-nonplant-reference-review.v1', 'rows': rows,
        'reference_proteins': [{**proteins[a], 'review_status': 'unreviewed'} for a in sorted(proteins)],
        'lookups': [lookup], 'source_sha256': {str(source): hashlib.sha256(source.read_bytes()).hexdigest()},
        'summary': {'equation_gaps': len(rows), 'unreviewed_nonplant_reference_records': len(proteins),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows)},
        'claim_boundary': 'Exact Rhea-family query of unreviewed nonfragment sequences outside Viridiplantae. No name-only joins, no Cannabis activity claim, no negative inference of enzyme absence. Primary literature may identify enzyme-independent chemistry; no model change is authorized by missing annotation alone.'}
    Path('data/reports/phase1-nonplant-reference-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
