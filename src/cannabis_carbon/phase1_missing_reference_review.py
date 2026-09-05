"""Explicitly unreviewed plant-reference follow-up for missing reviewed families."""
import copy
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from .phase1_new_references import attach
from .phase1_new_protein_search import run as screen, export_table


def prepare(discovery, search):
    absent = {r['reaction_id']: r for r in search['rows'] if r['search_status'] == 'no-reference-sequence' and not r['reference_matches']}
    rows = [copy.deepcopy(r) for r in discovery['rows'] if r['reaction_id'] in absent]
    for row in rows:
        row['prior_reviewed_search'] = copy.deepcopy(absent[row['reaction_id']])
    return rows


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-remaining-gap-references', 'phase1-remaining-gap-search')]
    discovery, search = [json.loads(p.read_text()) for p in paths]
    if hashlib.sha256(paths[0].read_bytes()).hexdigest() != search['source_discovery_sha256']:
        raise ValueError('Prior discovery lineage changed')
    rows = prepare(discovery, search)
    masters = sorted({f['RHEA_ID_MASTER'] for row in rows for f in row['rhea_families'].values()})
    if not masters:
        raise ValueError('No missing reviewed reaction families')
    query = '(' + ' OR '.join('cc_catalytic_activity:"' + m.lower() + '"' for m in masters) + ') AND taxonomy_id:33090 AND reviewed:false AND fragment:false'
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': query,
        'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
    raw = Path('data/raw/phase1-missing-reference-review'); raw.mkdir(parents=True, exist_ok=True)
    snapshot, cache = raw / 'plant-unreviewed.tsv', raw / 'lookup.json'
    if not cache.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        snapshot.write_bytes(payload)
        lookup = {'url': url, 'query': query, 'requested_master_ids': masters, 'status': 'retrieved',
            'snapshot': str(snapshot), 'sha256': hashlib.sha256(payload).hexdigest(),
            'retrieved_at': datetime.now(timezone.utc).isoformat(), 'review_status': 'unreviewed', 'taxonomy_scope': 'Viridiplantae (33090)'}
        cache.write_text(json.dumps(lookup) + '\n')
    lookup = json.loads(cache.read_text())
    if lookup['url'] != url:
        raise ValueError('Cached request mismatch')
    proteins = attach(rows, [lookup])
    for row in rows:
        row['lookup_status'] = 'unreviewed-plant-references-found' if row['reference_matches'] else 'no-unreviewed-plant-reference-returned'
        for ref in row['reference_matches']:
            ref['review_status'] = 'unreviewed'
        row['next_step'] = 'Screen exact-family references, retaining unreviewed annotation uncertainty.' if row['reference_matches'] else 'Review nonplant references and primary reaction literature; missing plant annotation is not enzyme absence.'
    report = {'schema': 'cannabis-carbon.phase1-missing-reference-review.v1', 'rows': rows,
        'reference_proteins': [{**proteins[a], 'review_status': 'unreviewed'} for a in sorted(proteins)],
        'lookups': [lookup], 'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'equation_gaps': len(rows), 'unreviewed_plant_reference_records': len(proteins),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows)},
        'claim_boundary': 'Unreviewed nonfragment Viridiplantae annotations only, joined through exact published Rhea families. Prior negative reviewed searches retained. Not characterized reference activity, Cannabis activity, physiological direction, or proof of absence. Atom tracing deferred.'}
    output = Path('data/reports/phase1-missing-reference-review.json')
    output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)
    if any(row['reference_matches'] for row in rows):
        result = Path('data/reports/phase1-missing-reference-search.json')
        screen(output, raw, result, evidence_class='unreviewed-plant-reference-homology-candidate',
            additional_blockers=('unreviewed-reference-activity-unverified',), claim_boundary=report['claim_boundary'])
        export_table(result, Path('data/derived/phase1-missing-reference-search.ndjson'))


if __name__ == '__main__':
    run()
