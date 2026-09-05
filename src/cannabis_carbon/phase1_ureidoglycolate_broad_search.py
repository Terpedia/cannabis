"""EC/name-only reference leads, never automatic cofactor-specific assignments."""
import copy
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from .phase1_new_protein_search import run as screen, export_table

QUERIES = {
    'ec154': 'ec:1.1.1.154 AND fragment:false',
    'name': 'protein_name:"ureidoglycolate dehydrogenase" AND fragment:false',
}
RAW = Path('data/raw/phase1-ureidoglycolate-broad-search')
BOUNDARY = ('EC/name-only discovery with no automatic exact-reaction assignment. '
    'Both reviewed and unreviewed records retained with source annotations. Cofactor specificity, '
    'reference assay identity, Cannabis activity, reaction direction and in-vivo flux are unverified. '
    'Any passing homology is an experimental lead only, ineligible for model integration. Atom tracing deferred.')


def lookup(name, query, raw=RAW):
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
        'query': query, 'format': 'tsv',
        'fields': 'accession,reviewed,organism_name,protein_name,ec,rhea'})
    snapshot, cache = raw / (name + '.tsv'), raw / (name + '-lookup.json')
    if not cache.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        records = list(csv.DictReader(io.StringIO(payload.decode()), delimiter='\t'))
        if not payload.startswith(b'Entry\tReviewed\t') or any(not r.get('Entry') for r in records):
            raise ValueError('Invalid annotation response')
        snapshot.write_bytes(payload)
        cache.write_text(json.dumps({'url': url, 'query': query, 'snapshot': str(snapshot),
            'sha256': hashlib.sha256(payload).hexdigest(), 'records': len(records),
            'retrieved_at': datetime.now(timezone.utc).isoformat()}) + '\n')
    result = json.loads(cache.read_text())
    if result['url'] != url or hashlib.sha256(snapshot.read_bytes()).hexdigest() != result['sha256']:
        raise ValueError('Annotation lookup mismatch')
    return result


def build(lookups):
    source = Path('data/reports/phase1-ureide-gap-references.json')
    review = Path('data/curation/ureidoglycolate-cofactor-review.json')
    cofactor = json.loads(review.read_text())
    parent = json.loads(source.read_text())
    row = copy.deepcopy(next(r for r in parent['rows'] if r['reaction_id'] == cofactor['nadp_reaction_id']))
    row['prior_exact_reference_matches'] = row.pop('reference_matches')
    row['related_nad_reaction_id'] = cofactor['nad_reaction_id']
    row['model_eligible'] = False
    refs = {}
    paths = [source, review]
    for item in lookups:
        path = Path(item['snapshot'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('Changed annotation snapshot')
        paths.append(path)
        records = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'))
        if len(records) != item['records']:
            raise ValueError('Annotation count mismatch')
        for record in records:
            match = refs.setdefault(record['Entry'], {'accession': record['Entry'],
                'model_eligible': False, 'exact_reaction_annotation_match': False,
                'match_type': 'EC-or-name-only; cofactor-and-substrate-specificity-unverified',
                'source_records': []})
            match['source_records'].append({'lookup_url': item['url'], 'query': item['query'], 'record': record})
    row['reference_matches'] = [refs[a] for a in sorted(refs)]
    row['lookup_status'] = 'broad-annotation-leads-only'
    row['next_step'] = 'Screen full proteome; verify substrate/cofactor activity before exact reaction assignment.'
    return {'schema': 'cannabis-ureidoglycolate-broad-discovery-v1', 'rows': [row],
        'model_eligible': False, 'lookups': lookups,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'distinct_reference_leads': len(refs), 'exact_enzyme_assignments': 0},
        'claim_boundary': BOUNDARY}


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    report = build([lookup(k, q) for k, q in QUERIES.items()])
    source = Path('data/reports/phase1-ureidoglycolate-broad-references.json')
    source.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)
    if report['rows'][0]['reference_matches']:
        output = Path('data/reports/phase1-ureidoglycolate-broad-search.json')
        screen(source, RAW, output, evidence_class='cofactor-unverified-broad-annotation-homology-lead',
            additional_blockers=('not-eligible-for-exact-reaction-model', 'EC-or-name-only-reference; cofactor-unverified'),
            claim_boundary=BOUNDARY)
        result = json.loads(output.read_text())
        result['model_eligible'] = False
        for row in result['rows']:
            row['model_eligible'] = False
        output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
        export_table(output, Path('data/derived/phase1-ureidoglycolate-broad-search.ndjson'))


if __name__ == '__main__':
    run()
