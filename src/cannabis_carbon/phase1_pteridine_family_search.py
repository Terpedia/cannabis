"""Broader reviewed pteridine references; never an exact biopterin assignment."""
import copy
import csv
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .phase1_new_protein_search import run as screen


def run():
    parent_path = Path('data/reports/phase1-biopterin-lead-references.json')
    parent = json.loads(parent_path.read_text())
    raw = Path('data/raw/phase1-pteridine-family-search')
    raw.mkdir(parents=True, exist_ok=True)
    query = '(cc_catalytic_activity:"rhea:17865" OR cc_catalytic_activity:"rhea:17869") AND reviewed:true AND fragment:false'
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
        'query': query, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
    snapshot, cache = raw / 'reviewed-masters.tsv', raw / 'master-lookup.json'
    if not cache.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        snapshot.write_bytes(payload)
        cache.write_text(json.dumps({'url': url, 'query': query, 'snapshot': str(snapshot),
            'sha256': hashlib.sha256(payload).hexdigest(),
            'retrieved_at': datetime.now(timezone.utc).isoformat()}) + '\n')
    lookup = json.loads(cache.read_text())
    if lookup['url'] != url or hashlib.sha256(snapshot.read_bytes()).hexdigest() != lookup['sha256']:
        raise ValueError('Reference lookup cache mismatch')
    row = copy.deepcopy(parent['rows'][0])
    row['reference_matches'] = []
    allowed = {f'RHEA:{n}' for n in range(17865, 17873)}
    for record in csv.DictReader(io.StringIO(snapshot.read_text()), delimiter='\t'):
        ids = set(re.findall(r'RHEA:\d+', record['Rhea ID']))
        matched = sorted(ids & allowed)
        if not matched:
            raise ValueError('Returned reference lacks requested generic reaction annotation')
        row['reference_matches'].append({
            'accession': record['Entry'], 'review_status': 'reviewed',
            'source_record': record, 'matched_generic_rhea_ids': matched,
            'match_type': 'generic-substrate-rhea-family-not-exact-biopterin',
            'exact_reaction_annotation_match': False, 'model_eligible': False})
    report = {'schema': 'cannabis-pteridine-generic-reference-discovery-v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (parent_path, snapshot, cache)},
        'lookups': [lookup], 'rows': [row], 'claim_boundary': parent['claim_boundary']}
    source = Path('data/reports/phase1-pteridine-family-references.json')
    source.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(f'Reviewed generic references: {len(row["reference_matches"])}', flush=True)
    if row['reference_matches']:
        screen(source, raw, Path('data/reports/phase1-pteridine-family-search.json'),
               evidence_class='generic-substrate-specificity-unverified-homology-lead',
               additional_blockers=('not-eligible-for-exact-reaction-model',
                                    'generic-reference-substrate-not-exact-biopterin-assignment'),
               claim_boundary=parent['claim_boundary'])


if __name__ == '__main__':
    run()
