"""Independent plant CPR homology screen; no inferred P450 partnerships."""
import csv
import hashlib
import json
from pathlib import Path
from .phase1_ureidoglycolate_broad_search import lookup
from .phase1_new_protein_search import run as screen, export_table

BOUNDARY = ('Plant CPR class homology only. Reviewed reference status does not establish '
    'Cannabis activity or FNS-II compatibility. Distinct protein-bound carrier states '
    'remain unresolved; no exact model assignment, free-flavin substitution or atom-tracing claim.')


def run():
    raw = Path('data/raw/cpr-search')
    raw.mkdir(parents=True, exist_ok=True)
    query = 'ec:1.6.2.4 AND taxonomy_id:33090 AND reviewed:true AND fragment:false'
    item = lookup('reviewed-plant-cpr', query, raw=raw)
    snapshot = Path(item['snapshot'])
    records = list(csv.DictReader(snapshot.open(), delimiter='\t'))
    if not records:
        raise ValueError('No reviewed plant CPR references recovered')
    review_path = Path('data/curation/cpr-fnsii-carrier-interface-review.json')
    audit_path = Path('data/reports/phase1-cpr-annotation-audit.json')
    review = json.loads(review_path.read_text())
    row = {'reaction_id': review['cpr_annotation_master'], 'model_eligible': False,
        'target_ids': [], 'priority_target_ids': [], 'hypothesis_ids': [],
        'carrier_interface_review': review, 'compatible_fnsii_partners': [],
        'reference_matches': [{'accession': r['Entry'], 'model_eligible': False,
            'exact_reaction_annotation_match': False, 'match_type': 'reviewed-plant-CPR-class-only',
            'source_records': [{'record': r, 'query': query, 'lookup_url': item['url']}]}
            for r in records]}
    source = Path('data/reports/phase1-cpr-references.json')
    paths = [snapshot, raw / 'reviewed-plant-cpr-lookup.json', review_path, audit_path]
    discovery = {'schema': 'cannabis-cpr-reference-discovery-v1', 'rows': [row],
        'model_eligible': False, 'lookups': [item], 'claim_boundary': BOUNDARY,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    output = Path('data/reports/phase1-cpr-search.json')
    screen(source, raw, output, evidence_class='plant-CPR-homology-experimental-lead',
        additional_blockers=('P450-partnership-unverified', 'carrier-interface-unresolved',
            'not-eligible-for-exact-reaction-model'), claim_boundary=BOUNDARY)
    report = json.loads(output.read_text())
    report['model_eligible'] = False
    for result in report['rows']:
        result.update(model_eligible=False, carrier_interface_review=review,
            compatible_fnsii_partners=[], partner_status='unverified')
    output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    export_table(output, Path('data/derived/phase1-cpr-search.ndjson'))
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
