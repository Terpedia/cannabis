"""Recover CHI references without assuming a populated Rhea cross-reference."""
import csv
import hashlib
import io
import json
from pathlib import Path
from .phase1_ureidoglycolate_broad_search import lookup
from .phase1_new_protein_search import run as screen, export_table

RID = 'balanced-equation:4f54b627e4bdef8e0ed5d73c4d6fafaef5222100c538665b7ffabbc9bb3c9568'
QUERY = 'ec:5.5.1.6 AND reviewed:true AND fragment:false'
RAW = Path('data/raw/phase1-chalcone-reference-search')
BOUNDARY = ('Reviewed EC-class discovery, not an exact substrate or Rhea annotation join. '
    'Individual references may differ in substrate specificity or lack catalysis. Preserve all '
    'search results as leads pending reference assay, catalytic-site and Cannabis sequence review. '
    'No physiological direction, net-route rescue or atom-tracing claim; model integration disabled.')


def build(item):
    parent = Path('data/reports/phase1-current-gap-priority.json')
    gap = next(r for r in json.loads(parent.read_text())['rows'] if r['reaction_id'] == RID)
    snapshot = Path(item['snapshot'])
    if hashlib.sha256(snapshot.read_bytes()).hexdigest() != item['sha256']:
        raise ValueError('Changed reference snapshot')
    records = list(csv.DictReader(io.StringIO(snapshot.read_text()), delimiter='\t'))
    if len(records) != item['records'] or len({r['Entry'] for r in records}) != len(records):
        raise ValueError('Missing or duplicated reference records')
    refs = [{'accession': r['Entry'], 'model_eligible': False,
        'exact_reaction_annotation_match': False, 'match_type': 'reviewed-EC-class-only',
        'source_records': [{'lookup_url': item['url'], 'query': item['query'], 'record': r}]} for r in records]
    selected = [{**u, 'target_ids': [t for t in u['target_ids'] if t in gap['remaining_target_ids']]}
        for u in gap['selected_uses']]
    return {'schema': 'cannabis-chalcone-reference-discovery-v1', 'model_eligible': False,
        'lookups': [item], 'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (parent, snapshot)},
        'rows': [{'reaction_id': RID, 'left': gap['reaction']['left'], 'right': gap['reaction']['right'],
            'sources': gap['reaction']['sources'], 'participants': gap['participants'],
            'target_ids': gap['remaining_target_ids'], 'priority_target_ids': gap['remaining_target_ids'],
            'hypothesis_ids': [], 'selected_uses': [u for u in selected if u['target_ids']],
            'historical_selected_uses': gap['selected_uses'], 'prior_searches': gap['prior_searches'],
            'reference_matches': sorted(refs, key=lambda r: r['accession']), 'model_eligible': False}],
        'summary': {'reference_leads': len(refs), 'references_without_rhea_field': sum(not r.get('Rhea ID', '').strip() for r in records),
            'current_target_records': len(gap['remaining_target_ids']), 'new_exact_enzyme_assignments': 0},
        'claim_boundary': BOUNDARY}


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    report = build(lookup('reviewed-ec5516', QUERY, RAW))
    source = Path('data/reports/phase1-chalcone-references.json')
    source.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)
    output = Path('data/reports/phase1-chalcone-search.json')
    screen(source, RAW, output, evidence_class='reviewed-EC-class-substrate-unverified-homology-lead',
        additional_blockers=('EC-class-not-exact-substrate-assignment', 'not-eligible-for-exact-reaction-model'),
        claim_boundary=BOUNDARY)
    result = json.loads(output.read_text())
    result['model_eligible'] = False
    for row in result['rows']:
        row['model_eligible'] = False
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    export_table(output, Path('data/derived/phase1-chalcone-search.ndjson'))


if __name__ == '__main__':
    run()
