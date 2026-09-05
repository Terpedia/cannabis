"""EC-class FNS-I leads for the exact remaining oxygenase equation."""
import csv
import hashlib
import json
from pathlib import Path
from .phase1_ureidoglycolate_broad_search import lookup
from .phase1_new_protein_search import run as screen, export_table

RAW = Path('data/raw/phase1-flavone-reference-search')
QUERY = 'ec:1.14.20.5 AND reviewed:true AND fragment:false'
BOUNDARY = 'Reviewed FNS-I EC-class leads only. Exact substrate, charge, specificity and Cannabis activity require review. FNS-II cofactor chemistry is not interchangeable. No model promotion, biological route rescue, or atom tracing.'


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    item = lookup('reviewed-fnsi', QUERY, RAW)
    source_path = Path('data/reports/phase1-chalcone-remaining-gaps.json')
    parent = json.loads(source_path.read_text())
    if len(parent['candidate_gaps']) != 1:
        raise ValueError('Expected one reviewed target gap')
    gap = parent['candidate_gaps'][0]
    with Path(item['snapshot']).open() as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    tids = sorted({u['target_id'] for u in gap['selected_uses']})
    row = {'reaction_id': gap['reaction_id'], 'left': gap['reaction']['left'], 'right': gap['reaction']['right'],
        'sources': gap['reaction']['sources'], 'target_ids': tids, 'priority_target_ids': tids,
        'hypothesis_ids': [], 'selected_uses': [{**u, 'target_ids': [u['target_id']]} for u in gap['selected_uses']],
        'model_eligible': False, 'reference_matches': [{'accession': r['Entry'], 'model_eligible': False,
            'exact_reaction_annotation_match': False, 'match_type': 'reviewed-EC-class-only',
            'source_record': r, 'source_url': item['url']} for r in records]}
    source = Path('data/reports/phase1-flavone-references.json')
    report = {'schema': 'cannabis-flavone-reference-discovery-v1', 'model_eligible': False,
        'rows': [row], 'lookups': [item], 'claim_boundary': BOUNDARY,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source_path, Path(item['snapshot']))}}
    source.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    output = Path('data/reports/phase1-flavone-search.json')
    screen(source, RAW, output, evidence_class='reviewed-EC-class-specificity-unverified-homology-lead',
        additional_blockers=('not-eligible-for-exact-reaction-model', 'reference-exact-substrate-review-pending'), claim_boundary=BOUNDARY)
    result = json.loads(output.read_text())
    result['model_eligible'] = False
    for r in result['rows']:
        r['model_eligible'] = False
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    export_table(output, Path('data/derived/phase1-flavone-search.ndjson'))


if __name__ == '__main__':
    run()
