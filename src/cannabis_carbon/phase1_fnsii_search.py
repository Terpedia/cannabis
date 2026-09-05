"""Full-proteome FNS-II leads, explicitly outside the exact reaction model."""
import hashlib
import json
from pathlib import Path
from .phase1_new_protein_search import run as screen, export_table

BOUNDARY = 'Generic FNS-II homology leads only. Exact naringenin/apigenin specificity, protein-bound reductase partnership and Cannabis activity remain unverified. No substitution for the exact FNS-I equation, free-flavin assumption, model promotion or atom-tracing claim.'


def run():
    parent_path = Path('data/reports/phase1-fnsii-alternative-audit.json')
    parent = json.loads(parent_path.read_text())
    for path, digest in parent['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Changed alternative audit source')
    source_record = next(r for r in parent['source_records'] if r['record']['rule_id'] == 'RHEA:57681')
    targets = sorted({u['target_id'] for u in parent['parent_fnsi_gap']['selected_uses']})
    row = {'reaction_id': 'RHEA:57681', 'model_eligible': False,
        'target_ids': targets, 'priority_target_ids': targets, 'hypothesis_ids': [],
        'generic_source_record': source_record, 'carrier_review': parent['review'],
        'related_fnsi_gap_id': parent['parent_fnsi_gap']['reaction_id'],
        'reference_matches': [{**r, 'exact_reaction_annotation_match': False,
            'match_type': 'generic-FNS-II-EC-class; not-exact-FNS-I-equation'} for r in parent['reference_leads']]}
    source = Path('data/reports/phase1-fnsii-references.json')
    discovery = {'schema': 'cannabis-fnsii-reference-discovery-v1', 'model_eligible': False,
        'rows': [row], 'claim_boundary': BOUNDARY,
        'source_sha256': {str(parent_path): hashlib.sha256(parent_path.read_bytes()).hexdigest()}}
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    output = Path('data/reports/phase1-fnsii-search.json')
    screen(source, Path('data/raw/fnsii-search'), output,
        evidence_class='generic-FNS-II-homology-experimental-lead',
        additional_blockers=('generic-flavonoid-substituents-unresolved', 'protein-bound-carrier-partner-unresolved',
            'not-eligible-for-exact-reaction-model'), claim_boundary=BOUNDARY)
    result = json.loads(output.read_text())
    result['model_eligible'] = False
    for r in result['rows']:
        r.update(model_eligible=False, generic_source_record=source_record, carrier_review=parent['review'],
            related_fnsi_gap_id=row['related_fnsi_gap_id'])
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    export_table(output, Path('data/derived/phase1-fnsii-search.ndjson'))


if __name__ == '__main__':
    run()
