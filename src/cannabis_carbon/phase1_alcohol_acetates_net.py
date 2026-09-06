"""Full-inventory sensitivity to six exact alcohol acetylation hypotheses."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build


def run():
    root = Path('data/reports')
    parent_path = root / 'phase1-c17-elongation-net.json'
    parent = json.loads(parent_path.read_bytes())
    paths = [root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             parent_path, root / 'phase1-alcohol-acetates.json'] + [root / n for n in parent['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale model input')
    previous = docs[4:]
    certs = {c['compound_id']: c for doc in previous for c in doc.get('certificates', []) + doc.get('new_certificates', [])}
    if len(certs) + len(parent['new_certificates']) != 2733:
        raise ValueError('Missing baseline certificate')
    original = {'certificates': list(certs.values()), 'external_exchange_compound_ids': parent['external_exchange_compound_ids']}
    layers = [{**doc, 'forbidden_step_ids': doc.get('forbidden_step_ids', [])} for doc in [*previous, parent] if 'added_reactions' in doc]
    report = build(docs[0], docs[1], original, parent, docs[3], prior_layers=layers,
                   extra_forward_types=('alcohol-acetylation',))
    report.update({'schema': 'cannabis-carbon.phase1-alcohol-acetates-net.v1',
        'baseline_certificate_reports': parent['baseline_certificate_reports'] + [parent_path.name],
        'lipid_evidence_report': 'phase1-alcohol-acetates.json',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'claim_boundary': 'Alcohol acetylation substrate-scope sensitivity, not target-specific enzyme or confirmed Cannabis activity. '
            'New equations are proposed forward-only; inherited exact identities, directions and permissive inorganic boundary unchanged. '
            'Positive net certificates allow regenerated pre-existing pools; startup, nutrient uptake, light, energetics, growth and compartments are not established.'})
    (root / 'phase1-alcohol-acetates-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
