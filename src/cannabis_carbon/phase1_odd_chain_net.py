"""Full-inventory sensitivity to proposed odd-chain elongation cycles."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build


def run(parent_name='alkane-net', hypothesis_name='odd-chain-elongation', output_name='odd-chain-net'):
    root = Path('data/reports')
    read = lambda p: json.loads(p.read_bytes())
    parent_path = root / ('phase1-' + parent_name + '.json')
    parent = read(parent_path)
    paths = [root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             parent_path, root / ('phase1-' + hypothesis_name + '.json')] + [root / n for n in parent['baseline_certificate_reports']]
    docs = [read(p) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale model input')
    previous = docs[4:]
    certs = {c['compound_id']: c for doc in previous for c in doc.get('certificates', []) + doc.get('new_certificates', [])}
    if len(certs) + len(parent['new_certificates']) != 2729:
        raise ValueError('Missing baseline certificates')
    original = {'certificates': list(certs.values()), 'external_exchange_compound_ids': parent['external_exchange_compound_ids']}
    layers = [{**doc, 'forbidden_step_ids': doc.get('forbidden_step_ids', [])}
              for doc in [*previous, parent] if 'added_reactions' in doc]
    report = build(docs[0], docs[1], original, parent, docs[3], prior_layers=layers,
                   extra_forward_types=('odd-chain-elongation',))
    report.update({'schema': 'cannabis-carbon.phase1-' + output_name + '.v1',
        'baseline_certificate_reports': parent['baseline_certificate_reports'] + [parent_path.name],
        'lipid_evidence_report': 'phase1-' + hypothesis_name + '.json',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'claim_boundary': 'Odd-chain elongation substrate-scope sensitivity, not established Cannabis activity. '
            'New equations are proposed forward-only. Exact inherited identities, directions and permissive inorganic boundary retained. '
            'Positive exact net certificates allow regenerated pre-existing pools; nutrient uptake, growth, light, energetics and startup are not validated.'})
    (root / ('phase1-' + output_name + '.json')).write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
