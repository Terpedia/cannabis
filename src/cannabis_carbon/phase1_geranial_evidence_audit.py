"""Full participant/direction/input evidence for all gains since the C17 Pages release."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from .phase1_c17_evidence_audit import annotate
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_reaction_completion_net import validate_certificate


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             root / 'phase1-geranial-reduction.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    reactions, compounds, model = assemble_current(docs[1], docs[2], current, docs[4:])
    steps = {s['id']: s for s in model.steps}
    parent = next(doc for doc in docs[4:] if doc['schema'] == 'cannabis-carbon.phase1-alcohol-acetates-net.v1')
    certs = parent['new_certificates'] + current['new_certificates']
    rows = []
    for cert in certs:
        validate_certificate(cert, steps, compounds, set(current['external_exchange_compound_ids']), current['co2_compound_id'])
        target = next(t for t in current['targets'] if t['compound_id'] == cert['compound_id'])
        annotated = []
        for s in cert['steps']:
            step = steps[s['step_id']]; reaction = reactions[s['reaction_id']]
            annotation = annotate(reaction, step, compounds)
            annotation['source_evidence_type'] = reaction.get('source_evidence_type')
            if reaction.get('hypothesis_type') == 'geranial-S-citronellal':
                annotation['evidence_class'] = 'substrate-specific-cascade-biochemistry-outside-Cannabis'
                annotation['biochemical_evidence'] = docs[3]['biochemical_evidence']
                annotation['review_flags'].append('external-organism-evidence-not-Cannabis-enzyme-assignment')
            annotated.append({**s, **annotation, 'required_inputs': step['required_inputs'], 'outputs': step['outputs']})
        inputs = [{'compound_id': c, 'smiles': compounds[c]['smiles'],
                   'amount_per_target': str(Fraction(n) / Fraction(cert['target_amount'])),
                   'carbon_count': compounds[c]['carbon_count'], 'nutrient_uptake_established': False}
                  for c, n in cert['external_net_consumption'].items()]
        rows.append({'cannabisdb_id': target['cannabisdb_id'], 'label': target['label'],
            'compound_id': cert['compound_id'], 'certificate': cert, 'steps': annotated, 'net_inputs': inputs,
            'summary': {'step_count': len(annotated),
                'evidence_class_counts': dict(Counter(s['evidence_class'] for s in annotated)),
                'review_flag_step_counts': dict(Counter(f for s in annotated for f in s['review_flags']))}})
        print(json.dumps({'cannabisdb_id': target['cannabisdb_id'], **rows[-1]['summary']}), flush=True)
    used = {p['compound_id'] for row in rows for s in row['steps'] for side in ('required_inputs', 'outputs') for p in s[side]}
    report = {'schema': 'cannabis-carbon.phase1-geranial-evidence-audit.v1', 'certificates': rows,
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'audited_certificates': len(rows), 'new_enzyme_assignments': 0,
                    'new_confirmed_Cannabis_pathways': 0, 'coverage_change': 0},
        'claim_boundary': 'Audit of four saved gains since the C17 release, not all alternatives or essential inputs. '
            'External-organism substrate-specific assay evidence is distinguished from reaction-class analogy and catalog links. '
            'No Cannabis enzyme, physiological direction, energetic feasibility, uptake or growth requirement is established. '
            'CO2 and peroxide/oxygen flags are participant-pattern review triggers, not thermodynamic proofs. '
            'Input amounts are exact net-certificate bookkeeping, not a medium recipe. All original bounds and pool assumptions remain.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-geranial-evidence-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
