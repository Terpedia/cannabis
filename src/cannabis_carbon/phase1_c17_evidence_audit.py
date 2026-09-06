"""Evidence and external-input annotations for saved C17-gain certificates."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_reaction_completion_net import validate_certificate


def annotate(reaction, step, compounds):
    inputs = {p['compound_id'] for p in step['required_inputs']}
    outputs = {p['compound_id'] for p in step['outputs']}
    input_smiles = {compounds[c]['smiles'] for c in inputs}
    output_smiles = {compounds[c]['smiles'] for c in outputs}
    flags = []
    if 'OO' in input_smiles and 'O=O' in output_smiles:
        flags.append('peroxide-consuming-oxygen-producing-direction-review')
    if 'O=C=O' in input_smiles:
        flags.append('co2-consuming-step-not-automatically-photosynthetic-fixation')
    if reaction.get('hypothesis_type'):
        evidence = 'proposed-gap-filling-reaction'
    elif reaction.get('sources'):
        evidence = 'catalog-linked-equation-not-Cannabis-activity-proof'
    else:
        evidence = 'source-link-or-stoichiometric-completion-review-required'
    return {'reaction_id': reaction['id'], 'evidence_class': evidence,
            'hypothesis_type': reaction.get('hypothesis_type'),
            'direction_mode': step['direction_mode'], 'review_flags': flags,
            'cannabis_physiological_direction_established_by_this_audit': False,
            'enzyme_assignment_established_by_this_audit': False,
            'source_record': reaction}


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-c17-elongation-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             root / 'phase1-odd-chain-supply.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    reactions, compounds, model = assemble_current(docs[1], docs[2], current, docs[4:])
    steps = {s['id']: s for s in model.steps}
    precursor = next(t for t in docs[3]['targets'] if t['acyl_carbon_count'] == 15)
    certificates = [('independent-C15-precursor-certificate', None, {'compound_id': precursor['compound_id'], **precursor['net_result']})]
    for c in current['new_certificates']:
        t = next(t for t in current['targets'] if t['compound_id'] == c['compound_id'])
        certificates.append(('new-inventory-target-certificate', t['cannabisdb_id'], c))
    rows = []
    for role, record_id, cert in certificates:
        validate_certificate(cert, steps, compounds, set(current['external_exchange_compound_ids']), current['co2_compound_id'])
        annotated = []
        for s in cert['steps']:
            step = steps[s['step_id']]
            annotated.append({**s, **annotate(reactions[s['reaction_id']], step, compounds),
                              'required_inputs': step['required_inputs'], 'outputs': step['outputs']})
        inputs = [{'compound_id': c, 'smiles': compounds[c]['smiles'],
                   'amount_per_target': str(Fraction(n) / Fraction(cert['target_amount'])),
                   'carbon_count': compounds[c]['carbon_count'],
                   'nutrient_uptake_established': False}
                  for c, n in cert['external_net_consumption'].items()]
        rows.append({'role': role, 'cannabisdb_id': record_id, 'compound_id': cert['compound_id'],
                     'certificate': cert, 'steps': annotated, 'net_inputs': inputs,
                     'summary': {'step_count': len(annotated),
                         'evidence_class_counts': dict(Counter(s['evidence_class'] for s in annotated)),
                         'review_flag_step_counts': dict(Counter(f for s in annotated for f in s['review_flags']))}})
        print(json.dumps({'role': role, 'cannabisdb_id': record_id, **rows[-1]['summary']}), flush=True)
    used = {p['compound_id'] for row in rows for s in row['steps'] for side in ('required_inputs', 'outputs') for p in s[side]}
    report = {'schema': 'cannabis-carbon.phase1-c17-evidence-audit.v1', 'certificates': rows,
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'audited_certificates': len(rows), 'new_inventory_target_certificates': len(current['new_certificates']),
                    'new_enzyme_assignments': 0, 'new_confirmed_Cannabis_pathways': 0, 'coverage_change': 0},
        'claim_boundary': 'Annotations describe the exact saved certificates, not all feasible alternatives or essential inputs. '
            'The independent C15 precursor certificate is not substituted for the upstream steps of any target certificate. '
            'A catalog source link does not establish target-organism activity, physiological direction or energy feasibility. '
            'Peroxide/oxygen flags are exact participant-pattern review triggers, not thermodynamic impossibility proofs. '
            'Net input amounts are certificate bookkeeping, not medium requirements or measured yields. '
            'All original input species and pre-existing-pool assumptions remain; no reaction, bound, coverage or enzyme assignment changes.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-c17-evidence-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
