"""Exact carrier-matched odd-chain precursor supply audit; adds no chemistry."""
import hashlib
import json
import re
from pathlib import Path
from rdkit import Chem
from .phase1_catalog import stable_id
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_reaction_completion_net import validate_certificate

STARTER = 'structure:f006b29e743aebc0dc6d337f839bdacc0816851d971dd990e2ddef1091ddd047'


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-odd-chain-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    reactions, compounds, model = assemble_current(docs[1], docs[2], current, docs[3:])
    steps = {s['id']: s for s in model.steps}
    sm = compounds[STARTER]['smiles']
    prefix = re.match(r'^C+', sm).group()
    assert len(prefix) == 17
    tail = sm[len(prefix):]
    targets = []
    used_reactions = set()
    for length in range(3, 24, 2):
        exact = Chem.MolToSmiles(Chem.MolFromSmiles('C' * length + tail), isomericSmiles=True)
        cid = stable_id('structure', exact)
        if cid in compounds and compounds[cid]['smiles'] != exact:
            raise ValueError('Exact carrier identity mismatch')
        producing = []
        for step in model.steps:
            output = sum(p['coefficient'] for p in step['outputs'] if p['compound_id'] == cid)
            consumed = sum(p['coefficient'] for p in step['required_inputs'] if p['compound_id'] == cid)
            if output > consumed:
                producing.append(step)
        result = model.solve(cid)
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate({'compound_id': cid, **result}, steps, compounds,
                                 set(current['external_exchange_compound_ids']), current['co2_compound_id'])
            used_reactions.update(s['reaction_id'] for s in result['steps'])
        used_reactions.update(s['reaction_id'] for s in producing)
        row = {'acyl_carbon_count': length, 'compound_id': cid, 'smiles': exact,
               'present_in_model': cid in compounds, 'producing_steps': producing, 'net_result': result}
        targets.append(row)
        print(json.dumps({'acyl_carbons': length, 'producer_steps': len(producing), 'status': result['status']}), flush=True)
    used = {p['compound_id'] for rid in used_reactions for side in ('left', 'right') for p in reactions[rid][side]}
    report = {'schema': 'cannabis-carbon.phase1-odd-chain-supply.v1',
        'targets': targets, 'reactions': [reactions[r] for r in sorted(used_reactions)],
        'compounds': [compounds[c] for c in sorted(used)],
        'external_exchange_compound_ids': current['external_exchange_compound_ids'],
        'forbidden_step_ids': current['forbidden_step_ids'], 'co2_compound_id': current['co2_compound_id'],
        'reference_carrier_compound_id': STARTER,
        'summary': {'assessed_exact_odd_chain_acyl_CoAs': len(targets),
            'exact_net_certificates': sum(t['net_result']['status'] == 'exact-net-conversion-hypothesis' for t in targets),
            'added_equations': 0, 'inventory_coverage_gain_claimed': 0},
        'claim_boundary': 'Diagnostic precursor panel, not a replacement for the 6220-record inventory. '
            'Each chain uses the exact same charged, stereospecified CoA carrier as the C17 gap. '
            'Other CoA protonation/stereo forms and ACP are not merged. All producers retain complete required inputs and direction identifiers. '
            'Net certificates are permissive chemistry hypotheses with regenerated pre-existing pools, not confirmed Cannabis synthesis, uptake, startup, light or minimum-medium evidence. '
            'No new equation or exchange species is added. A missing or infeasible precursor is a pinned-model gap, not biological absence.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-odd-chain-supply.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
