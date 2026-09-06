"""Exact supply and obstruction audit for unresolved alcohol-acetate routes."""
import hashlib
import json
from pathlib import Path
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_net_obstruction import solve as obstruction
from .phase1_reaction_completion_net import validate_certificate


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-alcohol-acetates-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             root / 'phase1-alcohol-acetates.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    reactions, compounds, model = assemble_current(docs[1], docs[2], current, docs[4:])
    steps = {s['id']: s for s in model.steps}
    rows = []; used = set()
    for proposal in docs[3]['targets']:
        if 'hypothesis_id' not in proposal:
            continue
        target = next(t for t in current['targets'] if t['cannabisdb_id'] == proposal['cannabisdb_id'])
        if target['net_status'] != 'solver-reported-infeasible':
            continue
        cid = proposal['alcohol_compound_id']
        producers = [s for s in model.steps if
            sum(p['coefficient'] for p in s['outputs'] if p['compound_id'] == cid) >
            sum(p['coefficient'] for p in s['required_inputs'] if p['compound_id'] == cid)]
        result = model.solve(cid)
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate({'compound_id': cid, **result}, steps, compounds,
                                 set(current['external_exchange_compound_ids']), current['co2_compound_id'])
            used.update(s['reaction_id'] for s in result['steps'])
        proof = obstruction(model, target['compound_id'])
        used.update(s['reaction_id'] for s in producers)
        used.update(steps[s['step_id']]['reaction_id'] for s in proof.get('weighted_steps', []))
        used.add(proposal['hypothesis_id'])
        rows.append({'cannabisdb_id': target['cannabisdb_id'], 'label': target['label'],
            'compound_id': target['compound_id'], 'acetylation_hypothesis_id': proposal['hypothesis_id'],
            'alcohol_compound_id': cid, 'alcohol_smiles': compounds[cid]['smiles'],
            'alcohol_producing_steps': producers, 'alcohol_net_result': result, 'target_obstruction': proof})
        print(json.dumps({'target': target['cannabisdb_id'], 'alcohol_producers': len(producers),
            'alcohol_net_status': result['status'], 'target_obstruction': proof['status'],
            'weighted_compounds': len(proof.get('weights', {}))}), flush=True)
    cids = {p['compound_id'] for rid in used for side in ('left', 'right') for p in reactions[rid][side]}
    cids.update(c for row in rows for c in row['target_obstruction'].get('weights', {}))
    report = {'schema': 'cannabis-carbon.phase1-acetate-precursor-gaps.v1', 'targets': rows,
        'reactions': [reactions[r] for r in sorted(used)], 'compounds': [compounds[c] for c in sorted(cids)],
        'external_exchange_compound_ids': current['external_exchange_compound_ids'],
        'forbidden_step_ids': current['forbidden_step_ids'],
        'summary': {'assessed_blocked_records': len(rows),
            'exact_target_obstructions': sum(r['target_obstruction']['status'] == 'exact-stoichiometric-obstruction' for r in rows),
            'alcohols_without_producing_steps': sum(not r['alcohol_producing_steps'] for r in rows),
            'balanced_equations': len(reactions), 'allowed_steps': len(model.steps), 'coverage_change': 0},
        'claim_boundary': 'Exact precursor and target-obstruction audit in the pinned permissive chemistry model. '
            'A nonnegative obstruction weight set cannot increase under any allowed directed reaction; positive target production '
            'without internal-pool depletion is therefore impossible in that model when a proof validates. '
            'A missing proof or solver failure establishes no impossibility. No biological absence, unique missing enzyme, '
            'minimum cut, startup synthesis, minimum medium or physiological direction is established. '
            'The four-record diagnostic does not replace the complete 6220-record objective or change any model bounds.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-acetate-precursor-gaps.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
