"""Current-model uptake restriction for four recent gains; not full-inventory coverage."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_medium_boundary import UptakeLimitedModel, blocked_species
from .phase1_reaction_completion_net import validate_certificate
from .phase1_net_obstruction import solve as obstruction


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
             root / 'phase1-geranial-evidence-audit.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    reactions, compounds, _ = assemble_current(docs[1], docs[2], current, docs[4:])
    exchange = set(current['external_exchange_compound_ids'])
    blocked = blocked_species(compounds, exchange)
    uptake = exchange - {c['compound_id'] for c in blocked}
    model = UptakeLimitedModel(list(reactions.values()), exchange, uptake, current['forbidden_step_ids'])
    steps = {s['id']: s for s in model.steps}
    rows = []
    for audited in docs[3]['certificates']:
        cid = audited['compound_id']
        result = {'compound_id': cid, **model.solve(cid)}
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate(result, steps, compounds, uptake, current['co2_compound_id'])
        row = {k: audited[k] for k in ('cannabisdb_id', 'label', 'compound_id')}
        row.update({'restricted_net_result': result,
            'prior_certificate_blocked_inputs': sorted(set(audited['certificate']['external_net_consumption']) - uptake)})
        if result['status'] != 'exact-net-conversion-hypothesis':
            row['obstruction_result'] = obstruction(model.model, cid)
        rows.append(row)
        print(json.dumps({'cannabisdb_id': row['cannabisdb_id'], 'status': result['status'],
                          'obstruction': row.get('obstruction_result', {}).get('status')}), flush=True)
    used = set(exchange)
    used_reactions = set()
    for row in rows:
        for s in row['restricted_net_result'].get('steps', []):
            used_reactions.add(s['reaction_id'])
            used.update(p['compound_id'] for side in ('required_inputs', 'outputs') for p in steps[s['step_id']][side])
        used.update(row.get('obstruction_result', {}).get('weights', {}))
        used.add(row['compound_id'])
    report = {'schema': 'cannabis-carbon.phase1-recent-gains-medium.v1', 'targets': rows,
        'compounds': [compounds[c] for c in sorted(used)],
        'certificate_reactions': [reactions[r] for r in sorted(used_reactions)],
        'blocked_inputs': blocked, 'allowed_uptake_compound_ids': sorted(uptake),
        'allowed_external_output_compound_ids': sorted(exchange),
        'forbidden_step_ids': current['forbidden_step_ids'], 'co2_compound_id': current['co2_compound_id'],
        'summary': {'assessed_target_records': len(rows), 'full_model_inventory_records': len(current['targets']),
            'balanced_equations': len(reactions), 'allowed_steps': len(steps),
            'net_status_counts': dict(Counter(r['restricted_net_result']['status'] for r in rows)),
            'exact_obstructions': sum(r.get('obstruction_result', {}).get('status') == 'exact-stoichiometric-obstruction' for r in rows)},
        'claim_boundary': 'Four-target diagnostic in the geranial-reduction scenario, not a full-inventory restriction run. '
            'No net import of iron-sulfur clusters or hydrogen peroxide; their synthesis and disposal remain allowed. '
            'Other permissive inputs, directions and recycled initial pools remain. A retained certificate demonstrates an '
            'alternative net conversion under these restrictions, not physiological nutrient uptake, growth or a minimum medium. '
            'An exact obstruction applies only to this pinned model; solver failure alone is not an impossibility proof.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-recent-gains-medium.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
