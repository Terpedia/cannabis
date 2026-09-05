"""Chemical witnesses for targets not rescued by the single CHI counterfactual."""
import hashlib
import json
from pathlib import Path
from .phase1_chalcone_reference_search import RID
from .phase1_net_flux import NetModel


def run():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in
        ('chalcone-addition-sensitivity', 'remaining-candidate-net', 'full-balanced-network')]
    sensitivity, candidate, network = [json.loads(p.read_text()) for p in paths]
    scenario = next(s for s in sensitivity['scenarios'] if s['id'] == 'one-forward-chalcone-addition')
    targets = [r for r in scenario['rows'] if r['status'] != 'exact-net-conversion-hypothesis']
    model = NetModel(network['reactions'], sensitivity['external_exchange_compound_ids'], scenario['forbidden_step_ids'])
    cheap = set(candidate['candidate_reaction_evidence_ids']) | {RID}
    costs = [1 if s['reaction_id'] in cheap else 1000 for s in model.steps]
    rows, gaps = [], {}
    for target in targets:
        result = model.solve(target['compound_id'], step_costs=costs)
        row = {k: target[k] for k in ('cannabisdb_id', 'compound_id', 'label')}
        rows.append({**row, **result})
        for step in result.get('steps', []):
            if step['reaction_id'] not in cheap:
                gap = gaps.setdefault(step['reaction_id'], {'reaction_id': step['reaction_id'], 'selected_uses': []})
                gap['selected_uses'].append({**step, 'target_id': target['cannabisdb_id']})
        print(target['cannabisdb_id'], result['status'], flush=True)
    used = {s['reaction_id'] for r in rows for s in r.get('steps', [])}
    reactions = [r for r in network['reactions'] if r['id'] in used]
    by_id = {r['id']: r for r in reactions}
    for rid, gap in gaps.items():
        gap['reaction'] = by_id[rid]
    compound_ids = set(sensitivity['external_exchange_compound_ids']) | {r['compound_id'] for r in rows} | {
        p['compound_id'] for r in reactions for side in ('left', 'right') for p in r[side]}
    report = {'schema': 'cannabis-chalcone-remaining-gap-witnesses-v1', 'model_eligible': False,
        'rows': rows, 'candidate_gaps': [gaps[rid] for rid in sorted(gaps)], 'reactions': reactions,
        'compounds': [c for c in network['compounds'] if c['id'] in compound_ids],
        'external_exchange_compound_ids': sensitivity['external_exchange_compound_ids'],
        'forbidden_step_ids': scenario['forbidden_step_ids'],
        'hypothetical_discounted_reaction_id': RID,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'target_records': len(rows), 'missing_equations_in_selected_witnesses': len(gaps)},
        'claim_boundary': 'Full balanced chemical catalog, not an enzyme-supported model. Positive extent cost is 1 for baseline candidate reactions and hypothetical forward CHI, 1000 for others. This does not minimize the number of missing enzymes or prove necessity. The CHI addition remains hypothetical; full net witnesses may require pre-existing pools. No new organic carbon inputs, evidence assignments, confirmed rescues or atom tracing.'}
    Path('data/reports/phase1-chalcone-remaining-gaps.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
