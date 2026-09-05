"""Reprioritize saved chemical witnesses against the current candidate model."""
import hashlib
import json
from pathlib import Path


def build():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in
             ('remaining-weighted-routes', 'remaining-candidate-net', 'full-balanced-network')]
    diagnostic, model, network = [json.loads(p.read_text()) for p in paths]
    scenario = next(s for s in model['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    targets = {t['cannabisdb_id']: t for t in scenario['targets']}
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    rows = []
    for gap in diagnostic['candidate_gaps']:
        rid = gap['reaction_id']
        remaining = [t for t in gap['target_ids'] if targets[t]['net_status'] not in
                     ('exact-net-conversion-hypothesis', 'explicit-exchange-species; not a synthesis target')]
        reaction = reactions[rid]
        rows.append({'reaction_id': rid, 'already_in_candidate_model': rid in model['candidate_reaction_evidence_ids'],
            'historical_target_ids': gap['target_ids'], 'remaining_target_ids': remaining,
            'remaining_target_statuses': {t: targets[t]['net_status'] for t in remaining},
            'remaining_target_count': len(remaining), 'selected_uses': gap['selected_uses'],
            'prior_searches': gap['prior_searches'], 'reaction': reaction,
            'participants': [compounds[cid] for cid in sorted({p['compound_id'] for side in ('left', 'right') for p in reaction[side]})]})
    rows.sort(key=lambda r: (r['already_in_candidate_model'], -r['remaining_target_count'], r['reaction_id']))
    return {'schema': 'cannabis-current-gap-priority-v1', 'rows': rows,
        'scenario_id': scenario['id'],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'saved_chemical_witness_gaps': len(rows),
                    'already_in_candidate_model': sum(r['already_in_candidate_model'] for r in rows)},
        'claim_boundary': 'Ranking counts unresolved targets using each reaction in saved chemical witnesses, not predicted rescues, essentiality or physiological flux. These witnesses are not rerun or proven optimal for the current model. Historical searches are preserved as historical, not a complete latest-search inventory. Exact participants and hypothetical directions remain explicit.'}


if __name__ == '__main__':
    report = build()
    Path('data/reports/phase1-current-gap-priority.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))
