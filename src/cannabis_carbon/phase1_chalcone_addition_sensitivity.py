"""Counterfactual impact of one forward CHI equation, not an enzyme promotion."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_chalcone_reference_search import RID
from .phase1_net_flux import NetModel


def run():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in
        ('remaining-candidate-net', 'full-balanced-network', 'current-gap-priority')]
    candidate, network, priority = [json.loads(p.read_text()) for p in paths]
    review = Path('data/curation/chalcone-primary-assay-review.json')
    if json.loads(review.read_text())['reaction_id'] != RID:
        raise ValueError('Wrong reviewed reaction')
    gap = next(r for r in priority['rows'] if r['reaction_id'] == RID)
    reaction = next(r for r in network['reactions'] if r['id'] == RID)
    baseline_reactions = [r for r in network['reactions'] if r['id'] in candidate['candidate_reaction_evidence_ids']]
    if len(baseline_reactions) != len(candidate['candidate_reaction_evidence_ids']):
        raise ValueError('Incomplete baseline reaction reconstruction')
    if RID in {r['id'] for r in baseline_reactions}:
        raise ValueError('Not an addition to baseline')
    baseline = next(s for s in candidate['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    exchanges = candidate['external_exchange_compound_ids']
    compounds = {c['id']: c for c in network['compounds']}
    if [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] != ['O=C=O']:
        raise ValueError('Unexpected carbon exchange')
    targets = {t['cannabisdb_id']: t for t in network['targets']}
    scenarios = []
    for name, reactions, forbidden in (
        ('baseline', baseline_reactions, baseline['forbidden_step_ids']),
        ('one-forward-chalcone-addition', baseline_reactions + [reaction],
            baseline['forbidden_step_ids'] + [RID + ':hypothetical-right-to-left'])):
        model = NetModel(reactions, exchanges, forbidden)
        cache, rows = {}, []
        for tid in gap['remaining_target_ids']:
            target = targets[tid]
            cid = target['compound_id']
            if cid not in cache:
                cache[cid] = model.solve(cid)
            rows.append({'cannabisdb_id': tid, 'compound_id': cid, 'label': target['label'], **cache[cid]})
        scenarios.append({'id': name, 'forbidden_step_ids': forbidden, 'rows': rows,
            'summary': dict(Counter(r['status'] for r in rows)), 'reaction_count': len(reactions)})
        print(name, json.dumps(scenarios[-1]['summary']), flush=True)
    report = {'schema': 'cannabis-chalcone-addition-sensitivity-v1', 'model_eligible': False,
        'hypothetical_added_reaction_id': RID, 'scenarios': scenarios,
        'reactions': baseline_reactions + [reaction],
        'external_exchange_compound_ids': exchanges, 'compounds': network['compounds'],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths + [review]},
        'claim_boundary': 'Only the seven preselected affected records are assessed, not all Cannabis metabolites. Baseline is rerun unchanged; the counterfactual adds exactly one balanced forward equation, without new enzyme evidence or reverse activity. Existing eight reverse exclusions and all explicit exchange assumptions retained. Net certificates may require pre-existing conserved pools and do not prove startup, expression, compartments, thermodynamics, or in-vivo flux. No published model or confirmed completeness change.'}
    Path('data/reports/phase1-chalcone-addition-sensitivity.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
