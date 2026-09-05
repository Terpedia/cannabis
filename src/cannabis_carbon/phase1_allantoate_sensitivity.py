"""Re-solve affected chemical targets without hypothetical allantoate condensation."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_net_flux import NetModel

RID = 'balanced-equation:641e06d2caccbcc147af3b55c3911a1e75bc4eb95cf7ebaa7999e4b4f2cca4f9'


def run(*, additional_forbidden=(), name='allantoate-sensitivity', review_path=None):
    paths = [Path('data/reports/phase1-' + n + '.json') for n in
             ('current-gap-priority', 'remaining-candidate-net', 'full-balanced-network')]
    priority, candidate, network = [json.loads(p.read_text()) for p in paths]
    reference_path = Path('data/raw/allantoate-direction-review/Q9W6S5.json')
    reference = json.loads(reference_path.read_text())
    gap = next(r for r in priority['rows'] if r['reaction_id'] == RID)
    scenario = next(s for s in candidate['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    forbidden = scenario['forbidden_step_ids'] + [RID + ':hypothetical-left-to-right'] + list(additional_forbidden)
    model = NetModel(network['reactions'], candidate['external_exchange_compound_ids'], forbidden)
    costs = [1 if s['reaction_id'] in candidate['candidate_reaction_evidence_ids'] else 1000 for s in model.steps]
    targets = {t['cannabisdb_id']: t for t in network['targets']}
    results = {}
    rows = []
    for tid in gap['remaining_target_ids']:
        target = targets[tid]
        cid = target['compound_id']
        if cid not in results:
            results[cid] = model.solve(cid, step_costs=costs)
        rows.append({'cannabisdb_id': tid, 'compound_id': cid, 'label': target['label'], **results[cid]})
        print(tid, results[cid]['status'], flush=True)
    used = {s['reaction_id'] for r in rows for s in r.get('steps', [])} | {RID}
    reactions = [r for r in network['reactions'] if r['id'] in used]
    ids = set(candidate['external_exchange_compound_ids']) | {r['compound_id'] for r in rows} | {
        p['compound_id'] for r in reactions for side in ('left', 'right') for p in r[side]}
    report = {'schema': 'cannabis-allantoate-condensation-sensitivity-v1', 'rows': rows,
        'forbidden_step_ids': forbidden, 'external_exchange_compound_ids': candidate['external_exchange_compound_ids'],
        'reactions': reactions, 'compounds': [c for c in network['compounds'] if c['id'] in ids],
        'reference_activity': [c for c in reference['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY'],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths + [reference_path]},
        'summary': dict(Counter(r['status'] for r in rows)),
        'claim_boundary': 'Affected-target chemical sensitivity only, not the candidate-enzyme model. Excludes hypothetical allantoate condensation without claiming irreversibility; source reaction notation alone is not physiological direction. All other eight exclusions and explicit exchange inputs retained. Weighted extents prioritize candidate reactions but do not establish physiology. Atom tracing deferred.'}
    if review_path is not None:
        review_path = Path(review_path)
        report['direction_review'] = json.loads(review_path.read_text())
        report['source_sha256'][str(review_path)] = hashlib.sha256(review_path.read_bytes()).hexdigest()
        report['claim_boundary'] += ' Additional direction exclusions are sensitivity assumptions documented in direction_review, not experimentally established irreversible bounds.'
        report['schema'] = 'cannabis-ureide-condensation-sensitivity-v1'
    Path('data/reports/phase1-' + name + '.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
