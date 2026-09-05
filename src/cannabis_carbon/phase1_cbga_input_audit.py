"""Exact participant supply probes for the three candidate-linked CBGA steps."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows

SOURCE_IDS = ('RHEA:34120', 'RHEA:34124', 'RHEA:34128')


def build(network, candidate):
    selected = []
    for sid in SOURCE_IDS:
        matches = [r for r in network['reactions'] if any(s['source_reaction_id'] == sid for s in r['sources'])]
        if len(matches) != 1:
            raise ValueError('Ambiguous reaction source')
        r = matches[0]
        if r['id'] not in candidate['candidate_reaction_evidence_ids']:
            raise ValueError('Expected upstream candidate link absent')
        selected.append({**r, 'candidate_evidence_ids': candidate['candidate_reaction_evidence_ids'][r['id']]})
    ids = {m['compound_id'] for r in selected for side in ('left', 'right') for m in r[side]}
    compounds = {c['id']: c for c in network['compounds']}
    exchange = candidate['external_exchange_compound_ids']
    assert {c for c in exchange if compounds[c]['carbon_count']} == {candidate['co2_compound_id']}
    rows = []
    reactions = [r for r in network['reactions'] if r['id'] in candidate['candidate_reaction_evidence_ids']]
    for scenario in candidate['scenarios']:
        model = NetModel(reactions, exchange, scenario['forbidden_step_ids'])
        for cid in sorted(ids):
            rows.append({'id': scenario['id'] + ':' + cid, 'scenario_id': scenario['id'],
                'compound_id': cid, **model.solve(cid)})
    used = {s['reaction_id'] for r in rows for s in r.get('steps', [])} | {r['id'] for r in selected}
    needed = set(exchange) | ids | {m['compound_id'] for r in reactions if r['id'] in used for side in ('left', 'right') for m in r[side]}
    return {'schema': 'cannabis-carbon.phase1-cbga-input-audit.v1', 'selected_source_ids': list(SOURCE_IDS),
        'selected_reactions': selected, 'probe_compound_ids': sorted(ids), 'results': rows,
        'reactions': [r for r in reactions if r['id'] in used], 'compounds': [compounds[c] for c in sorted(needed)],
        'external_exchange_compound_ids': exchange, 'co2_compound_id': candidate['co2_compound_id'],
        'scenario_constraints': {s['id']: s['forbidden_step_ids'] for s in candidate['scenarios']},
        'summary': {'selected_candidate_equations': len(selected), 'exact_participants': len(ids),
            'status_counts_by_scenario': {s['id']: dict(Counter(r['status'] for r in rows if r['scenario_id'] == s['id'])) for s in candidate['scenarios']}},
        'claim_boundary': 'All inputs and products of three exact candidate-linked steps retained. Independent net-production probes are not joint supply, startup or physiological direction proof. A failure to accumulate a cofactor as net output does not show a regenerated pool is unavailable. Solver infeasibility is model-limited, not biological absence. No organic carbon seeds, identity merges, enzyme promotions or target-count changes. Atom tracing remains deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-synthase-candidate-net')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-cbga-input-audit.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('result', 'results'), ('reaction', 'reactions'), ('compound', 'compounds')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection in groups:
        rows.extend((kind, r['id'], r) for r in report[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-cbga-input-audit.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
