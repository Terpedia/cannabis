"""Prioritize exact catalog gaps for every unresolved historical target."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows


def build(network, candidate, searches):
    scenario = next(s for s in candidate['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    exchange = candidate['external_exchange_compound_ids']
    evidence = candidate['candidate_reaction_evidence_ids']
    compounds = {c['id']: c for c in network['compounds']}
    if {c for c in exchange if compounds[c]['carbon_count']} != {candidate['co2_compound_id']}:
        raise ValueError('Organic carbon exchange forbidden')
    model = NetModel(network['reactions'], exchange, scenario['forbidden_step_ids'])
    costs = [1 if s['reaction_id'] in evidence else 1000 for s in model.steps]
    cache, targets = {}, []
    for t in scenario['targets']:
        cid = t['compound_id']
        if t['net_status'] in ('exact-net-conversion-hypothesis', 'explicit-exchange-species; not a synthesis target'):
            targets.append({**t, 'diagnostic_status': 'retained-baseline-status; not rescreened'})
            continue
        if cid not in cache:
            solved = model.solve(cid, step_costs=costs)
            cache[cid] = {'compound_id': cid, **solved,
                'selected_missing_reaction_ids': sorted({s['reaction_id'] for s in solved.get('steps', [])} - evidence.keys())}
            if solved.get('steps'):
                print(cid, 'gaps', len(cache[cid]['selected_missing_reaction_ids']), flush=True)
        targets.append({**t, 'diagnostic_status': cache[cid]['status'], 'diagnostic_compound_id': cid})
    uses = defaultdict(list)
    records = defaultdict(list)
    for t in targets:
        if 'diagnostic_compound_id' in t:
            records[t['compound_id']].append(t['cannabisdb_id'])
    for r in cache.values():
        for step in r.get('steps', []):
            if step['reaction_id'] not in evidence:
                uses[step['reaction_id']].append({'compound_id': r['compound_id'],
                    'target_ids': records[r['compound_id']], **step})
    reactions = {r['id']: r for r in network['reactions']}
    gaps = [{'reaction_id': rid, 'source_joins': reactions[rid]['sources'], 'selected_uses': use,
        'target_ids': sorted({t for u in use for t in u['target_ids']}),
        'prior_searches': [{'report': p, 'row': row} for p, search in sorted(searches.items()) for row in search['rows'] if row['reaction_id'] == rid]}
        for rid, use in sorted(uses.items(), key=lambda item: (-len(item[1]), item[0]))]
    used = {s['reaction_id'] for r in cache.values() for s in r.get('steps', [])}
    needed = set(cache) | set(exchange) | {m['compound_id'] for rid in used for side in ('left', 'right') for m in reactions[rid][side]}
    return {'schema': 'cannabis-carbon.phase1-remaining-weighted-routes.v1', 'targets': targets,
        'results': list(cache.values()), 'candidate_gaps': gaps,
        'reactions': [reactions[rid] for rid in sorted(used)], 'compounds': [compounds[c] for c in sorted(needed)],
        'forbidden_step_ids': scenario['forbidden_step_ids'], 'external_exchange_compound_ids': exchange,
        'co2_compound_id': candidate['co2_compound_id'],
        'summary': {'historical_target_records': len(targets), 'diagnosed_exact_structures': len(cache),
            'target_diagnostic_status_counts': dict(Counter(t['diagnostic_status'] for t in targets)),
            'selected_gap_equations': len(gaps), 'gaps_without_prior_search': sum(not g['prior_searches'] for g in gaps)},
        'claim_boundary': 'All historical target records retained. Every target lacking a baseline net certificate or explicit exchange status is diagnosed against the full balanced catalog, preserving eight reverse exclusions. Candidate-linked directed extents cost 1; others cost 1000. This heuristic is not minimum gap cardinality or biological likelihood. Selected equations are alternative-route hypotheses, not proven necessary gaps. Only CO2 supplies net carbon; regenerated pools are allowed. No baseline status, chemical identity, or enzyme evidence is upgraded; missing annotation is not biological absence. Atom tracing remains deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-thiolase-candidate-net')]
    search_paths = sorted(Path('data/reports').glob('phase1-*search.json'))
    inputs = [json.loads(p.read_text()) for p in paths]
    report = build(*inputs, {str(p): json.loads(p.read_text()) for p in search_paths})
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*paths, *search_paths]}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-remaining-weighted-routes.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('target', 'targets', 'cannabisdb_id'), ('result', 'results', 'compound_id'),
        ('gap', 'candidate_gaps', 'reaction_id'), ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, id_key in groups:
        rows.extend((kind, r[id_key], r) for r in report[key])
    count = write_rows(rows, sha, Path('data/derived/phase1-remaining-weighted-routes.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
