"""Evidence-cost sensitivity for cannabinoid chemical witnesses, not biology."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows


def build(network, candidate, precursor, searches):
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    evidence = candidate['candidate_reaction_evidence_ids']
    exchange = candidate['external_exchange_compound_ids']
    if {c for c in exchange if compounds[c]['carbon_count']} != {candidate['co2_compound_id']}:
        raise ValueError('Organic carbon exchange forbidden')
    forbidden = candidate['scenarios'][1]['forbidden_step_ids']
    model = NetModel(network['reactions'], exchange, forbidden)
    old = {r['compound_id']: r for r in precursor['results'] if r['scenario_id'] == 'full-catalog-chemistry:eight-reverse-steps-forbidden'}
    results = []
    for penalty in (10, 1000):
        costs = [1 if s['reaction_id'] in evidence else penalty for s in model.steps]
        for probe in precursor['probes']:
            cid = probe['id']
            if cid in exchange:
                continue
            solved = model.solve(cid, step_costs=costs)
            missing = sorted({s['reaction_id'] for s in solved.get('steps', [])} - evidence.keys())
            missing_extent = sum((Fraction(s['extent']) for s in solved.get('steps', []) if s['reaction_id'] not in evidence), Fraction())
            total_extent = sum((Fraction(s['extent']) for s in solved.get('steps', [])), Fraction())
            results.append({'id': str(penalty) + ':' + cid, 'compound_id': cid, 'unsupported_step_cost': penalty,
                **solved, 'selected_missing_reaction_ids': missing,
                'missing_directed_extent': str(missing_extent), 'total_directed_extent': str(total_extent),
                'weighted_extent': str(total_extent + (penalty - 1) * missing_extent),
                'previous_missing_reaction_ids': sorted({s['reaction_id'] for s in old[cid].get('steps', [])} - evidence.keys())})
            print(penalty, cid, solved['status'], 'missing', len(missing), flush=True)
    used = {s['reaction_id'] for r in results for s in r.get('steps', [])}
    gaps = [{'reaction_id': rid, 'source_joins': reactions[rid]['sources'],
        'selected_uses': [{'probe_id': r['id'], **s} for r in results for s in r.get('steps', []) if s['reaction_id'] == rid],
        'prior_searches': [{'report': path, 'row': row} for path, report in sorted(searches.items()) for row in report['rows'] if row['reaction_id'] == rid]}
        for rid in sorted(used - evidence.keys())]
    cids = set(exchange) | {r['compound_id'] for r in results} | {m['compound_id'] for rid in used for side in ('left', 'right') for m in reactions[rid][side]}
    return {'schema': 'cannabis-carbon.phase1-evidence-weighted-routes.v1', 'results': results,
        'candidate_gaps': gaps, 'reactions': [reactions[rid] for rid in sorted(used)],
        'compounds': [compounds[c] for c in sorted(cids)], 'forbidden_step_ids': forbidden,
        'external_exchange_compound_ids': exchange, 'co2_compound_id': candidate['co2_compound_id'],
        'cost_definition': 'Candidate-linked directed extents cost 1; every other directed extent costs 10 or 1000 in separate analyst-selected sensitivities. No evidence class is upgraded.',
        'summary': {'probe_sensitivity_tests': len(results), 'exact_certificates': sum(r['status'] == 'exact-net-conversion-hypothesis' for r in results),
            'selected_gap_equations': len(gaps), 'gaps_without_prior_search_row': sum(not g['prior_searches'] for g in gaps)},
        'claim_boundary': 'Weighted reaction extents are a prioritization heuristic, not a minimum count of missing enzymes or a likelihood of plant production. Equation scaling affects costs. No feasibility, charge, isotope, stereochemical or exchange constraint is relaxed. Eight reverse exclusions are preserved; other directions remain hypothetical. Net certificates allow regenerated pre-existing pools and do not establish startup, compartmentation or flux. Atom tracing and whole-metabolome completeness claims remain deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-synthase-candidate-net', 'phase1-synthase-precursor-audit')]
    search_paths = sorted(Path('data/reports').glob('phase1-*search.json'))
    report = build(*(json.loads(p.read_text()) for p in paths), {str(p): json.loads(p.read_text()) for p in search_paths})
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*paths, *search_paths]}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-evidence-weighted-routes.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('result', 'results', 'id'), ('gap', 'candidate_gaps', 'reaction_id'), ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-evidence-weighted-routes.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}))


if __name__ == '__main__':
    run()
