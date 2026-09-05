"""Probe exact synthase participants and retain missing enzyme steps in witnesses."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows


def build(network, candidates, links, searches):
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    evidence = candidates['candidate_reaction_evidence_ids']
    exchange = candidates['external_exchange_compound_ids']
    if {c for c in exchange if compounds[c]['carbon_count']} != {candidates['co2_compound_id']}:
        raise ValueError('Organic carbon exchange forbidden')
    if compounds[candidates['co2_compound_id']]['smiles'] != 'O=C=O':
        raise ValueError('CO2 identity mismatch')
    results, scenarios = [], []
    for level in ('candidate', 'full-catalog-chemistry'):
        selected = [r for r in network['reactions'] if level != 'candidate' or r['id'] in evidence]
        for scenario in candidates['scenarios']:
            sid = level + ':' + scenario['id']
            model = NetModel(selected, exchange, scenario['forbidden_step_ids'])
            scenarios.append({'id': sid, 'equation_count': len(selected), 'forbidden_step_ids': scenario['forbidden_step_ids']})
            for compound in links['compounds']:
                cid = compound['id']
                result = {'id': sid + ':' + cid, 'scenario_id': sid, 'compound_id': cid, **model.solve(cid)}
                results.append(result)
                print(sid, cid, result['status'], flush=True)
    used = {s['reaction_id'] for r in results for s in r.get('steps', [])}
    gaps = []
    for rid in sorted(used - evidence.keys()):
        prior = [{'report': path, 'row': row} for path, report in sorted(searches.items()) for row in report['rows'] if row['reaction_id'] == rid]
        gaps.append({'reaction_id': rid, 'source_joins': reactions[rid]['sources'],
            'selected_uses': [{'probe_id': r['id'], **s} for r in results for s in r.get('steps', []) if s['reaction_id'] == rid],
            'prior_searches': prior, 'next_action': 'Review existing evidence and reference gaps before rescreening.' if prior else 'Discover exact reference enzymes and screen the pinned whole proteome.'})
    needed = set(exchange) | {c['id'] for c in links['compounds']} | {m['compound_id'] for rid in used for side in ('left', 'right') for m in reactions[rid][side]}
    return {'schema': 'cannabis-carbon.phase1-synthase-precursor-audit.v1', 'results': results,
        'scenarios': scenarios, 'candidate_gaps': gaps, 'probes': links['compounds'],
        'external_exchange_compound_ids': exchange, 'co2_compound_id': candidates['co2_compound_id'],
        'reactions': [reactions[rid] for rid in sorted(used)], 'compounds': [compounds[c] for c in sorted(needed)],
        'summary': {'exact_participant_probes': len(links['compounds']), 'scenario_tests': len(results),
            'status_counts_by_scenario': {s['id']: dict(Counter(r['status'] for r in results if r['scenario_id'] == s['id'])) for s in scenarios},
            'selected_missing_candidate_equations': len(gaps), 'gaps_without_prior_search_row': sum(not g['prior_searches'] for g in gaps)},
        'claim_boundary': 'Exact participant probes retain neutral/charged and unspecified/specified geometry separately. Full-catalog witnesses are chemical possibilities without complete Cannabis enzyme support; chosen gaps are not necessary or sufficient in every alternative route. Each probe is solved independently, not a joint supply or startup model. Reference directions remain a sensitivity and other orientations hypothetical. Only CO2 supplies net carbon. No historical target, enzyme status or completeness count is promoted. Atom tracing is deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-synthase-candidate-net', 'phase1-synthase-reaction-links')]
    search_paths = sorted(Path('data/reports').glob('phase1-*search.json'))
    inputs = [json.loads(p.read_text()) for p in paths]
    report = build(*inputs, {str(p): json.loads(p.read_text()) for p in search_paths})
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [*paths, *search_paths]}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-synthase-precursor-audit.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('result', 'results', 'id'), ('gap', 'candidate_gaps', 'reaction_id'), ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-synthase-precursor-audit.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}))


if __name__ == '__main__':
    run()
