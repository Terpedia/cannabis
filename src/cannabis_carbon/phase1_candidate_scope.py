"""Evidence-constrained scope and explicit single-gap rescue hypotheses."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_scope import expand, orientations, write_rows


def build(network, overlay, scope, search):
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    added = {e['reaction_id']: e for e in overlay['enzyme_evidence']}
    if not added.keys() <= reactions.keys():
        raise ValueError('Overlay reaction missing from network')
    evidence = {rid: r['enzyme_evidence_ids'] + ([added[rid]['id']] if rid in added else []) for rid, r in reactions.items()}
    supported = [r for rid, r in reactions.items() if evidence[rid]]
    unsupported = [r for rid, r in reactions.items() if not evidence[rid]]
    search_rows = {r['reaction_id']: r for r in search['rows']}
    scenarios = []
    for original in scope['scenarios']:
        seeds = set(original['seed_compound_ids'])
        result = expand(supported, seeds)
        available = result['available']
        targets = [{k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'carbon_count')} |
            {'status': 'explicit-seed' if t['compound_id'] in seeds else 'candidate-scope-reachable' if t['compound_id'] in available else 'blocked',
             'chemistry_only_status': t['status']} for t in original['targets']]
        frontiers = []
        for step in orientations(unsupported):
            if not all(m['compound_id'] in available for m in step['required_inputs']):
                continue
            new_outputs = sorted({m['compound_id'] for m in step['outputs']} - available)
            if not new_outputs:
                continue
            reaction = reactions[step['reaction_id']]
            # The single exceptional equation is admitted in both hypothetical
            # directions, matching the existing upper-bound direction scenario.
            rescued = expand(supported + [reaction], seeds)
            new_compounds = rescued['available'] - available
            new_targets = [t['cannabisdb_id'] for t in targets if t['compound_id'] in new_compounds]
            entry = search_rows.get(reaction['id'])
            frontiers.append({**step, 'status': 'single-gap-counterfactual; not candidate-supported',
                'sources': reaction['sources'], 'new_output_compound_ids': new_outputs,
                'prior_search_status': entry['search_status'] if entry else 'not-in-selected-route-protein-screen',
                'rescued_target_ids': new_targets, 'rescued_compound_ids': sorted(new_compounds),
                'rescued_witnesses': rescued['witnesses'],
                'rescued_target_count': len(new_targets),
                'next_test': 'First curate whether the full equation is enzyme-catalyzed, spontaneous, or only a transformation rule, and establish a viable direction. If enzymatic, resolve exact reaction references and investigate Cannabis candidates before biochemical assays. Do not promote this counterfactual to enzyme evidence.'})
        frontiers.sort(key=lambda f: (-f['rescued_target_count'], -len(f['rescued_compound_ids']), f['id']))
        scenarios.append({'id': original['id'], 'seed_compound_ids': sorted(seeds),
            'seed_boundary': original['seed_boundary'], 'direction_boundary': original['direction_boundary'],
            'summary': {'candidate_equations': len(supported), 'excluded_equations_without_candidate_evidence': len(unsupported),
                'available_compounds': len(available), 'newly_available_compounds': len(available - seeds),
                'target_status_counts': dict(Counter(t['status'] for t in targets)),
                'carbon_bearing_target_status_counts': dict(Counter(t['status'] for t in targets if t['carbon_count'])),
                'candidate_reactant_denominator': len({m['compound_id'] for r in supported for side in ('left', 'right') for m in r[side]}),
                'full_network_reactant_denominator': len({m['compound_id'] for r in reactions.values() for side in ('left', 'right') for m in r[side]}),
                'single_gap_frontier_steps': len(frontiers),
                'frontier_steps_rescuing_targets': sum(bool(f['rescued_target_ids']) for f in frontiers)},
            'targets': targets, 'witnesses': result['witnesses'], 'frontiers': frontiers})
    return {'schema': 'cannabis-carbon.phase1-candidate-scope.v1', 'scenarios': scenarios,
        'candidate_reaction_evidence_ids': {rid: ids for rid, ids in evidence.items() if ids},
        'compound_details': [compounds[c] for c in sorted({c for s in scenarios for c in s['witnesses']} |
            {c for s in scenarios for f in s['frontiers'] for c in f['rescued_compound_ids']})],
        'claim_boundary': 'Candidate annotations/homology are not direct enzyme activity. Both directions and the highly permissive carbon-free reservoir remain hypothetical. Qualitative zero-organic-inventory startup is not steady-state assimilation: pre-existing cofactors and autocatalytic pools are not silently seeded. Failure to start does not establish inability of a living plant to assimilate CO2. Single-gap rescues remain unsupported hypotheses, not evidence. Atom tracing remains deferred.'}


def run():
    names = ['phase1-full-balanced-network', 'phase1-route-enzyme-overlay', 'phase1-all-reactants-scope', 'phase1-route-protein-search']
    paths = [Path('data/reports', name + '.json') for name in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    network, overlay, scope, search = inputs
    if scope['source_network_sha256'] != hashes[str(paths[0])]:
        raise ValueError('Scope network checksum mismatch')
    for p in (paths[0], paths[3]):
        if overlay['source_sha256'][str(p)] != hashes[str(p)]:
            raise ValueError('Overlay source checksum mismatch')
    report = build(*inputs)
    report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    for folder in ('data/reports', 'docs/data'):
        Path(folder, 'phase1-candidate-scope.json').write_text(payload)
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in ('scenarios', 'compound_details')})]
    rows.extend(('compound', c['id'], c) for c in report['compound_details'])
    for scenario in report['scenarios']:
        sid = scenario['id']
        rows.append(('scenario', sid, {k: v for k, v in scenario.items() if k not in ('targets', 'witnesses', 'frontiers')}))
        rows.extend(('target', sid + ':' + t['cannabisdb_id'], {'scenario_id': sid, **t}) for t in scenario['targets'])
        rows.extend(('witness', sid + ':' + c, {'scenario_id': sid, 'compound_id': c, **w}) for c, w in scenario['witnesses'].items())
        rows.extend(('frontier', sid + ':' + f['id'], {'scenario_id': sid, **f}) for f in scenario['frontiers'])
    write_rows(rows, hashlib.sha256(payload.encode()).hexdigest(), Path('data/reports/phase1-candidate-scope.ndjson'))
    print(json.dumps({s['id']: s['summary'] for s in report['scenarios']}), flush=True)


if __name__ == '__main__':
    run()
