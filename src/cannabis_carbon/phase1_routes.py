"""All-prerequisite route witnesses and exact rational execution certificates.

No atom mapping, thermodynamics or physiological direction is inferred.
"""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


def route(target, scenario, reactions):
    witnesses = scenario['witnesses']
    seeds = set(scenario['seed_compound_ids'])
    cid = target['compound_id']
    if cid not in witnesses or cid in seeds:
        raise ValueError('Route requires a reachable non-seed target')
    visited, steps = set(), {}

    def visit(compound):
        if compound in visited:
            return
        w = witnesses[compound]
        if compound in seeds:
            if w['level'] != 0 or w['source'] != 'explicit-seed':
                raise ValueError('Invalid seed witness')
            visited.add(compound)
            return
        r = reactions[w['reaction_id']]
        if w['direction_mode'] not in ('hypothetical-left-to-right', 'hypothetical-right-to-left'):
            raise ValueError('Unexpected direction')
        side = 'left' if w['direction_mode'] == 'hypothetical-left-to-right' else 'right'
        other = 'right' if side == 'left' else 'left'
        if w['required_inputs'] != r[side] or w['outputs'] != r[other]:
            raise ValueError('Witness equation mismatch')
        if not any(m['compound_id'] == compound for m in w['outputs']):
            raise ValueError('Witness does not produce compound')
        for m in w['required_inputs']:
            if witnesses[m['compound_id']]['level'] >= w['level']:
                raise ValueError('Cyclic or unordered prerequisite')
            visit(m['compound_id'])
        steps[w['id']] = {k: w[k] for k in ('id', 'level', 'reaction_id', 'direction_mode', 'required_inputs', 'outputs')}
        visited.add(compound)

    visit(cid)
    ordered = sorted(steps.values(), key=lambda s: (s['level'], s['id']))
    demand = defaultdict(Fraction, {cid: Fraction(1)})
    for step in reversed(ordered):
        # Credit only outputs assigned to this witness. Ignoring incidental
        # coproduct credits deliberately yields a conservative, nonminimal seed budget.
        owned = [m for m in step['outputs'] if m['compound_id'] in visited
                 and witnesses[m['compound_id']].get('id') == step['id']]
        scale = max(demand[m['compound_id']] / m['coefficient'] for m in owned)
        if scale <= 0:
            raise ValueError('Unneeded route step')
        step['extent'] = str(scale)
        for m in owned:
            demand[m['compound_id']] = Fraction(0)
        for m in step['required_inputs']:
            demand[m['compound_id']] += scale * m['coefficient']
        step['enzyme_evidence_ids'] = reactions[step['reaction_id']].get('enzyme_evidence_ids', [])
        step['enzyme_status'] = 'candidate-evidence-attached' if step['enzyme_evidence_ids'] else 'no-candidate-evidence-attached'
        step['direction_status'] = 'hypothetical; physiology unresolved'
    if any(amount and compound not in seeds for compound, amount in demand.items()):
        raise ValueError('Unresolved non-seed demand')
    seed_amounts = {c: str(n) for c, n in sorted(demand.items()) if n}
    inventory = replay(ordered, seed_amounts)
    if inventory.get(cid, Fraction(0)) < 1:
        raise ValueError('Target not delivered')
    missing = [s['id'] for s in ordered if not s['enzyme_evidence_ids']]
    return {'cannabisdb_id': target['cannabisdb_id'], 'label': target['label'], 'compound_id': cid,
        'scenario_id': scenario['id'], 'status': 'structurally_possible',
        'biological_status': 'unestablished', 'selection': 'Deterministic first scope witness; not shortest, optimal, unique or exhaustive.',
        'steps': ordered, 'seed_amounts': seed_amounts,
        'final_inventory': {c: str(n) for c, n in sorted(inventory.items()) if n},
        'target_amount_requested': '1', 'target_amount_delivered': str(inventory[cid]),
        'missing_enzyme_step_ids': missing,
        'first_missing_enzyme_step_id': missing[0] if missing else None,
        'blockers': ['physiological-direction-unresolved', 'thermodynamics-unassessed',
            'compartments-and-transport-unassessed', 'Cannabis-route-activity-unestablished',
            'nonphysiological-seed-scenario'] + (['steps-without-candidate-enzyme-evidence'] if missing else []),
        'claim_boundary': 'Exact rational amounts certify nonnegative sequential inventory for these balanced equations and explicit seeds only. They do not establish energetics, enzyme activity, flux in Cannabis, minimal requirements or atom provenance.'}


def replay(steps, seed_amounts):
    inventory = defaultdict(Fraction, {c: Fraction(n) for c, n in seed_amounts.items()})
    if any(n < 0 for n in inventory.values()):
        raise ValueError('Negative initial inventory')
    for step in steps:
        extent = Fraction(step['extent'])
        if extent <= 0:
            raise ValueError('Nonpositive extent')
        # Consume all inputs before adding any outputs, including unchanged species.
        for m in step['required_inputs']:
            amount = extent * m['coefficient']
            if m['coefficient'] <= 0 or inventory[m['compound_id']] < amount:
                raise ValueError('Insufficient reactant inventory')
            inventory[m['compound_id']] -= amount
        for m in step['outputs']:
            if m['coefficient'] <= 0:
                raise ValueError('Nonpositive coefficient')
            inventory[m['compound_id']] += extent * m['coefficient']
    return inventory


def build(network, scope):
    reactions = {r['id']: r for r in network['reactions']}
    routes, targets, gaps = [], [], {}
    for scenario in scope['scenarios']:
        for t in scenario['targets']:
            item = {k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'carbon_count', 'status')}
            item['scenario_id'] = scenario['id']
            item['route_index'] = None
            if t['status'] == 'structural-scope-reachable':
                r = route(t, scenario, reactions)
                item['route_index'] = len(routes)
                routes.append(r)
                for step in r['steps']:
                    if step['enzyme_evidence_ids']:
                        continue
                    rid = step['reaction_id']
                    gap = gaps.setdefault(rid, {'reaction_id': rid, 'target_ids': set(), 'route_indices': [],
                        'sources': reactions[rid]['sources'],
                        'next_test': 'Retrieve reaction-specific characterized reference proteins, screen the full Cannabis proteome, and assay shortlisted proteins with all listed substrates. Homology alone does not establish activity.'})
                    gap['target_ids'].add(t['cannabisdb_id'])
                    gap['route_indices'].append(item['route_index'])
            targets.append(item)
    gap_rows = []
    for gap in gaps.values():
        gap['target_ids'] = sorted(gap['target_ids'])
        gap['route_indices'] = sorted(set(gap['route_indices']))
        gap['selected_route_target_count'] = len(gap['target_ids'])
        gap_rows.append(gap)
    gap_rows.sort(key=lambda g: (-g['selected_route_target_count'], g['reaction_id']))
    return {'schema': 'cannabis-carbon.phase1-route-certificates.v1',
        'summary': {'routes': len(routes), 'target_scenario_records': len(targets),
            'unique_selected_equations': len({s['reaction_id'] for r in routes for s in r['steps']}),
            'selected_equations_without_candidate_evidence': len(gap_rows),
            'routes_with_missing_candidate_evidence': sum(bool(r['missing_enzyme_step_ids']) for r in routes),
            'max_route_steps': max((len(r['steps']) for r in routes), default=0),
            'target_status_counts': dict(Counter(t['status'] for t in targets))},
        'targets': targets, 'routes': routes, 'enzyme_gap_queue': gap_rows,
        'evidence_boundary': 'Gap priorities count selected deterministic route witnesses, not unavoidable bottlenecks across all possible pathways. Candidate evidence IDs resolve through the checksummed network provenance. No atom tracing is performed.',
        'scenario_boundaries': [{k: s[k] for k in ('id', 'seed_boundary', 'direction_boundary', 'seed_compound_ids')} for s in scope['scenarios']]}


def run():
    network_path = Path('data/reports/phase1-full-balanced-network.json')
    scope_path = Path('data/reports/phase1-all-reactants-scope.json')
    scope = json.loads(scope_path.read_text())
    digest = hashlib.sha256(network_path.read_bytes()).hexdigest()
    if digest != scope['source_network_sha256']:
        raise ValueError('Network checksum mismatch')
    report = build(json.loads(network_path.read_text()), scope)
    report['source_sha256'] = {str(network_path): digest, str(scope_path): hashlib.sha256(scope_path.read_bytes()).hexdigest()}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    for folder in ('data/reports', 'docs/data'):
        Path(folder, 'phase1-route-certificates.json').write_text(payload)
    from .phase1_scope import write_rows
    rows = [('metadata', 'routes', {k: v for k, v in report.items() if k not in ('targets', 'routes', 'enzyme_gap_queue')})]
    rows.extend(('target', t['scenario_id'] + ':' + t['cannabisdb_id'], t) for t in report['targets'])
    rows.extend(('route', r['scenario_id'] + ':' + r['cannabisdb_id'], r) for r in report['routes'])
    rows.extend(('enzyme_gap', g['reaction_id'], g) for g in report['enzyme_gap_queue'])
    write_rows(rows, hashlib.sha256(payload.encode()).hexdigest(), Path('data/reports/phase1-route-certificates.ndjson'))
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
