"""Candidate-linked net CO2 conversion with explicit conserved-pool assumptions."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from .phase1_scope import orientations, write_rows


def exact_net(steps, extents):
    if len(steps) != len(extents):
        raise ValueError('Every step needs an extent')
    net = defaultdict(Fraction)
    for step, extent in zip(steps, extents):
        amount = Fraction(extent)
        if amount < 0:
            raise ValueError('Negative directed extent')
        for side, sign in [('required_inputs', -1), ('outputs', 1)]:
            for m in step[side]:
                if m['coefficient'] <= 0:
                    raise ValueError('Nonpositive reaction coefficient')
                net[m['compound_id']] += sign * amount * m['coefficient']
    return dict(net)


class NetModel:
    def __init__(self, reactions, exchange_ids, forbidden_step_ids=()):
        self.steps = orientations(reactions)
        forbidden = set(forbidden_step_ids)
        if not forbidden <= {s['id'] for s in self.steps}:
            raise ValueError('Unknown forbidden directed step')
        self.steps = [s for s in self.steps if s['id'] not in forbidden]
        self.exchange_ids = set(exchange_ids)
        self.internal_ids = sorted({m['compound_id'] for s in self.steps for side in ('required_inputs', 'outputs') for m in s[side]} - self.exchange_ids)
        self.index = {c: i for i, c in enumerate(self.internal_ids)}
        row, col, values = [], [], []
        self.producible = set()
        for j, step in enumerate(self.steps):
            net = exact_net([step], ['1'])
            self.producible.update(c for c, n in net.items() if n > 0)
            for c, amount in net.items():
                if c in self.index and amount:
                    row.append(self.index[c]); col.append(j); values.append(-float(amount))
        self.matrix = coo_matrix((values, (row, col)), shape=(len(self.index), len(self.steps))).tocsr()

    def solve(self, target):
        if target in self.exchange_ids:
            return {'status': 'explicit-exchange-species; not a synthesis target'}
        if target not in self.producible or target not in self.index:
            return {'status': 'no-net-producing-candidate-equation'}
        bound = np.zeros(len(self.index)); bound[self.index[target]] = -1
        result = linprog(np.ones(len(self.steps)), A_ub=self.matrix, b_ub=bound,
            bounds=(0, None), method='highs', options={'time_limit': 30,
                'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
        record = {'solver_status': int(result.status), 'solver_message': result.message}
        if not result.success:
            return {**record, 'status': 'solver-reported-infeasible' if result.status == 2 else 'solver-incomplete-or-failed'}
        extents = [Fraction(float(x)).limit_denominator(1000000) for x in result.x]
        net = exact_net(self.steps, extents)
        if net.get(target, 0) < 1 or any(n < 0 for c, n in net.items() if c not in self.exchange_ids):
            return {**record, 'status': 'numerical-solution-failed-exact-validation'}
        used = [(s, n) for s, n in zip(self.steps, extents) if n]
        participants = {m['compound_id'] for s, _ in used for side in ('required_inputs', 'outputs') for m in s[side]}
        return {**record, 'status': 'exact-net-conversion-hypothesis',
            'steps': [{'step_id': s['id'], 'reaction_id': s['reaction_id'], 'direction_mode': s['direction_mode'], 'extent': str(n)} for s, n in used],
            'target_amount': str(net[target]),
            'external_net_consumption': {c: str(-n) for c, n in sorted(net.items()) if n < 0},
            'net_exports': {c: str(n) for c, n in sorted(net.items()) if n > 0},
            'zero_net_internal_participants': sorted(c for c in participants if c not in self.exchange_ids and not net[c]),
            'selection_boundary': 'One minimum-total-directed-extent numerical solution, rationally reconstructed and checked exactly. Not a shortest route, unique solution, functional ranking or thermodynamic model.',
            'startup_status': 'Pre-existing internal pools may be required; synthesis and minimum pool sizes are not established by this certificate.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ['phase1-full-balanced-network', 'phase1-route-enzyme-overlay', 'phase1-candidate-scope']]
    network, overlay, startup = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for p in paths[:2]:
        if startup['source_sha256'][str(p)] != hashes[str(p)]:
            raise ValueError('Startup source mismatch')
    evidence = startup['candidate_reaction_evidence_ids']
    reactions = [r for r in network['reactions'] if r['id'] in evidence]
    compounds = {c['id']: c for c in network['compounds']}
    scenario = next(s for s in startup['scenarios'] if s['id'] == 'CO2-plus-all-carbon-free-species')
    exchanges = set(scenario['seed_compound_ids'])
    carbon_exchanges = [c for c in exchanges if compounds[c]['carbon_count']]
    if len(carbon_exchanges) != 1 or compounds[carbon_exchanges[0]]['smiles'] != 'O=C=O':
        raise ValueError('CO2 must be the sole carbon exchange')
    co2 = carbon_exchanges[0]
    model = NetModel(reactions, exchanges)
    targets, certificates = [], {}
    cache = {}
    for i, target in enumerate(scenario['targets'], 1):
        cid = target['compound_id']
        if cid not in cache:
            cache[cid] = model.solve(cid)
        result = cache[cid]
        if result['status'] == 'exact-net-conversion-hypothesis':
            if cid not in certificates:
                carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
                carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
                if carbon_in != carbon_out or carbon_in <= 0 or Fraction(result['external_net_consumption'].get(co2, '0')) != carbon_in:
                    raise ValueError('Exact carbon exchange validation failed')
                certificates[cid] = {'compound_id': cid, **result, 'net_carbon_in': str(carbon_in), 'net_carbon_out': str(carbon_out)}
        targets.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'carbon_count')} |
            {'startup_status': target['status'], 'net_status': result['status'],
             'certificate_compound_id': cid if cid in certificates else None} |
            {k: result[k] for k in ('solver_status', 'solver_message') if k in result})
        if i % 500 == 0:
            print(f'Net conversion targets {i}/{len(scenario["targets"])}; certificates {len(certificates)}', flush=True)
    used_ids = {s['reaction_id'] for cert in certificates.values() for s in cert['steps']}
    used_compounds = {m['compound_id'] for r in reactions if r['id'] in used_ids for side in ('left', 'right') for m in r[side]} | exchanges
    report = {'schema': 'cannabis-carbon.phase1-candidate-net-flux.v1',
        'summary': {'target_records': len(targets), 'candidate_equations': len(reactions), 'exact_structure_certificates': len(certificates),
            'target_status_counts': dict(Counter(t['net_status'] for t in targets)),
            'net_feasible_targets_not_startup_reachable': sum(t['net_status'] == 'exact-net-conversion-hypothesis' and t['startup_status'] == 'blocked' for t in targets)},
        'targets': targets, 'certificates': list(certificates.values()),
        'reactions': [{**r, 'enzyme_evidence_ids': evidence[r['id']]} for r in reactions if r['id'] in used_ids],
        'compounds': [compounds[c] for c in sorted(used_compounds)],
        'external_exchange_compound_ids': sorted(exchanges), 'co2_compound_id': co2,
        'source_sha256': hashes, 'scipy_version': scipy.__version__,
        'solver': {'method': 'highs', 'objective': 'minimize sum of nonnegative directed reaction extents',
            'constraints': 'Every nonexchange species has net production >= 0; target net production >= 1. All positive net products are explicitly exported. Only CO2 and listed carbon-free species may have negative net production.',
            'time_limit_seconds_per_target': 30, 'rational_reconstruction_max_denominator': 1000000,
            'documentation': 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html'},
        'claim_boundary': 'Net stoichiometric hypotheses with pre-existing internal pools allowed, not zero-pool startup certificates. No internal compound is depleted overall and all net carbon input is CO2. Full pool regeneration does not establish pool origin or atom-wise provenance. Both reaction directions and all carbon-free exchanges are permissive assumptions; energy/redox feasibility, thermodynamic loops, compartments, enzyme activity, specificity and physiological direction are unverified. Numerical infeasibility is solver-reported, not proof of biological absence. Atom tracing remains deferred.'}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    for folder in ('data/reports', 'docs/data'):
        Path(folder, 'phase1-candidate-net-flux.json').write_text(payload)
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in ('targets', 'certificates', 'reactions', 'compounds')})]
    for kind, collection, key in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id'), ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id')]:
        rows.extend((kind, row[key], row) for row in report[collection])
    write_rows(rows, hashlib.sha256(payload.encode()).hexdigest(), Path('data/reports/phase1-candidate-net-flux.ndjson'))
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
