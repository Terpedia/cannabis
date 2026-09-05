import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import pytest
from cannabis_carbon.phase1_net_flux import NetModel, exact_net
from cannabis_carbon.phase1_scope import expand, orientations


def reaction(rid, left, right):
    return {'id': rid, 'left': [{'compound_id': c, 'coefficient': n} for c, n in left],
            'right': [{'compound_id': c, 'coefficient': n} for c, n in right]}


def test_regenerated_pool_allows_net_conversion_but_not_zero_pool_startup():
    reactions = [reaction('fix', [('CO2', 1), ('pool', 1)], [('intermediate', 1)]),
                 reaction('release', [('intermediate', 1)], [('pool', 1), ('target', 1)])]
    assert expand(reactions, {'CO2'})['available'] == {'CO2'}
    result = NetModel(reactions, {'CO2'}).solve('target')
    assert result['status'] == 'exact-net-conversion-hypothesis'
    assert result['external_net_consumption'] == {'CO2': '1'}
    assert result['net_exports'] == {'target': '1'}
    assert result['zero_net_internal_participants'] == ['intermediate', 'pool']


def test_unregenerated_currency_and_missing_carbon_cannot_supply_target():
    reactions = [reaction('consume', [('CO2', 1), ('cofactor', 1)], [('target', 1)])]
    model = NetModel(reactions, {'CO2'})
    assert model.solve('target')['status'] == 'solver-reported-infeasible'
    assert model.solve('absent')['status'] == 'no-net-producing-candidate-equation'
    assert model.solve('CO2')['status'].startswith('explicit-exchange')
    balanced_cycle = [reaction('one', [('A', 1)], [('B', 1)]), reaction('two', [('B', 1)], [('A', 1)])]
    assert NetModel(balanced_cycle, {'CO2'}).solve('A')['status'] == 'solver-reported-infeasible'


def test_all_equation_coefficients_and_byproducts_are_retained():
    reactions = [reaction('split', [('CO2', 3)], [('target', 2), ('byproduct', 1)])]
    result = NetModel(reactions, {'CO2'}).solve('target')
    assert result['steps'][0]['extent'] == '1/2'
    assert result['external_net_consumption'] == {'CO2': '3/2'}
    assert result['net_exports'] == {'byproduct': '1/2', 'target': '1'}
    with pytest.raises(ValueError, match='Every step'):
        exact_net(orientations(reactions), ['1'])


def test_all_published_net_certificates_have_exact_closed_internal_balances():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-candidate-net-flux.json'
    report = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-candidate-net-flux.json').read_bytes()
    for source, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    network = json.loads((root / 'data/reports/phase1-full-balanced-network.json').read_text())
    startup = json.loads((root / 'data/reports/phase1-candidate-scope.json').read_text())
    scenario = next(s for s in startup['scenarios'] if s['id'] == 'CO2-plus-all-carbon-free-species')
    exchanges = set(report['external_exchange_compound_ids'])
    assert exchanges == set(scenario['seed_compound_ids'])
    assert [t['cannabisdb_id'] for t in report['targets']] == [t['cannabisdb_id'] for t in scenario['targets']]
    original = {r['id']: r for r in network['reactions']}
    reactions = {r['id']: r for r in report['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    assert [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] == ['O=C=O']
    for rid, r in reactions.items():
        assert r['enzyme_evidence_ids'] == startup['candidate_reaction_evidence_ids'][rid]
        assert {k: v for k, v in r.items() if k != 'enzyme_evidence_ids'} == {k: v for k, v in original[rid].items() if k != 'enzyme_evidence_ids'}
    certs = {c['compound_id']: c for c in report['certificates']}
    for cert in certs.values():
        # Independent exact replay, without trusting linprog or certificate net fields.
        net = defaultdict(Fraction)
        for step in cert['steps']:
            r = reactions[step['reaction_id']]
            assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
            left = 'left' if step['direction_mode'] == 'hypothetical-left-to-right' else 'right'
            right = 'right' if left == 'left' else 'left'
            extent = Fraction(step['extent'])
            assert extent > 0
            for side, sign in [(left, -1), (right, 1)]:
                for m in r[side]:
                    net[m['compound_id']] += sign * extent * m['coefficient']
        assert net[cert['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == cert['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == cert['net_exports']
        assert sorted(c for c, n in net.items() if not n and c not in exchanges) == cert['zero_net_internal_participants']
        for prop in ('carbon_count', 'formal_charge'):
            assert sum(n * compounds[c][prop] for c, n in net.items()) == 0
        assert -net[report['co2_compound_id']] == Fraction(cert['net_carbon_in']) == Fraction(cert['net_carbon_out'])
    for old, target in zip(scenario['targets'], report['targets']):
        assert target['startup_status'] == old['status']
        assert (target['certificate_compound_id'] in certs) == (target['net_status'] == 'exact-net-conversion-hypothesis')
    # Independently repeat every numerically infeasible solve; this remains
    # numerical evidence, not an exact infeasibility or biological proof.
    model = NetModel([r for r in original.values() if r['id'] in startup['candidate_reaction_evidence_ids']], exchanges)
    for cid in {t['compound_id'] for t in report['targets'] if t['net_status'] == 'solver-reported-infeasible'}:
        assert model.solve(cid)['status'] == 'solver-reported-infeasible'
