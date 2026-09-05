import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import pytest
from cannabis_carbon.phase1_routes import route, replay, build
from cannabis_carbon.phase1_scope import expand


def reaction(rid, inputs, outputs):
    return {'id': rid, 'left': [{'compound_id': c, 'coefficient': n} for c, n in inputs],
            'right': [{'compound_id': c, 'coefficient': n} for c, n in outputs]}


def test_shared_prerequisites_coproducts_and_exact_fractional_extents():
    rs = [reaction('one', [('A', 3)], [('B', 2), ('C', 1)]),
          reaction('two', [('B', 1), ('C', 2)], [('D', 3)])]
    scope = expand(rs, {'A'})
    scenario = {'id': 'synthetic', 'seed_compound_ids': ['A'], 'witnesses': scope['witnesses']}
    target = {'compound_id': 'D', 'label': 'D', 'cannabisdb_id': 'D'}
    r = route(target, scenario, {s['id']: s for s in rs})
    assert r['seed_amounts'] == {'A': '2'}
    assert [s['extent'] for s in r['steps']] == ['2/3', '1/3']
    assert r['final_inventory'] == {'B': '1', 'D': '1'}
    assert r['first_missing_enzyme_step_id'] == r['steps'][0]['id']
    scenario['witnesses']['B']['required_inputs'] = []
    with pytest.raises(ValueError, match='equation mismatch'):
        route(target, scenario, {s['id']: s for s in rs})


def test_replay_does_not_bootstrap_unchanged_catalyst_or_accept_insufficient_seed():
    step = {'extent': '1', 'required_inputs': [{'compound_id': 'A', 'coefficient': 1}, {'compound_id': 'B', 'coefficient': 1}],
            'outputs': [{'compound_id': 'C', 'coefficient': 1}, {'compound_id': 'B', 'coefficient': 1}]}
    with pytest.raises(ValueError, match='Insufficient'):
        replay([step], {'A': '1'})
    with pytest.raises(ValueError, match='Insufficient'):
        replay([step], {'A': '1/2', 'B': '1'})
    assert replay([step], {'A': '1', 'B': '1'})['C'] == 1


def test_all_published_routes_preserve_equations_provenance_and_finite_carbon_inventory():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-route-certificates.json'
    report = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-route-certificates.json').read_bytes()
    for source, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    network = json.loads((root / 'data/reports/phase1-full-balanced-network.json').read_text())
    scope = json.loads((root / 'data/reports/phase1-all-reactants-scope.json').read_text())
    rebuilt = build(network, scope)
    assert rebuilt == {k: v for k, v in report.items() if k != 'source_sha256'}
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    scenarios = {s['id']: s for s in scope['scenarios']}
    assert len(report['targets']) == sum(len(s['targets']) for s in scope['scenarios'])
    for route_record in report['routes']:
        scenario = scenarios[route_record['scenario_id']]
        inventory = defaultdict(Fraction, {c: Fraction(n) for c, n in route_record['seed_amounts'].items()})
        initial = dict(inventory)
        assert set(initial) <= set(scenario['seed_compound_ids'])
        assert all(compounds[c]['smiles'] == 'O=C=O' for c in initial if compounds[c]['carbon_count'])
        for step in route_record['steps']:
            reaction = reactions[step['reaction_id']]
            side = 'left' if step['direction_mode'] == 'hypothetical-left-to-right' else 'right'
            assert step['required_inputs'] == reaction[side]
            assert step['outputs'] == reaction['right' if side == 'left' else 'left']
            assert step['enzyme_evidence_ids'] == reaction['enzyme_evidence_ids']
            extent = Fraction(step['extent'])
            assert extent > 0
            for m in step['required_inputs']:
                inventory[m['compound_id']] -= extent * m['coefficient']
                assert inventory[m['compound_id']] >= 0
            for m in step['outputs']:
                inventory[m['compound_id']] += extent * m['coefficient']
        assert {c: str(n) for c, n in inventory.items() if n} == route_record['final_inventory']
        assert inventory[route_record['compound_id']] >= 1
        for property_name in ('carbon_count', 'formal_charge'):
            assert sum(n * compounds[c][property_name] for c, n in initial.items()) == sum(n * compounds[c][property_name] for c, n in inventory.items())
    for t in report['targets']:
        if t['route_index'] is not None:
            assert report['routes'][t['route_index']]['cannabisdb_id'] == t['cannabisdb_id']
        else:
            assert t['status'] in ('blocked', 'explicit-seed')
