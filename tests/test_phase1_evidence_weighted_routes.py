import hashlib
import json
from fractions import Fraction
from pathlib import Path
import pytest

from cannabis_carbon.phase1_net_flux import NetModel, exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def test_positive_costs_can_select_longer_route_without_changing_net_constraints():
    def r(rid, a, b):
        return {'id': rid, 'left': [{'compound_id': a, 'coefficient': 1}], 'right': [{'compound_id': b, 'coefficient': 1}]}
    model = NetModel([r('direct', 'CO2', 'B'), r('one', 'CO2', 'A'), r('two', 'A', 'B')], ['CO2'])
    assert {s['reaction_id'] for s in model.solve('B')['steps']} == {'direct'}
    result = model.solve('B', step_costs=[10 if s['reaction_id'] == 'direct' else 1 for s in model.steps])
    assert {s['reaction_id'] for s in result['steps']} == {'one', 'two'}
    assert result['external_net_consumption'] == {'CO2': '1'}
    for costs in ([1], [0] * 6, [-1] * 6, [float('inf')] * 6, [float('nan')] * 6):
        with pytest.raises(ValueError):
            model.solve('B', step_costs=costs)


def test_weighted_witnesses_replay_exactly_without_promoting_gaps():
    def read(name):
        return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())
    report = read('phase1-evidence-weighted-routes')
    candidate = read('phase1-synthase-candidate-net')
    network = read('phase1-full-balanced-network')
    compounds = {c['id']: c for c in network['compounds']}
    steps = {s['id']: s for s in orientations(network['reactions'])}
    evidence = candidate['candidate_reaction_evidence_ids']
    exchanges = set(report['external_exchange_compound_ids'])
    assert len(report['results']) == 12
    assert report['forbidden_step_ids'] == candidate['scenarios'][1]['forbidden_step_ids']
    assert {c for c in exchanges if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    for row in report['results']:
        if row['status'] != 'exact-net-conversion-hypothesis':
            continue
        assert not {s['step_id'] for s in row['steps']} & set(report['forbidden_step_ids'])
        net = exact_net([steps[s['step_id']] for s in row['steps']], [s['extent'] for s in row['steps']])
        assert net[row['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == row['net_exports']
        incoming = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
        outgoing = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
        assert incoming == outgoing > 0
        missing = sum((Fraction(s['extent']) for s in row['steps'] if s['reaction_id'] not in evidence), Fraction())
        total = sum((Fraction(s['extent']) for s in row['steps']), Fraction())
        assert Fraction(row['weighted_extent']) == total + (row['unsupported_step_cost'] - 1) * missing
        assert set(row['selected_missing_reaction_ids']) == {s['reaction_id'] for s in row['steps']} - evidence.keys()
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
