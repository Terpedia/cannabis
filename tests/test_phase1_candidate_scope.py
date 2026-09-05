import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_candidate_scope import build
from cannabis_carbon.phase1_scope import expand


def test_single_gap_rescue_requires_all_inputs_and_does_not_promote_evidence():
    def r(rid, a, b, ev):
        return {'id': rid, 'left': [{'compound_id': c, 'coefficient': 1} for c in a],
                'right': [{'compound_id': b, 'coefficient': 1}], 'enzyme_evidence_ids': ev, 'sources': []}
    reactions = [r('supported', ['A'], 'B', ['e']), r('gap', ['B'], 'C', []),
                 r('downstream', ['C'], 'D', ['e']), r('blocked', ['D', 'X'], 'E', ['e'])]
    targets = [{'cannabisdb_id': c, 'compound_id': c, 'carbon_count': 1, 'label': c, 'status': 'structural-scope-reachable'} for c in ['D', 'E']]
    network = {'compounds': [{'id': c} for c in 'ABCDE'], 'reactions': reactions}
    scope = {'scenarios': [{'id': 'test', 'seed_compound_ids': ['A'], 'seed_boundary': 'test', 'direction_boundary': 'test', 'targets': targets}]}
    report = build(network, {'enzyme_evidence': []}, scope, {'rows': []})
    scenario = report['scenarios'][0]
    assert set(scenario['witnesses']) == {'A', 'B'}
    assert len(scenario['frontiers']) == 1
    frontier = scenario['frontiers'][0]
    assert frontier['rescued_target_ids'] == ['D']
    assert set(frontier['rescued_witnesses']) == {'A', 'B', 'C', 'D'}
    assert 'gap' not in report['candidate_reaction_evidence_ids']
    assert all(t['status'] == 'blocked' for t in scenario['targets'])


def test_published_candidate_scope_and_every_rescue_are_replayable():
    root = Path(__file__).resolve().parents[1]
    def read(name):
        return json.loads((root / 'data/reports' / (name + '.json')).read_text())
    report = read('phase1-candidate-scope')
    assert (root / 'data/reports/phase1-candidate-scope.json').read_bytes() == (root / 'docs/data/phase1-candidate-scope.json').read_bytes()
    for source, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    network = read('phase1-full-balanced-network')
    overlay = read('phase1-route-enzyme-overlay')
    original = read('phase1-all-reactants-scope')
    search = read('phase1-route-protein-search')
    assert build(network, overlay, original, search) == {k: v for k, v in report.items() if k != 'source_sha256'}
    supported_ids = set(report['candidate_reaction_evidence_ids'])
    reactions = {r['id']: r for r in network['reactions']}
    selected = [r for r in reactions.values() if r['id'] in supported_ids]
    for source_scenario, scenario in zip(original['scenarios'], report['scenarios']):
        seeds = set(scenario['seed_compound_ids'])
        assert seeds == set(source_scenario['seed_compound_ids'])
        result = expand(selected, seeds)
        assert result['witnesses'] == scenario['witnesses']
        assert [t['cannabisdb_id'] for t in scenario['targets']] == [t['cannabisdb_id'] for t in source_scenario['targets']]
        for frontier in scenario['frontiers']:
            assert frontier['reaction_id'] not in supported_ids
            assert all(m['compound_id'] in result['available'] for m in frontier['required_inputs'])
            rescued = expand(selected + [reactions[frontier['reaction_id']]], seeds)
            assert rescued['witnesses'] == frontier['rescued_witnesses']
            difference = rescued['available'] - result['available']
            assert set(frontier['rescued_compound_ids']) == difference
            assert frontier['rescued_target_ids'] == [t['cannabisdb_id'] for t in scenario['targets'] if t['compound_id'] in difference]
            for c, w in frontier['rescued_witnesses'].items():
                if c not in seeds:
                    assert w['reaction_id'] in supported_ids | {frontier['reaction_id']}
                    assert all(frontier['rescued_witnesses'][m['compound_id']]['level'] < w['level'] for m in w['required_inputs'])
