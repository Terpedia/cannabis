import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_route_overlay import combine, annotate_routes
from cannabis_carbon.phase1_screened_overlay import apply_overlay, build_overlay


def test_rejects_colliding_screen_evidence():
    overlay = {'enzyme_evidence': [{'id': 'same'}]}
    with pytest.raises(ValueError, match='collision'):
        combine({'reactions': []}, overlay, overlay)


def test_published_route_evidence_and_map_preserve_all_chemistry():
    root = Path(__file__).resolve().parents[1]
    def read(name):
        return json.loads((root / 'data/reports' / (name + '.json')).read_text())
    parent = read('phase1-target-hypotheses')
    previous = read('phase1-screened-enzyme-overlay')
    network = read('phase1-full-balanced-network')
    search = read('phase1-route-protein-search')
    cert = read('phase1-route-certificates')
    added = read('phase1-route-enzyme-overlay')
    combined = read('phase1-combined-enzyme-overlay')
    status = read('phase1-route-evidence-status')
    for name in ['phase1-route-enzyme-overlay', 'phase1-combined-enzyme-overlay', 'phase1-route-evidence-status']:
        path = root / 'data/reports' / (name + '.json')
        assert path.read_bytes() == (root / 'docs/data' / path.name).read_bytes()
        for source, digest in read(name)['source_sha256'].items():
            assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    rebuilt = build_overlay({'reactions': network['reactions'], 'hypotheses': []}, search, 'phase1-route-protein-search.json')
    assert rebuilt['enzyme_evidence'] == added['enzyme_evidence']
    expected = {r['reaction_id'] for r in search['rows'] if r['screened_cannabis_proteins']}
    assert {e['reaction_id'] for e in added['enzyme_evidence']} == expected
    assert all(e['full_search_report'] == 'phase1-route-protein-search.json' for e in added['enzyme_evidence'])
    rebuilt_combined = combine(parent, previous, added)
    assert rebuilt_combined['enzyme_evidence'] == combined['enzyme_evidence']
    assert rebuilt_combined['summary'] == combined['summary']
    integrated = apply_overlay(parent, combined)
    assert integrated['summary'] == combined['integrated_summary']
    assert integrated['compounds'] == parent['compounds']
    assert integrated['reactions'] == parent['reactions']
    assert integrated['targets'] == parent['targets']
    for old, new in zip(parent['hypotheses'], integrated['hypotheses']):
        for key in ('id', 'reaction_id', 'cannabisdb_id', 'compound_id', 'required_inputs', 'outputs', 'net_target_coefficient', 'status', 'direction_mode'):
            assert old[key] == new[key]
        assert set(old['blockers']) - set(new['blockers']) <= {'no-candidate-enzyme-evidence-attached'}
    assert annotate_routes(cert, added) == {k: v for k, v in status.items() if k != 'source_sha256'}
    added_by_reaction = {e['reaction_id']: e for e in added['enzyme_evidence']}
    for old, new in zip(cert['routes'], status['routes']):
        assert (new['cannabisdb_id'], new['scenario_id'], new['compound_id']) == (old['cannabisdb_id'], old['scenario_id'], old['compound_id'])
        assert new['status'] == old['status'] and new['biological_status'] == old['biological_status']
        for old_step, step in zip(old['steps'], new['steps']):
            extra = added_by_reaction.get(old_step['reaction_id'])
            assert step['step_id'] == old_step['id']
            assert step['enzyme_evidence_ids'] == old_step['enzyme_evidence_ids'] + ([extra['id']] if extra else [])
            assert step['has_candidate_enzyme_evidence'] == bool(step['enzyme_evidence_ids'])
        assert new['missing_enzyme_step_ids'] == [s['step_id'] for s in new['steps'] if not s['has_candidate_enzyme_evidence']]
