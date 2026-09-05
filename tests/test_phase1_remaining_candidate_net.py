import hashlib
import json
from pathlib import Path
import test_phase1_thiolase_candidate_net as previous_tests
from cannabis_carbon.phase1_marts_completions import balanced

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_remaining_model_retains_all_records_and_only_adds_screened_equations():
    report = read('phase1-remaining-candidate-net')
    parent = read('phase1-thiolase-candidate-net')
    search = read('phase1-remaining-gap-search')
    network = read('phase1-full-balanced-network')
    expected = {r['reaction_id'] for r in search['rows'] if r['search_status'] == 'screened-candidates'}
    assert len(expected) == 4
    ids = set(report['candidate_reaction_evidence_ids'])
    assert len(ids) == 1609
    assert ids - parent['candidate_reaction_evidence_ids'].keys() == expected
    for rid, evidence in parent['candidate_reaction_evidence_ids'].items():
        assert report['candidate_reaction_evidence_ids'][rid] == evidence
    assert report['constraints'] == parent['constraints']
    assert report['synthase_reference_links'] == parent['synthase_reference_links']
    assert report['probes'] == parent['probes']
    assert {(r['scenario_id'], r['compound_id']) for r in report['probe_results']} == {(r['scenario_id'], r['compound_id']) for r in parent['probe_results']}
    compounds = {c['id']: c for c in network['compounds']}
    assert all(balanced([r['left'], r['right']], compounds) for r in network['reactions'] if r['id'] in ids)
    for scenario, old in zip(report['scenarios'], parent['scenarios']):
        assert len(scenario['targets']) == 6220
        assert [(t['cannabisdb_id'], t['compound_id']) for t in scenario['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in old['targets']]
        assert {c['compound_id'] for c in old['certificates']} <= {c['compound_id'] for c in scenario['certificates']}
        assert scenario['forbidden_step_ids'] == old['forbidden_step_ids']
    searches = {r['reaction_id']: r for r in search['rows']}
    assert {e['reaction_id'] for e in report['enzyme_evidence']} == expected
    for evidence in report['enzyme_evidence']:
        row = searches[evidence['reaction_id']]
        assert evidence['evidence_class'] == row['evidence_class']
        assert {p['accession'] for p in evidence['screened_proteins']} == set(row['screened_cannabis_proteins'])
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_every_remaining_model_target_and_probe_certificate_has_exact_co2_balance(monkeypatch):
    # Replay the same exhaustive certificate verifier against the new report.
    monkeypatch.setattr(previous_tests, 'read', lambda name: read('phase1-remaining-candidate-net' if name == 'phase1-thiolase-candidate-net' else name))
    previous_tests.test_all_thiolase_target_and_probe_certificates_balance_exactly()
