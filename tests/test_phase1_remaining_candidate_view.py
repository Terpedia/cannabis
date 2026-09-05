import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_remaining_candidate_view import build

ROOT = Path(__file__).resolve().parents[1]


def test_remaining_view_replays_all_evidence_sources_and_exact_target_inventory():
    folder = ROOT / 'docs/data/remaining-net-view'
    manifest = json.loads((folder / 'index.json').read_text())
    payload = (folder / 'bundle.json').read_bytes()
    assert len(payload) == manifest['bytes']
    assert hashlib.sha256(payload).hexdigest() == manifest['sha256']
    reports = []
    for path, sha in manifest['source_sha256'].items():
        data = (ROOT / path).read_bytes()
        assert hashlib.sha256(data).hexdigest() == sha
        reports.append(json.loads(data))
    bundle = json.loads(payload)
    assert build(reports[0], reports[1:]) == bundle
    for scenario, original in [(bundle, reports[0]['scenarios'][0]), (bundle['restricted_scenario'], reports[0]['scenarios'][1])]:
        assert len(scenario['targets']) == 6220
        assert scenario['certificates'] == original['certificates']
        assert scenario['forbidden_step_ids'] == original['forbidden_step_ids']
    evidence = {e['id']: e for e in bundle['enzyme_evidence']}
    for e in reports[0]['enzyme_evidence']:
        assert evidence[e['id']] == e
    assert bundle['probe_results'] == reports[0]['probe_results']
