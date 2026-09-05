import copy
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_new_references import attach

ROOT = Path(__file__).resolve().parents[1]


def test_nonplant_negative_reference_query_is_replayable_and_not_activity_evidence(monkeypatch):
    monkeypatch.chdir(ROOT)
    report = json.loads(Path('data/reports/phase1-nonplant-reference-review.json').read_text())
    parent = json.loads(Path('data/reports/phase1-missing-reference-review.json').read_text())
    rows = copy.deepcopy(parent['rows'])
    lookup = report['lookups'][0]
    assert 'NOT taxonomy_id:33090' in lookup['query']
    assert 'reviewed:false' in lookup['query'] and 'fragment:false' in lookup['query']
    assert lookup['requested_master_ids'] == ['RHEA:30767', 'RHEA:32247']
    assert attach(rows, report['lookups']) == {}
    assert len(report['rows']) == len(rows) == 2
    for before, after in zip(rows, report['rows']):
        for key in ('reaction_id', 'left', 'right', 'rhea_families', 'target_ids', 'prior_reviewed_search'):
            assert before[key] == after[key]
        assert after['lookup_status'] == 'no-unreviewed-nonplant-reference-returned'
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    review = json.loads(Path('data/curation/glucuronolactone-literature-review.json').read_text())
    assert review['reaction_id'] in {r['reaction_id'] for r in rows}
    assert review['candidate_model_changed'] is False
    assert all(s['read_scope'] == 'publisher abstract' for s in review['sources'])
