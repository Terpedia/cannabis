import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_missing_reference_review import prepare
from cannabis_carbon.phase1_new_references import attach

ROOT = Path(__file__).resolve().parents[1]


def test_missing_reference_review_preserves_exact_gaps_and_negative_query(monkeypatch):
    monkeypatch.chdir(ROOT)
    report = json.loads(Path('data/reports/phase1-missing-reference-review.json').read_text())
    discovery = json.loads(Path('data/reports/phase1-remaining-gap-references.json').read_text())
    search = json.loads(Path('data/reports/phase1-remaining-gap-search.json').read_text())
    rows = prepare(discovery, search)
    assert len(rows) == 2
    assert {r['reaction_id'] for r in rows} == {r['reaction_id'] for r in report['rows']}
    lookup = report['lookups'][0]
    assert 'taxonomy_id:33090' in lookup['query']
    assert 'reviewed:false' in lookup['query']
    assert lookup['requested_master_ids'] == ['RHEA:30767', 'RHEA:32247']
    assert attach(rows, report['lookups']) == {}
    for before, after in zip(rows, report['rows']):
        assert before['prior_reviewed_search'] == after['prior_reviewed_search']
        assert before['left'] == after['left'] and before['right'] == after['right']
        assert before['reference_matches'] == after['reference_matches'] == []
        assert after['lookup_status'] == 'no-unreviewed-plant-reference-returned'
    assert report['summary']['unreviewed_plant_reference_records'] == 0
    assert hashlib.sha256(Path(lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
