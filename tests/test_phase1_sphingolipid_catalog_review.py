import json
from pathlib import Path
from cannabis_carbon.phase1_sphingolipid_catalog_review import build


def test_alternative_references_join_catalog_exclusions_without_promotion(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-sphingolipid-catalog-review.json').read_text())
    assert report == build()
    assert {r['rhea_id'] for r in report['rows']} == {'RHEA:46268', 'RHEA:46272'}
    excluded = {e['source_reaction_id'] for r in report['rows'] for e in r['excluded_catalog_records']}
    assert excluded == {'RHEA:46269', 'RHEA:46270', 'RHEA:46273', 'RHEA:46274'}
    for row in report['rows']:
        assert row['model_eligible'] is False
        assert not row['accepted_balanced_reaction_ids']
        assert row['passing_alignments']
        assert all(e['status'] == 'not-auditable' for e in row['excluded_catalog_records'])
