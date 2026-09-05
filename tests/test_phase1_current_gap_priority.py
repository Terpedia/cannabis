import json
from pathlib import Path
from cannabis_carbon.phase1_current_gap_priority import build


def test_priority_retains_all_gaps_and_excludes_current_route_targets(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-current-gap-priority.json').read_text())
    assert report == build()
    diagnostic = json.loads(Path('data/reports/phase1-remaining-weighted-routes.json').read_text())
    assert {r['reaction_id'] for r in report['rows']} == {r['reaction_id'] for r in diagnostic['candidate_gaps']}
    for row in report['rows']:
        assert row['remaining_target_count'] == len(row['remaining_target_ids'])
        assert all(s not in ('exact-net-conversion-hypothesis', 'explicit-exchange-species; not a synthesis target') for s in row['remaining_target_statuses'].values())
        assert {p['id'] for p in row['participants']} == {p['compound_id'] for side in ('left', 'right') for p in row['reaction'][side]}
