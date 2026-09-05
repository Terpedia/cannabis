import json
from pathlib import Path
from cannabis_carbon.phase1_allantoate_alternative_gaps import build


def test_all_alternative_route_gaps_have_pinned_search_history(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-allantoate-alternative-gaps.json').read_text())
    search_paths = [p for p in report['source_sha256'] if p.endswith('search.json')]
    assert report == build(search_paths)
    routes = json.loads(Path('data/reports/phase1-allantoate-sensitivity.json').read_text())
    model = json.loads(Path('data/reports/phase1-remaining-candidate-net.json').read_text())
    expected = {s['reaction_id'] for r in routes['rows'] for s in r['steps']} - model['candidate_reaction_evidence_ids'].keys()
    assert {g['reaction_id'] for g in report['candidate_gaps']} == expected
    assert len(expected) == 10
    assert report['summary']['previously_unsearched_equations'] == 0
    for gap in report['candidate_gaps']:
        assert set(gap['target_ids']) == {r['cannabisdb_id'] for r in routes['rows'] if any(s['reaction_id'] == gap['reaction_id'] for s in r['steps'])}
        assert gap['prior_searches']
        for prior in gap['prior_searches']:
            actual = json.loads(Path(prior['report']).read_text())
            assert prior['row'] in actual['rows']
