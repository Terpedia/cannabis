import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_reviewed_weak_search import build
from cannabis_carbon.phase1_screened_overlay import build_overlay


def test_review_rejection_is_enforced_without_erasing_original_search(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-reviewed-weak-search.json').read_text())
    assert report == build()
    original = json.loads(Path('data/reports/phase1-weak-nonplant-search.json').read_text())
    assert report['rows'] == original['rows']
    assert report['summary'] == original['summary']
    assert report['passing_alignments'][0]['passes_screen'] is True
    assert report['passing_alignments'][0]['model_eligible'] is False
    assert 'model_eligible' not in original['passing_alignments'][0]
    parent = {'reactions': [{'id': r['reaction_id']} for r in report['rows']],
              'hypotheses': [{'id': h, 'reaction_id': r['reaction_id']} for r in report['rows'] for h in r['hypothesis_ids']]}
    with pytest.raises(ValueError, match='Alignment evidence is explicitly ineligible'):
        build_overlay(parent, report)
