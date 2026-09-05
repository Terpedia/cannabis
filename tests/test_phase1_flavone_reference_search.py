import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_flavone_screen_preserves_exact_oxygenase_gap_and_is_not_promoted(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-flavone-search')
    search = json.loads(Path('data/reports/phase1-flavone-search.json').read_text())
    discovery = json.loads(Path('data/reports/phase1-flavone-references.json').read_text())
    parent = json.loads(Path('data/reports/phase1-chalcone-remaining-gaps.json').read_text())
    row = discovery['rows'][0]
    gap = parent['candidate_gaps'][0]
    assert row['reaction_id'] == gap['reaction_id']
    for side in ('left', 'right'):
        assert row[side] == gap['reaction'][side]
    assert row['sources'] == gap['reaction']['sources']
    assert set(row['target_ids']) == {'CDB005071', 'CDB005072'}
    assert {u['direction_mode'] for u in row['selected_uses']} == {'hypothetical-right-to-left'}
    for path, sha in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['requested_references'] == search['summary']['retrieved_references'] == 1
    assert search['summary']['raw_alignments'] == 105
    assert search['summary']['passing_alignments'] == 32
    assert search['model_eligible'] is discovery['model_eligible'] is False
    assert all(r['model_eligible'] is False for r in search['rows'])
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay({'reactions': [{'id': row['reaction_id']}], 'hypotheses': []}, search)
