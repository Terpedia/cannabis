import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_dopa_literature_leads_preserve_exact_scope_and_cannot_promote(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-dopa-lyase-search')
    search = json.loads(Path('data/reports/phase1-dopa-lyase-search.json').read_text())
    discovery = json.loads(Path(search['source_discovery']).read_text())
    review = json.loads(Path('data/curation/dopa-lyase-reference-review.json').read_text())
    priority = json.loads(Path('data/reports/phase1-current-gap-priority.json').read_text())
    gap = next(r for r in priority['rows'] if r['reaction_id'] == review['reaction_id'])
    for path, sha in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    row = discovery['rows'][0]
    assert row['target_ids'] == gap['remaining_target_ids']
    assert len(row['target_ids']) == 12
    assert {t for use in row['selected_uses'] for t in use['target_ids']} == set(row['target_ids'])
    assert row['historical_selected_uses'] == gap['selected_uses']
    for side in ('left', 'right'):
        assert row[side] == gap['reaction'][side]
    match = row['reference_matches'][0]
    assert match['accession'] == 'Q3IWB0'
    assert match['exact_reaction_annotation_match'] is False
    assert match['model_eligible'] is False
    assert discovery['model_eligible'] is search['model_eligible'] is False
    ref = json.loads(Path('data/raw/phase1-dopa-lyase-search/Q3IWB0.json').read_text())
    assert search['reference_sequences'][0]['sequence'] == ref['sequence']['value']
    assert len(ref['sequence']['value']) == 523
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['raw_alignments'] == 9
    assert search['summary']['passing_alignments'] == 7
    for name in ('proteome', 'reference', 'hits'):
        assert Path(search[name + '_path']).is_file()
    parent = {'reactions': [{'id': row['reaction_id']}], 'hypotheses': []}
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay(parent, search)
