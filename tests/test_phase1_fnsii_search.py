import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_generic_fnsii_screen_never_becomes_exact_fnsi_edge(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-fnsii-search')
    report = json.loads(Path('data/reports/phase1-fnsii-search.json').read_text())
    discovery = json.loads(Path('data/reports/phase1-fnsii-references.json').read_text())
    audit = json.loads(Path('data/reports/phase1-fnsii-alternative-audit.json').read_text())
    row = report['rows'][0]
    assert row['reaction_id'] == 'RHEA:57681' != row['related_fnsi_gap_id']
    assert row['related_fnsi_gap_id'] == audit['parent_fnsi_gap']['reaction_id']
    assert row['generic_source_record'] == next(r for r in audit['source_records'] if r['record']['rule_id'] == 'RHEA:57681')
    assert row['carrier_review'] == audit['review']
    assert set(row['target_ids']) == {'CDB005071', 'CDB005072'}
    assert report['model_eligible'] is discovery['model_eligible'] is row['model_eligible'] is False
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    assert all(r['model_eligible'] is False and not r['exact_reaction_annotation_match'] for r in row['reference_matches'])
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['retrieved_references'] == 3
    assert report['summary']['raw_alignments'] == 647
    assert report['summary']['passing_alignments'] == 258
    assert report['summary']['distinct_cannabis_candidates'] == 95
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay({'reactions': [{'id': row['reaction_id']}], 'hypotheses': []}, report)
