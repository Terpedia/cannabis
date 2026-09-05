import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_reference_discovery import direction_families
from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_weighted_gap_search import queue
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search

ROOT = Path(__file__).resolve().parents[1]


def test_remaining_gap_whole_proteome_and_exact_discovery(monkeypatch):
    monkeypatch.chdir(ROOT)
    verify_search('phase1-remaining-gap-search')
    report = json.loads(Path('data/reports/phase1-remaining-gap-search.json').read_text())
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['equation_gaps'] == 10
    for name in ('proteome', 'reference', 'hits'):
        assert Path(report[name + '_path']).is_file()
    discovery = json.loads(Path(report['source_discovery']).read_text())
    source = json.loads(Path('data/reports/phase1-remaining-weighted-routes.json').read_text())
    rows = queue(source, direction_families(Path('data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text()))
    proteins = attach(rows, discovery['lookups'])
    assert rows == discovery['rows']
    assert len(proteins) == len(discovery['reference_proteins'])
    assert proteins == {p['accession']: p for p in discovery['reference_proteins']}
    expected = {g['reaction_id']: g for g in source['candidate_gaps'] if not g['prior_searches']}
    assert {r['reaction_id'] for r in rows} == expected.keys()
    for row in rows:
        assert row['target_ids'] == expected[row['reaction_id']]['target_ids']
        assert row['selected_uses'] == expected[row['reaction_id']]['selected_uses']
    for path, sha in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
