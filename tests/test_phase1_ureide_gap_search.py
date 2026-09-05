import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_reference_discovery import direction_families
from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_weighted_gap_search import queue
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_ureide_gap_search_preserves_all_proteins_and_exact_reaction_joins(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-ureide-gap-search')
    report = json.loads(Path('data/reports/phase1-ureide-gap-search.json').read_text())
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['equation_gaps'] == 3
    assert report['summary']['retrieved_references'] == 3
    assert report['summary']['raw_alignments'] == 1
    assert report['summary']['passing_alignments'] == 0
    for name in ('proteome', 'reference', 'hits'):
        assert Path(report[name + '_path']).is_file()
    discovery = json.loads(Path(report['source_discovery']).read_text())
    source = json.loads(Path('data/reports/phase1-ureide-alternative-gaps.json').read_text())
    rows = queue(source, direction_families(Path('data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text()))
    proteins = attach(rows, discovery['lookups'])
    assert rows == discovery['rows']
    assert proteins == {p['accession']: p for p in discovery['reference_proteins']}
    assert {r['reaction_id'] for r in rows} == {g['reaction_id'] for g in source['candidate_gaps'] if not g['prior_searches']}
    for path, sha in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
