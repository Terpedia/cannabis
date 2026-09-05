import copy
import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_reference_discovery import direction_families
from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_weighted_gap_search import queue
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search

ROOT = Path(__file__).resolve().parents[1]


def test_weighted_gap_full_proteome_search_replays():
    verify_search('phase1-weighted-gap-search')
    report = json.loads((ROOT / 'data/reports/phase1-weighted-gap-search.json').read_text())
    for name in ('proteome', 'reference', 'hits'):
        assert (ROOT / report[name + '_path']).is_file()
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['retrieved_references'] == 75
    assert report['summary']['passing_alignments'] == 355
    assert report['summary']['distinct_cannabis_candidates'] == 5
    for row in report['rows']:
        assert row['evidence_class'] == 'direction-unresolved-reference-homology-candidate'
        assert 'physiological-direction-unverified' in row['validation_blockers']


def test_discovery_replays_exact_family_and_excludes_previously_searched_gaps(monkeypatch):
    monkeypatch.chdir(ROOT)
    source = json.loads((ROOT / 'data/reports/phase1-evidence-weighted-routes.json').read_text())
    report = json.loads((ROOT / 'data/reports/phase1-weighted-gap-references.json').read_text())
    families = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    rows = queue(source, families)
    assert len(rows) == 1
    assert set(rows[0]['source_reaction_ids']) == {'RHEA:21037', 'RHEA:21038'}
    assert {f['RHEA_ID_MASTER'] for f in rows[0]['rhea_families'].values()} == {'RHEA:21036'}
    proteins = attach(rows, report['lookups'])
    assert rows == report['rows']
    # Discovery iterates sets of Rhea IDs; accession, not iteration order, is identity.
    assert len(proteins) == len(report['reference_proteins'])
    assert proteins == {p['accession']: p for p in report['reference_proteins']}
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    for lookup in report['lookups']:
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']
    changed = copy.deepcopy(source)
    for gap in changed['candidate_gaps']:
        gap['prior_searches'] = [{'search_status': 'no-hits'}]
    assert queue(changed, families) == []
