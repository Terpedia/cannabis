import hashlib
import json
from pathlib import Path

import pytest

from cannabis_carbon.phase1_catalog_references import queue
from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def test_prior_unsuccessful_screens_are_preserved_not_resubmitted():
    reaction = {'id': 'r', 'enzyme_evidence_ids': [], 'missing_candidate_evidence': True,
        'left': [], 'right': [], 'sources': [{'source_reaction_id': 'RHEA:10000'}]}
    gap = {'id': 'r', 'reaction_id': 'r', 'selected_certificate_target_ids': ['target'],
        'selected_certificate_target_count': 1}
    catalog = {'reactions': [reaction], 'gap_priorities': [gap]}
    rows, audit = queue(catalog, {}, {})
    assert len(rows) == 1 and rows[0]['target_ids'] == ['target']
    assert rows[0]['rhea_families'] == {}  # never manufacture numeric families
    prior = {'old': {'rows': [{'reaction_id': 'r', 'search_status': 'weak-hits-only',
        'reference_sequences_present': ['ref'], 'passing_alignment_ids': []}]}}
    rows, audit = queue(catalog, prior, {})
    assert rows == []
    assert audit[0]['prior_screens'][0]['search_status'] == 'weak-hits-only'
    prior['old']['rows'][0]['passing_alignment_ids'] = ['passing']
    with pytest.raises(ValueError, match='Prior passing'):
        queue(catalog, prior, {})
    reaction['enzyme_evidence_ids'] = ['candidate']
    with pytest.raises(ValueError, match='already has'):
        queue(catalog, {}, {})


def test_catalog_reference_snapshot_replays_exact_sources_and_prior_partition(monkeypatch):
    monkeypatch.chdir(ROOT)
    read = lambda path: json.loads(Path(path).read_text())
    report = read('data/reports/phase1-catalog-references.json')
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    catalog = read('data/reports/phase1-catalog-net-gaps.json')
    prior = {p: read(p) for p in report['source_sha256'] if p.endswith('protein-search.json')}
    families = direction_families(Path(report['rhea_direction_source']['snapshot']).read_text())
    rows, audit = queue(catalog, prior, families)
    assert audit == report['prior_screen_audit']
    assert len(audit) == 465 and len(rows) == 199
    assert len({r['reaction_id'] for r in rows}) == len(rows)
    proteins = attach(rows, report['lookups'])
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    assert rows == report['rows']
    used = {m['accession'] for row in rows for m in row['reference_matches']}
    assert [proteins[a] for a in sorted(used)] == report['reference_proteins']
    for row in rows:
        family_ids = {rid for family in row['rhea_families'].values() for rid in family.values()}
        for match in row['reference_matches']:
            annotated = set(proteins[match['accession']]['annotated_rhea_ids'])
            assert match['family_annotation_matches'] == sorted(annotated & family_ids)
            assert match['family_annotation_matches']
            assert match['exact_source_id_matches'] == sorted(annotated & set(row['source_reaction_ids']))


def test_new_candidates_cover_only_one_additional_selected_certificate():
    read = lambda name: json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())
    catalog, screen = read('phase1-catalog-net-gaps'), read('phase1-catalog-protein-search')
    candidates = {r['reaction_id'] for r in screen['rows'] if r['screened_cannabis_proteins']}
    assert len(candidates) == 97
    assert candidates <= {g['reaction_id'] for g in catalog['gap_priorities']}
    closed = [t['cannabisdb_id'] for t in catalog['targets'] if t['missing_candidate_reaction_ids']
        and set(t['missing_candidate_reaction_ids']) <= candidates]
    assert closed == ['CDB004839']
    assert sum(bool(set(t['missing_candidate_reaction_ids']) & candidates) for t in catalog['targets']) == 181
    # Evidence comparison only: frozen catalog still exposes all original gaps.
    assert catalog['summary']['selected_missing_candidate_equations'] == 465
