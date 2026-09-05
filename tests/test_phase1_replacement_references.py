import copy
import hashlib
import json
from pathlib import Path

import pytest

from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_reference_discovery import direction_families
from cannabis_carbon.phase1_replacement_references import queue

ROOT = Path(__file__).resolve().parents[1]


def test_new_gap_discovery_preserves_prior_negatives_and_withheld_candidates():
    source = {'candidate_gaps': [{'reaction_id': 'R', 'selected_uses': [], 'source_joins': []}],
              'reactions': [{'id': 'R', 'left': [], 'right': []}], 'focused_targets': []}
    assert len(queue(source, {}, {})[0]) == 1
    for status, hits in [('no-hits', []), ('no-reference-sequence', []), ('weak-hits-only', []), ('screened-candidates', ['H'])]:
        prior = {'prior': {'rows': [{'reaction_id': 'R', 'search_status': status, 'passing_alignment_ids': hits}]}}
        rows, audit = queue(source, prior, {})
        assert rows == []
        assert audit[0]['prior_screens'][0]['row'] == prior['prior']['rows'][0]
        assert audit[0]['disposition'] == 'retain-prior-evidence-without-promotion'
    source['candidate_gaps'] *= 2
    with pytest.raises(ValueError, match='Duplicate gap'):
        queue(source, {}, {})


def test_complete_inventory_and_reviewed_reference_replay(monkeypatch):
    monkeypatch.chdir(ROOT)
    report = json.loads((ROOT / 'data/reports/phase1-replacement-references.json').read_text())
    source = json.loads((ROOT / 'data/reports/phase1-decay-sensitivity.json').read_text())
    searches = {p: json.loads((ROOT / p).read_text()) for p in report['search_report_paths']}
    families = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    before = copy.deepcopy((source, searches))
    rows, audit = queue(source, searches, families)
    proteins = attach(rows, report['lookups'])
    assert rows == report['rows']
    assert audit == report['prior_screen_audit']
    assert len(rows) == 5 and len(audit) == 66
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    assert [proteins[a] for a in sorted(used)] == report['reference_proteins']
    assert before == (source, searches)
    assert all(not a['prior_screens'] for a in audit if a['disposition'] == 'new-reference-discovery')
    assert sum(any(p['row'].get('passing_alignment_ids') for p in a['prior_screens']) for a in audit) == 1
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    for lookup in report['lookups']:
        assert lookup['status'] == 'retrieved'
        assert 'reviewed%3Atrue' in lookup['url']
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']


def test_new_candidate_architecture_coverage_does_not_resolve_reverse_direction():
    discovery = json.loads((ROOT / 'data/reports/phase1-replacement-references.json').read_text())
    search = json.loads((ROOT / 'data/reports/phase1-replacement-search.json').read_text())
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['requested_references'] == search['summary']['retrieved_references'] == 524
    assert search['summary']['equations_with_screened_candidates'] == 2
    assert search['summary']['distinct_cannabis_candidates'] == 4
    citrate = next(r for r in discovery['rows'] if 'RHEA:16846' in r['source_reaction_ids'])
    assert next(s for s in citrate['sources'] if s['source_reaction_id'] == 'RHEA:16846')['source_left_corresponds_to'] == 'right'
    assert {u['direction_mode'] for u in citrate['selected_uses']} == {'hypothetical-left-to-right'}
    for candidate in search['cannabis_candidates']:
        hit = max((h for h in search['passing_alignments'] if h['cannabis_accession'] == candidate['accession']), key=lambda h: h['bitscore'])
        assert hit['reference_coverage_percent'] >= 93.3
        assert hit['query_coverage_percent'] >= 81
        assert hit['identity_percent'] >= 65
    model = json.loads((ROOT / 'data/reports/phase1-purine-candidate-net.json').read_text())
    assert all(r['reaction_id'] not in model['candidate_reaction_evidence_ids'] for r in search['rows'])
