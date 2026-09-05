import copy
import hashlib
import json
from pathlib import Path
import pytest

from cannabis_carbon.phase1_deferred_references import queue
from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_deferred_discovery_requires_every_prior_attempt_to_be_explicitly_skipped():
    gap = {'reaction_id': 'R', 'left': [], 'right': [], 'source_joins': [], 'selected_uses': []}
    precursors = {'catalog_candidate_gaps': [gap], 'focused_targets': []}
    candidates = {'candidate_reaction_evidence_ids': {}}
    searches = {'prior': {'source_discovery': 'discovery', 'rows': [{'reaction_id': 'R',
        'search_status': 'no-reference-sequence', 'passing_alignment_ids': [], 'reference_sequences_present': []}]}}
    discoveries = {'discovery': {'rows': [{'reaction_id': 'R', 'lookup_status': 'not-searched-in-priority-pass', 'reference_matches': []}]}}
    assert len(queue(precursors, candidates, searches, discoveries, {})[0]) == 1
    for status in ('no-reviewed-reference-returned', 'lookup-incomplete-or-failed', None):
        discoveries['discovery']['rows'][0]['lookup_status'] = status
        rows, audit = queue(precursors, candidates, searches, discoveries, {})
        assert rows == [] and audit[0]['disposition'] == 'retain-other-prior-evidence'
    searches['prior']['rows'][0]['passing_alignment_ids'] = ['alignment']
    with pytest.raises(ValueError, match='Prior candidate omitted'):
        queue(precursors, candidates, searches, discoveries, {})


def test_complete_remaining_gap_audit_and_reference_attachment_replay(monkeypatch):
    monkeypatch.chdir(ROOT)
    report = read('phase1-deferred-references')
    searches = {p: json.loads((ROOT / p).read_text()) for p in report['search_report_paths']}
    discoveries = {s['source_discovery']: json.loads((ROOT / s['source_discovery']).read_text()) for s in searches.values()}
    source, candidates = read('phase1-purine-precursor-audit'), read('phase1-purine-candidate-net')
    inputs = (source, candidates, searches, discoveries)
    before = copy.deepcopy(inputs)
    families = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    rows, audit = queue(*inputs, families)
    assert audit == report['prior_attempt_audit']
    proteins = attach(rows, report['lookups'])
    assert rows == report['rows']
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    assert [proteins[a] for a in sorted(used)] == report['reference_proteins']
    assert inputs == before
    assert len(audit) == 60 and len(rows) == 7
    assert len(used) == 29 and sum(bool(r['reference_matches']) for r in rows) == 5
    assert all(r['reaction_id'] not in candidates['candidate_reaction_evidence_ids'] for r in rows)
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    for lookup in report['lookups']:
        assert lookup['status'] == 'retrieved'
        assert 'reviewed%3Atrue' in lookup['url']
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']


def test_candidate_leads_have_only_partial_reference_coverage():
    search = read('phase1-deferred-search')
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['retrieved_references'] == search['summary']['requested_references'] == 29
    assert search['summary']['distinct_cannabis_candidates'] == 9
    assert search['summary']['equations_with_screened_candidates'] == 1
    for candidate in search['cannabis_candidates']:
        hits = [h for h in search['passing_alignments'] if h['cannabis_accession'] == candidate['accession']]
        representative = max(hits, key=lambda h: h['bitscore'])
        assert 466 <= len(candidate['sequence']) <= 543
        assert representative['reference_length'] == 902
        assert 52.4 <= representative['reference_coverage_percent'] <= 54.0
