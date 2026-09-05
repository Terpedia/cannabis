import csv
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_ureidoglycolate_broad_search import build, QUERIES
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_broad_reference_discovery_retains_both_queries_and_all_annotations(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-ureidoglycolate-broad-references.json').read_text())
    assert report == build(report['lookups'])
    assert {item['query'] for item in report['lookups']} == set(QUERIES.values())
    expected = {}
    for item in report['lookups']:
        with Path(item['snapshot']).open() as stream:
            records = list(csv.DictReader(stream, delimiter='\t'))
        for record in records:
            expected.setdefault(record['Entry'], []).append({'lookup_url': item['url'], 'query': item['query'], 'record': record})
    row = report['rows'][0]
    assert row['prior_exact_reference_matches'] == []
    assert row['model_eligible'] is report['model_eligible'] is False
    assert len(expected) == report['summary']['distinct_reference_leads'] == 2124
    assert {m['accession']: m['source_records'] for m in row['reference_matches']} == expected
    assert all(m['model_eligible'] is False and m['exact_reaction_annotation_match'] is False for m in row['reference_matches'])


def test_full_proteome_broad_search_is_not_exact_reaction_evidence(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-ureidoglycolate-broad-search')
    report = json.loads(Path('data/reports/phase1-ureidoglycolate-broad-search.json').read_text())
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['requested_references'] == report['summary']['retrieved_references'] == 2124
    assert report['summary']['failed_retrieval_batches'] == 0
    for name in ('proteome', 'reference', 'hits'):
        assert Path(report[name + '_path']).is_file()
    assert report['model_eligible'] is False
    assert all(r['model_eligible'] is False for r in report['rows'])
    parent = {'reactions': [{'id': r['reaction_id']} for r in report['rows']],
              'hypotheses': [{'id': h, 'reaction_id': r['reaction_id']} for r in report['rows'] for h in r['hypothesis_ids']]}
    if report['passing_alignments']:
        with pytest.raises(ValueError, match='ineligible'):
            build_overlay(parent, report)
    else:
        assert build_overlay(parent, report)['enzyme_evidence'] == []
