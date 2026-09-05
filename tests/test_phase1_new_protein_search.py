import hashlib
import json
from collections import defaultdict
from pathlib import Path
import pytest
from cannabis_carbon.phase1_new_protein_search import annotate
from cannabis_carbon.phase1_family_search import parse_hits
from cannabis_carbon.genome import _fasta


def test_retains_missing_and_weak_hits_without_enzyme_confirmation():
    rows = [{'reaction_id': str(i), 'target_ids': [str(i)], 'priority_target_ids': [str(i)],
             'hypothesis_ids': [str(i)], 'reference_matches': [{'accession': ref}]} for i, ref in enumerate(['R1', 'R2', 'R3', 'R4'])]
    hits = parse_hits('Q1\tR1\t30\t50\t100\t100\t1e-5\t50\t50\t50\nQ2\tR2\t90\t40\t100\t100\t1e-10\t50\t40\t40', {'Q1', 'Q2'}, {'R1', 'R2'})
    annotated, passing = annotate({'rows': rows}, {'R1': {}, 'R2': {}, 'R3': {}}, hits)
    assert [r['search_status'] for r in annotated] == ['screened-candidates', 'weak-hits-only', 'no-hits', 'no-reference-sequence']
    assert annotated[0]['screened_cannabis_proteins'] == ['Q1']
    assert annotated[3]['reference_sequences_missing'] == ['R4']
    assert len(passing) == 1
    assert all('physiological-direction-unverified' in r['validation_blockers'] for r in annotated)


@pytest.mark.parametrize('report_name', ['phase1-new-protein-search', 'phase1-route-protein-search', 'phase1-completion-protein-search', 'phase1-archived-protein-search', 'phase1-catalog-protein-search', 'phase1-backfill-protein-search', 'phase1-plant-purine-search', 'phase1-purine-gap-search', 'phase1-deferred-search', 'phase1-replacement-search'])
def test_published_search_preserves_gap_scope_thresholds_sequences_and_joins(report_name):
    root = Path(__file__).resolve().parents[1]
    path = root / f'data/reports/{report_name}.json'
    report = json.loads(path.read_text())
    discovery = root / report['source_discovery']
    assert hashlib.sha256(discovery.read_bytes()).hexdigest() == report['source_discovery_sha256']
    source_rows = json.loads(discovery.read_text())['rows']
    assert [r['reaction_id'] for r in source_rows] == [r['reaction_id'] for r in report['rows']]
    refs = {r['accession']: r for r in report['reference_sequences']}
    candidates = {r['accession']: r for r in report['cannabis_candidates']}
    for sequence in [*refs.values(), *candidates.values()]:
        assert hashlib.sha256(sequence['sequence'].encode()).hexdigest() == sequence['sequence_sha256']
    alignments = {a['id']: a for a in report['passing_alignments']}
    by_reference = defaultdict(set)
    assert len(alignments) == len(report['passing_alignments'])
    for alignment in alignments.values():
        by_reference[alignment['reference_accession']].add(alignment['id'])
        assert alignment['identity_percent'] >= 30
        assert alignment['query_coverage_percent'] >= 50
        assert alignment['reference_coverage_percent'] >= 50
        assert 0 <= alignment['evalue'] <= 1e-5
        assert alignment['query_length'] == len(candidates[alignment['cannabis_accession']]['sequence'])
        assert alignment['reference_length'] == len(refs[alignment['reference_accession']]['sequence'])
    for before, row in zip(source_rows, report['rows']):
        assert row['reference_matches'] == before['reference_matches']
        allowed = {m['accession'] for m in row['reference_matches']}
        assert set(row['reference_sequences_present']) == allowed & refs.keys()
        expected = {aid for ref in allowed for aid in by_reference[ref]}
        assert set(row['passing_alignment_ids']) == expected
        assert set(row['screened_cannabis_proteins']) == {alignments[aid]['cannabis_accession'] for aid in expected}
    for name in ['proteome', 'reference', 'hits']:
        raw = root / report[name + '_path']
        if raw.exists():
            assert hashlib.sha256(raw.read_bytes()).hexdigest() == report[name + '_sha256']
    if (root / report['hits_path']).exists() and (root / report['proteome_path']).exists():
        queries = _fasta(root / report['proteome_path'])
        assert len(queries) == report['summary']['proteome_sequences']
        raw_hits = parse_hits((root / report['hits_path']).read_text(), queries.keys(), refs.keys())
        assert sum(map(len, raw_hits.values())) == report['summary']['raw_alignments']
        recovered = {hashlib.sha256(json.dumps(h, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                     for group in raw_hits.values() for h in group if h['passes_screen']}
        assert recovered == alignments.keys()
        replayed, _ = annotate({'rows': source_rows}, refs, raw_hits)
        for actual, expected in zip(report['rows'], replayed):
            for key in ('search_status', 'raw_alignment_count', 'passing_alignment_ids', 'screened_cannabis_proteins'):
                assert actual[key] == expected[key]
