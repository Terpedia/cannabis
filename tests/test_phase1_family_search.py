import pytest
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_family_search import parse_references, parse_hits, annotate_rows


def test_sequence_response_cannot_substitute_another_accession():
    assert parse_references(b'>sp|P1|NAME protein\nMALW\n', {'P1'})['P1']['sequence'] == 'MALW'
    for data in [b'>sp|P2|NAME protein\nMALW\n', b'<html>error</html>', b'>P1 protein\nMALW\n>P1 protein\nMALW\n']:
        with pytest.raises(ValueError):
            parse_references(data, {'P1'})


def test_screen_requires_both_coverages_and_evalue():
    lines = ['Q1\tP1\t35\t100\t200\t200\t1e-8\t100\t50\t50',
             'Q2\tP1\t90\t100\t200\t400\t1e-8\t100\t50\t25',
             'Q3\tP1\t90\t100\t200\t200\t1e-3\t100\t50\t50']
    hits = parse_hits('\n'.join(lines), {'Q1', 'Q2', 'Q3'}, {'P1'})
    assert [h['passes_screen'] for h in hits['P1']] == [True, False, False]
    with pytest.raises(ValueError, match='absent'):
        parse_hits(lines[0], {'Q2'}, {'P1'})
    rows = annotate_rows([
        {'reaction_id': 'RHEA:1', 'family_reference_annotations': [{'accession': 'P1'}]},
        {'reaction_id': 'RHEA:2', 'family_reference_annotations': [{'accession': 'P2'}]},
        {'reaction_id': 'RHEA:3', 'family_reference_annotations': []}], {'P1': {}}, hits)
    assert rows[0]['screened_cannabis_proteins'] == ['Q1']
    assert rows[0]['direction_status'].startswith('unresolved')
    assert rows[1]['reference_sequences_missing'] == ['P2']
    assert len(rows) == 3


def test_published_family_search_retains_all_gaps_and_screening_boundaries():
    root = Path(__file__).resolve().parents[1]
    discovery_path = root / 'data/reports/phase1-reference-discovery.json'
    discovery = json.loads(discovery_path.read_text())
    report = json.loads((root / 'docs/data/phase1-family-protein-search.json').read_text())
    assert report['source_discovery_sha256'] == hashlib.sha256(discovery_path.read_bytes()).hexdigest()
    key = lambda row: (row['reaction_id'], row['reaction_smarts'])
    assert [key(r) for r in report['rows']] == [key(r) for r in discovery['rows']]
    proteins = {p['accession']: p for p in report['cannabis_candidates']}
    for protein in proteins.values():
        assert protein['sequence_sha256'] == hashlib.sha256(protein['sequence'].encode()).hexdigest()
    for row in report['rows']:
        refs = {r['accession'] for r in row['family_reference_annotations']}
        passing = set()
        for hit in row['sequence_hits']:
            assert hit['reference_accession'] in refs
            expected = (hit['identity_percent'] >= 30 and hit['query_coverage_percent'] >= 50
                        and hit['reference_coverage_percent'] >= 50 and 0 <= hit['evalue'] <= 1e-5)
            assert hit['passes_screen'] == expected
            if expected:
                passing.add(hit['cannabis_accession'])
        assert set(row['screened_cannabis_proteins']) == passing
        assert passing <= proteins.keys()
        assert row['direction_status'] == 'unresolved-for-requested-reaction-direction'
