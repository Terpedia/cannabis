import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_family_search import parse_hits
from cannabis_carbon.phase1_flavone_fht_comparison import compare


def test_comparison_retains_weak_missing_and_rejects_duplicate_hits():
    source = {'cannabis_candidates': [{'accession': 'q'}], 'passing_alignments': []}
    hit = {'cannabis_accession': 'q', 'reference_accession': 'r', 'passes_screen': False}
    row = compare(source, {'r': [hit]}, ['r', 's'])[0]
    assert row['model_eligible'] is False
    assert [c['status'] for c in row['comparators']] == ['weak', 'no-reported-alignment']
    assert row['comparators'][1]['alignment'] is None
    with pytest.raises(ValueError, match='Duplicate'):
        compare(source, {'r': [hit, hit]}, ['r'])


def test_published_comparison_preserves_full_screen_and_all_original_leads(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-flavone-fht-comparison.json').read_text())
    source = json.loads(Path('data/reports/phase1-flavone-search.json').read_text())
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    refs = {r['accession']: r for r in report['references']}
    assert set(refs) == {'Q7XZQ6', 'Q7XZQ7'}
    fasta = _fasta(Path('data/raw/flavone-fht-comparison/references.fasta'))
    for acc, ref in refs.items():
        assert ref['annotation'] == json.loads(Path(ref['snapshot']).read_text())
        assert ref['annotation']['primaryAccession'] == acc
        assert ref['annotation']['sequence']['value'] == fasta[acc]
    queries = _fasta(Path(source['proteome_path']))
    hits = parse_hits(Path('data/raw/flavone-fht-comparison/hits.tsv').read_text(), queries, refs)
    assert report['all_comparator_alignments'] == hits
    assert report['rows'] == compare(source, hits, refs)
    assert report['model_eligible'] is False
    assert report['screen'] == source['screen']
    assert len(queries) == report['summary']['proteome_sequences'] == 30304
    assert sum(map(len, hits.values())) == 218
    assert len(report['rows']) == len(source['cannabis_candidates']) == 32
    assert report['summary']['comparator_status_counts'] == {'passing': 58, 'weak': 6}
    for group in hits.values():
        for hit in group:
            assert hit['query_length'] == len(queries[hit['cannabis_accession']])
            assert hit['reference_length'] == len(fasta[hit['reference_accession']])
    review = json.loads(Path('data/curation/flavone-fht-comparison-review.json').read_text())
    for row in report['rows']:
        if row['accession'] in review['leading_candidates']:
            fht = next(c['alignment'] for c in row['comparators'] if c['reference_accession'] == 'Q7XZQ7')
            fnsi = row['original_fnsi_alignments'][0]
            assert fht['identity_percent'] > fnsi['identity_percent']
            assert fht['bitscore'] > fnsi['bitscore']
    assert review['model_eligible'] is False
    assert review['primary_study']['cannabis_activity_established'] is False
