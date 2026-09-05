import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_dopa_site_review import map_sites
from cannabis_carbon.phase1_flavone_site_review import verify_numbering, restore_explicit_masks


def test_numbering_and_mask_validation_fail_closed():
    assert verify_numbering({'tested_substitutions': ['M1T']}, 'M', 'T') == {1: {'fht_residue': 'M', 'fnsi_residue': 'T'}}
    for substitutions in [['M1T', 'M1T'], ['M2T'], ['M1S'], ['invalid']]:
        with pytest.raises(ValueError):
            verify_numbering({'tested_substitutions': substitutions}, 'M', 'T')
    fields = [''] * 10 + ['1', '3', '1', '3', 'AX-C', 'A-BC']
    restored, masks = restore_explicit_masks(fields, 'ABC', 'ABC')
    assert restored[14:] == ['AB-C', 'A-BC']
    assert masks == [{'side': 'query', 'position': 2, 'reported_residue': 'X', 'pinned_residue': 'B'}]
    assert fields[14] == 'AX-C'
    with pytest.raises(ValueError, match='Non-mask'):
        restore_explicit_masks(fields[:14] + ['AZ-C', 'A-BC'], 'ABC', 'ABC')
    with pytest.raises(ValueError, match='shorter'):
        restore_explicit_masks(fields[:14] + ['AB--', 'A-BC'], 'ABC', 'ABC')


def test_full_coordinate_replay_and_every_lead_site(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-flavone-site-review.json').read_text())
    comparison = json.loads(Path('data/reports/phase1-flavone-fht-comparison.json').read_text())
    original = json.loads(Path('data/reports/phase1-flavone-search.json').read_text())
    refs = {r['accession']: r['annotation']['sequence']['value'] for r in comparison['references']}
    fnsi = next(r['sequence'] for r in original['reference_sequences'] if r['accession'] == 'Q7XZQ8')
    positions = verify_numbering(report['review'], refs['Q7XZQ7'], fnsi)
    assert set(positions) == {106, 115, 116, 131, 195, 215, 216}
    assert report['review'] == json.loads(Path('data/curation/flavone-seven-site-review.json').read_text())
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    lines = [l.split('\t') for l in Path('data/raw/flavone-site-review/coordinate-hits.tsv').read_text().splitlines()]
    assert len(lines) == report['coordinate_replay_alignment_count'] == 218
    assert Counter('\t'.join(f[:10]) for f in lines) == Counter(Path('data/raw/flavone-fht-comparison/hits.tsv').read_text().splitlines())
    queries = _fasta(Path(original['proteome_path']))
    masked = []
    for fields in lines:
        acc, ref = fields[0].split('|')[1], fields[1]
        restored, masks = restore_explicit_masks(fields, queries[acc], refs[ref])
        map_sites(restored, queries[acc], refs[ref], [])
        if masks:
            masked.append({'accession': acc, 'reference_accession': ref, 'masks': masks})
    assert masked == report['masked_alignments']
    assert len(masked) == 7
    assert {r['accession'] for r in report['rows']} == {r['accession'] for r in comparison['rows']}
    status = Counter()
    for row in report['rows']:
        assert row['model_eligible'] is False
        mapped = row['mapped_alignment']
        fields = mapped['alignment_columns']
        assert fields in lines
        restored, masks = restore_explicit_masks(fields, queries[row['accession']], refs['Q7XZQ7'])
        assert mapped['reported_masks'] == masks
        sites = map_sites(restored, queries[row['accession']], refs['Q7XZQ7'], positions)
        for actual, expected in zip(mapped['sites'], sites, strict=True):
            assert {k: actual[k] for k in expected} == expected
            assert actual['fht_residue'] == positions[expected['reference_position']]['fht_residue']
            assert actual['fnsi_residue'] == positions[expected['reference_position']]['fnsi_residue']
            assert actual['reported_masking'] == [m for m in masks if
                (m['side'] == 'query' and m['position'] == actual['query_position']) or
                (m['side'] == 'reference' and m['position'] == actual['reference_position'])]
            status[actual['status']] += 1
    assert status == {'aligned-residue': 215, 'outside-local-alignment': 6, 'query-gap': 3}
    assert report['model_eligible'] is False
