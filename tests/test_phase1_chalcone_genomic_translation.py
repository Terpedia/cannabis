import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_chalcone_genomic_translation import reconstruct, CODE


def test_reverse_exon_join_crosses_codon_boundary_and_rejects_bad_bases():
    # Reverse-complement exon order yields AT + GAAA = ATGAAA.
    exons = [{'start': 7, 'end': 8, 'strand': -1}, {'start': 1, 'end': 4, 'strand': -1}]
    cds, codons = reconstruct('TTTCCCAT', 1, exons)
    assert cds == 'ATGAAA'
    assert [c['amino_acid'] for c in codons] == ['M', 'K']
    assert codons[0]['genomic_positions'] == [8, 7, 4]
    with pytest.raises(ValueError, match='Ambiguous'):
        reconstruct('NTTCCCAT', 1, exons)
    with pytest.raises(ValueError, match='outside'):
        reconstruct('TTTCCCA', 1, exons)
    assert len(CODE) == 64
    assert {c for c, a in CODE.items() if a == '*'} == {'TAA', 'TAG', 'TGA'}


def test_full_chalcone_genomic_translation_and_all_coordinates(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-genomic-translation.json').read_text())
    model = json.loads(Path('data/reports/phase1-chalcone-gene-model.json').read_text())['source_record']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    fasta = next(Path(p) for p in report['source_sha256'] if p.endswith('.fasta'))
    sequence = ''.join(fasta.read_text().splitlines()[1:])
    cds, codons = reconstruct(sequence, report['genomic_region']['start'], model['transcript_order_exons'])
    assert report['spliced_cds'] == cds
    assert report['codons'] == codons
    assert len(cds) == 1041 and len(sequence) == 2132
    assert report['translation_including_stop'] == model['sequence'] + '*'
    assert ''.join(c['amino_acid'] for c in codons) == model['sequence'] + '*'
    expected_coords = [p for e in model['transcript_order_exons'] for p in range(e['end'], e['start']-1, -1)]
    assert [p for c in codons for p in c['genomic_positions']] == expected_coords
    assert len(set(expected_coords)) == 1041
    assert codons[-1]['codon'] == 'TGA'
    assert report['model_eligible'] is False
