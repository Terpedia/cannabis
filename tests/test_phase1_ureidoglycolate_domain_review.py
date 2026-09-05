import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.genome import _fasta


def test_weak_hit_coordinates_and_domain_annotation_are_sequence_linked(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-ureidoglycolate-domain-review.json').read_text())
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    search = json.loads(Path('data/reports/phase1-ureidoglycolate-broad-search.json').read_text())
    sequences = _fasta(Path(search['proteome_path']))
    refs = {r['accession']: r['sequence'] for r in search['reference_sequences']}
    for acc, snapshot in report['annotation_snapshots'].items():
        assert snapshot['sequence']['value'] == (sequences[acc] if acc == report['query_accession'] else refs[acc])
    coords = Path(report['coordinate_command'][report['coordinate_command'].index('--out') + 1])
    lines = [line.split('\t') for line in coords.read_text().splitlines()]
    assert len(lines) == report['coordinate_replay_alignment_count'] == 13
    assert Counter('\t'.join(f[:10]) for f in lines) == Counter(Path(search['hits_path']).read_text().splitlines())
    assert report['selected_alignment_columns'] in lines
    assert report['selected_alignment_columns'][10:] == ['50', '409', '4', '352']
    assert len(report['domain_overlaps']) == 2
    for row in report['domain_overlaps']:
        domain = row['domain']
        assert domain in report['annotation_snapshots'][row['accession']]['features']
        assert domain['description'] == 'Isopropylmalate dehydrogenase-like'
        overlap = min(row['alignment_end'], domain['location']['end']['value']) - max(row['alignment_start'], domain['location']['start']['value']) + 1
        assert row['overlapping_residues'] == max(0, overlap)
        assert overlap / (row['alignment_end'] - row['alignment_start'] + 1) > 0.95
    assert report['model_eligible'] is False
    assert search['summary']['passing_alignments'] == 0
