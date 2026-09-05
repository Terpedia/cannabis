import hashlib
import json
from pathlib import Path


def test_every_chalcone_reference_and_lead_annotation_matches_searched_sequence(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-annotations.json').read_text())
    search = json.loads(Path('data/reports/phase1-chalcone-search.json').read_text())
    sequences = {p['accession']: p['sequence'] for p in search['reference_sequences'] + search['cannabis_candidates']}
    refs = {p['accession'] for p in search['reference_sequences']}
    assert len(report['rows']) == len(sequences) == 55
    assert {r['accession'] for r in report['rows']} == set(sequences)
    assert report['model_eligible'] is False
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    for row in report['rows']:
        acc = row['accession']
        data = json.loads(Path(row['snapshot']).read_text())
        assert data == row['annotation']
        assert data['primaryAccession'] == acc
        assert data['sequence']['value'] == sequences[acc]
        assert row['role'] == ('reference' if acc in refs else 'Cannabis-lead')
        assert row['model_eligible'] is False
        assert set(row['passing_alignment_ids']) == {a['id'] for a in search['passing_alignments'] if acc in (a['reference_accession'], a['cannabis_accession'])}
