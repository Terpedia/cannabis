import hashlib
import json
from pathlib import Path


def test_cpr_annotations_preserve_full_source_evidence(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-cpr-annotations.json').read_text())
    search = json.loads(Path('data/reports/phase1-cpr-search.json').read_text())
    sequences = {p['accession']: p['sequence'] for p in search['reference_sequences'] + search['cannabis_candidates']}
    assert len(report['rows']) == len(sequences) == 14
    assert {r['accession'] for r in report['rows']} == set(sequences)
    assert report['model_eligible'] is False
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    for row in report['rows']:
        data = json.loads(Path(row['snapshot']).read_text())
        assert row['annotation'] == data
        assert data['primaryAccession'] == row['accession']
        assert data['sequence']['value'] == sequences[row['accession']]
        assert row['model_eligible'] is False
        assert set(row['passing_alignment_ids']) == {a['id'] for a in search['passing_alignments'] if row['accession'] in (a['cannabis_accession'], a['reference_accession'])}
