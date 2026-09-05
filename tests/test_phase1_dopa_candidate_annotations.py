import hashlib
import json
from pathlib import Path


def test_all_dopa_leads_have_sequence_linked_gene_domain_context(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-dopa-candidate-annotations.json').read_text())
    search = json.loads(Path('data/reports/phase1-dopa-lyase-search.json').read_text())
    sequences = {p['accession']: p['sequence'] for p in search['cannabis_candidates']}
    assert len(report['rows']) == len(sequences) == 7
    assert {r['accession'] for r in report['rows']} == sequences.keys()
    assert report['model_eligible'] is False
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    for row in report['rows']:
        snapshot = json.loads(Path(row['snapshot']).read_text())
        assert row['annotation'] == snapshot
        assert snapshot['sequence']['value'] == sequences[row['accession']]
        assert set(row['passing_alignment_ids']) == {a['id'] for a in search['passing_alignments'] if a['cannabis_accession'] == row['accession']}
        assert row['model_eligible'] is False
