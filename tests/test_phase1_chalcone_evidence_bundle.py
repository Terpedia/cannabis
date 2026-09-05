import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_chalcone_evidence_bundle import SOURCES
from cannabis_carbon.phase1_reference_gap_bundle import build


def test_chalcone_export_preserves_complete_source_documents(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    path = Path('data/reports/phase1-chalcone-evidence-bundle.json')
    report = json.loads(path.read_text())
    assert report == build(SOURCES)
    assert len(report['documents']) == 5
    assert all(d['document']['model_eligible'] is False for d in report['documents'])
    assert report['summary']['candidate_model_changed'] is False
    assert Path('docs/data/chalcone-evidence-bundle.json').read_bytes() == path.read_bytes()
    rows = [json.loads(line) for line in Path('data/derived/phase1-chalcone-evidence-bundle.ndjson').read_text().splitlines()]
    assert len(rows) == 6
    assert {r['report_sha256'] for r in rows} == {hashlib.sha256(path.read_bytes()).hexdigest()}
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'source_document'] == report['documents']
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'metadata'] == [{k: v for k, v in report.items() if k != 'documents'}]
