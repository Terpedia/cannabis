import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_ureide_evidence_bundle import SOURCES
from cannabis_carbon.phase1_reference_gap_bundle import build


def test_ureide_bundle_preserves_all_sources_and_exact_export(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    path = Path('data/reports/phase1-ureide-evidence-bundle.json')
    report = json.loads(path.read_text())
    assert report == build(SOURCES)
    assert len(report['documents']) == 8
    assert Path('docs/data/ureide-evidence-bundle.json').read_bytes() == path.read_bytes()
    rows = [json.loads(line) for line in Path('data/derived/phase1-ureide-evidence-bundle.ndjson').read_text().splitlines()]
    assert len(rows) == 9
    assert {r['report_sha256'] for r in rows} == {hashlib.sha256(path.read_bytes()).hexdigest()}
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'source_document'] == report['documents']
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'metadata'] == [{k: v for k, v in report.items() if k != 'documents'}]
