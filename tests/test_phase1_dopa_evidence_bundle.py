import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_dopa_evidence_bundle import SOURCES
from cannabis_carbon.phase1_reference_gap_bundle import build


def test_dopa_bundle_preserves_all_leads_and_evidence_boundaries(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    path = Path('data/reports/phase1-dopa-evidence-bundle.json')
    report = json.loads(path.read_text())
    assert report == build(SOURCES)
    assert len(report['documents']) == 5
    assert report['summary']['new_exact_enzyme_assignments'] == 0
    assert report['summary']['candidate_model_changed'] is False
    assert all(d['document']['model_eligible'] is False for d in report['documents'])
    assert Path('docs/data/dopa-evidence-bundle.json').read_bytes() == path.read_bytes()
    rows = [json.loads(line) for line in Path('data/derived/phase1-dopa-evidence-bundle.ndjson').read_text().splitlines()]
    assert len(rows) == 6
    assert {r['report_sha256'] for r in rows} == {hashlib.sha256(path.read_bytes()).hexdigest()}
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'source_document'] == report['documents']
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'metadata'] == [{k: v for k, v in report.items() if k != 'documents'}]
