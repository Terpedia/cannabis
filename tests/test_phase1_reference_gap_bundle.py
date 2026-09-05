import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_reference_gap_bundle import SOURCES, build


def test_gap_export_is_lossless_and_does_not_promote_evidence(monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.chdir(root)
    path = Path('data/reports/phase1-reference-gap-bundle.json')
    report = json.loads(path.read_text())
    assert Path('docs/data/reference-gap-bundle.json').read_bytes() == path.read_bytes()
    assert report == build()
    assert {d['source_path'] for d in report['documents']} == set(SOURCES)
    for item in report['documents']:
        assert item['document'] == json.loads(Path(item['source_path']).read_text())
        assert item['source_sha256'] == hashlib.sha256(Path(item['source_path']).read_bytes()).hexdigest()
    rows = [json.loads(line) for line in Path('data/derived/phase1-reference-gap-bundle.ndjson').read_text().splitlines()]
    assert len(rows) == len(SOURCES) + 1
    assert len({(r['record_kind'], r['record_id']) for r in rows}) == len(rows)
    assert {r['report_sha256'] for r in rows} == {hashlib.sha256(path.read_bytes()).hexdigest()}
    assert [json.loads(r['record_json']) for r in rows if r['record_kind'] == 'source_document'] == report['documents']
    assert report['summary']['candidate_model_changed'] is False
