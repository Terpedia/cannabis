import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_net_view import build


def test_missing_or_conflicting_evidence_fails_closed():
    report = {'reactions': [{'enzyme_evidence_ids': ['e']}], 'targets': [], 'compounds': []}
    with pytest.raises(ValueError, match='Missing'):
        build(report, [])
    with pytest.raises(ValueError, match='Conflicting'):
        build(report, [{'enzyme_evidence': [{'id': 'e', 'v': 1}]}, {'enzyme_evidence': [{'id': 'e', 'v': 2}]}])


def test_published_net_bundle_preserves_every_certificate_and_evidence_object():
    root = Path(__file__).resolve().parents[1]
    folder = root / 'docs/data/net-view'
    manifest = json.loads((folder / 'index.json').read_text())
    raw = (folder / manifest['file']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == manifest['sha256']
    assert len(raw) == manifest['bytes']
    sources = []
    for source, digest in manifest['source_sha256'].items():
        source_raw = (root / source).read_bytes()
        assert hashlib.sha256(source_raw).hexdigest() == digest
        sources.append(json.loads(source_raw))
    bundle = json.loads(raw)
    assert bundle == build(sources[0], sources[1:])
    assert bundle['certificates'] == sources[0]['certificates']
    assert bundle['targets'] == sources[0]['targets']
    assert bundle['reactions'] == sources[0]['reactions']
    assert len(raw) < 8_000_000
