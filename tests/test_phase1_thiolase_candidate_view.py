import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_thiolase_candidate_view import build

ROOT = Path(__file__).resolve().parents[1]


def test_view_replays_sources_and_retains_every_synthase_protein_link():
    folder = ROOT / 'docs/data/thiolase-net-view'
    manifest = json.loads((folder / 'index.json').read_text())
    payload = (folder / 'bundle.json').read_bytes()
    assert hashlib.sha256(payload).hexdigest() == manifest['sha256']
    assert len(payload) == manifest['bytes']
    reports = []
    for path, sha in manifest['source_sha256'].items():
        data = (ROOT / path).read_bytes()
        assert hashlib.sha256(data).hexdigest() == sha
        reports.append(json.loads(data))
    bundle = json.loads(payload)
    assert build(reports[0], reports[1:]) == bundle
    evidence = {e['id']: e for e in bundle['enzyme_evidence']}
    for link in reports[0]['synthase_reference_links']:
        row = evidence['synthase-reference-link:' + link['id']]
        assert row['source_link'] == link
        assert [p['accession'] for p in row['screened_proteins']] == [p['candidate_accession'] for p in link['protein_links']]
    assert bundle['probe_results'] == reports[0]['probe_results']
