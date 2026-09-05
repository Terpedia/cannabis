import copy
import hashlib
import json
from pathlib import Path

import pytest

from cannabis_carbon.phase1_gap_annotations import assemble, queue
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_gap_annotation_replay_and_evidence_boundary(monkeypatch):
    monkeypatch.chdir(ROOT)
    catalog = read('phase1-catalog-net-gaps')
    evidence = read('phase1-combined-catalog-evidence')
    report = read('phase1-gap-annotations')
    before = copy.deepcopy((catalog, evidence))
    families = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    rows = queue(catalog, evidence, families)
    assert assemble(rows, report['lookups']) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert (catalog, evidence) == before
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    assert len(rows) == len({r['id'] for r in rows}) == 349
    added = {e['reaction_id'] for e in evidence['enzyme_evidence']}
    assert {r['id'] for r in rows} == {g['id'] for g in catalog['gap_priorities']} - added
    for row in report['rows']:
        assert row['spontaneous_status'].startswith('unresolved')
        assert (row['annotation_status'] == 'source-ec-linked-Cannabis-function-unresolved') == bool(row['ec_annotation_urls'])
    top = report['rows'][0]
    assert top['annotation_master_ids'] == ['RHEA:46952']
    assert 'http://purl.uniprot.org/enzyme/1.17.1.11' in top['ec_annotation_urls']
    assert 'http://rdf.ncbi.nlm.nih.gov/pubmed/23872566' in top['citation_urls']
    assert top['selected_certificate_target_count'] == 181
    payload = (ROOT / 'data/reports/phase1-gap-annotations.json').read_bytes()
    exported = [json.loads(line) for line in (ROOT / 'data/derived/phase1-gap-annotations.ndjson').read_text().splitlines()]
    assert all(r['report_sha256'] == hashlib.sha256(payload).hexdigest() for r in exported)
    assert len(exported) == len({(r['record_kind'], r['record_id']) for r in exported})


def test_no_ec_or_missing_family_never_means_spontaneous(tmp_path):
    payload = json.dumps({'results': {'bindings': [{'s': {'type': 'uri', 'value': 'http://rdf.rhea-db.org/10'},
        'p': {'type': 'uri', 'value': 'http://rdf.rhea-db.org/equation'},
        'o': {'type': 'literal', 'value': 'A = B'}}]}}).encode()
    snapshot = tmp_path / 'source.json'
    snapshot.write_bytes(payload)
    lookup = {'snapshot': str(snapshot), 'sha256': hashlib.sha256(payload).hexdigest(), 'requested_master_ids': ['RHEA:10']}
    rows = [{'id': 'mapped', 'source_joins': [{'master_id': 'RHEA:10', 'source': {'source_reaction_id': 'RHEA:11'}}]},
            {'id': 'unmapped', 'source_joins': [{'master_id': None, 'source': {'source_reaction_id': 'unknown'}}]},
            {'id': 'missing', 'source_joins': [{'master_id': 'RHEA:20', 'source': {'source_reaction_id': 'RHEA:21'}}]}]
    output = assemble(rows, [lookup])
    assert [r['annotation_status'] for r in output['rows']] == ['no-source-ec-link-catalysis-unresolved', 'source-annotation-incomplete', 'source-annotation-incomplete']
    assert output['summary']['new_spontaneous_claims'] == output['summary']['new_Cannabis_activity_claims'] == 0
    snapshot.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='checksum'):
        assemble(rows, [lookup])
