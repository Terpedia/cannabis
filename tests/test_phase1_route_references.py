import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_route_references import queue
from cannabis_carbon.phase1_new_references import attach


def test_queue_preserves_full_equation_and_missing_family_without_id_arithmetic():
    reaction = {'id': 'r', 'left': [{'compound_id': 'a', 'coefficient': 2}],
        'right': [{'compound_id': 'b', 'coefficient': 1}], 'enzyme_evidence_ids': [],
        'sources': [{'source_reaction_id': 'RHEA:42'}]}
    gap = {'reaction_id': 'r', 'sources': reaction['sources'], 'target_ids': ['t'],
           'route_indices': [0], 'selected_route_target_count': 1}
    rows = queue({'enzyme_gap_queue': [gap]}, {'reactions': [reaction]}, {})
    assert rows[0]['left'] == reaction['left']
    assert rows[0]['rhea_families'] == {}
    assert rows[0]['hypothesis_ids'] == []
    assert rows[0]['route_indices'] == [0]
    reaction['enzyme_evidence_ids'] = ['evidence']
    with pytest.raises(ValueError, match='evidence mismatch'):
        queue({'enzyme_gap_queue': [gap]}, {'reactions': [reaction]}, {})


def test_published_route_reference_provenance_and_every_gap():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-route-references.json'
    report = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-route-references.json').read_bytes()
    cert = json.loads((root / 'data/reports/phase1-route-certificates.json').read_text())
    network = json.loads((root / 'data/reports/phase1-full-balanced-network.json').read_text())
    reactions = {r['id']: r for r in network['reactions']}
    assert {r['reaction_id'] for r in report['rows']} == {g['reaction_id'] for g in cert['enzyme_gap_queue']}
    proteins = {p['accession']: p for p in report['reference_proteins']}
    for row in report['rows']:
        reaction = reactions[row['reaction_id']]
        assert row['left'] == reaction['left'] and row['right'] == reaction['right']
        assert row['sources'] == reaction['sources']
        families = {rid for f in row['rhea_families'].values() for rid in f.values()}
        for m in row['reference_matches']:
            ids = set(proteins[m['accession']]['annotated_rhea_ids'])
            assert set(m['family_annotation_matches']) == ids & families
            assert set(m['exact_source_id_matches']) == ids & set(row['source_reaction_ids'])
            assert m['direction_status'] == 'not-established-for-hypothesis'
    for source, digest in report['source_sha256'].items():
        source_path = root / source
        if source_path.exists():
            assert hashlib.sha256(source_path.read_bytes()).hexdigest() == digest
    if all((root / l['snapshot']).exists() for l in report['lookups'] if l['status'] == 'retrieved'):
        rows = json.loads(json.dumps(report['rows']))
        recovered = attach(rows, [{**l, 'snapshot': str(root / l['snapshot'])} if l['status'] == 'retrieved' else l for l in report['lookups']])
        for row in rows:
            if not row['rhea_families']:
                row['lookup_status'] = 'no-published-Rhea-family-mapping'
        assert rows == report['rows']
        assert all(recovered[a]['annotated_rhea_ids'] == p['annotated_rhea_ids'] for a, p in proteins.items())
