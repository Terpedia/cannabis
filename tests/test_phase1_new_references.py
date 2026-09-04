import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_new_references import queue, attach


def test_queue_retains_all_carbon_gaps_and_prioritizes_targets_without_alternatives():
    family = {'RHEA_ID_MASTER': 'RHEA:1', 'RHEA_ID_LR': 'RHEA:2', 'RHEA_ID_RL': 'RHEA:3', 'RHEA_ID_BI': 'RHEA:4'}
    report = {'compounds': [{'id': 'carbon', 'carbon_count': 1}, {'id': 'water', 'carbon_count': 0}],
        'reactions': [{'id': 'r', 'left': [], 'right': [], 'sources': [{'source_reaction_id': 'RHEA:2'}]}],
        'hypotheses': [{'id': str(i), 'reaction_id': 'r', 'compound_id': c, 'cannabisdb_id': t, 'has_candidate_enzyme_evidence': support}
                      for i, (c, t, support) in enumerate([('carbon', 'A', False), ('carbon', 'B', False), ('carbon', 'B', True), ('water', 'W', False)])]}
    rows = queue(report, {'RHEA:2': family})
    assert len(rows) == 1
    assert rows[0]['target_ids'] == ['A', 'B']
    assert rows[0]['priority_target_ids'] == ['A']
    assert rows[0]['hypothesis_ids'] == ['0', '1']


def test_reference_family_join_keeps_direction_unresolved_and_failures_separate(tmp_path):
    raw = b'Entry\tRhea ID\tOrganism\tProtein names\nP1\tRHEA:1\tOther plant\tEnzyme\n'
    path = tmp_path / 'references.tsv'; path.write_bytes(raw)
    lookups = [{'status': 'retrieved', 'requested_master_ids': ['RHEA:1', 'RHEA:5'], 'snapshot': str(path), 'sha256': hashlib.sha256(raw).hexdigest(), 'url': 'https://rest.uniprot.org/query'},
               {'status': 'retrieval-failed', 'requested_master_ids': ['RHEA:9']}]
    rows = [{'source_reaction_ids': [f'RHEA:{n+1}'], 'rhea_families': {str(n): {'RHEA_ID_MASTER': f'RHEA:{n}', 'RHEA_ID_LR': f'RHEA:{n+1}'}}} for n in [1, 5, 9, 13]]
    proteins = attach(rows, lookups)
    assert proteins['P1']['retrieval_evidence'][0]['sha256'] == lookups[0]['sha256']
    assert rows[0]['reference_matches'][0]['exact_source_id_matches'] == []
    assert rows[0]['reference_matches'][0]['direction_status'] == 'not-established-for-hypothesis'
    assert [r['lookup_status'] for r in rows] == ['references-found', 'no-reviewed-reference-returned', 'lookup-incomplete-or-failed', 'not-searched-in-priority-pass']


def test_published_discovery_matches_every_gap_and_explicit_family():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / 'data/reports/phase1-new-references.json').read_text())
    parent = json.loads((root / 'data/reports/phase1-target-hypotheses.json').read_text())
    carbon = {c['id'] for c in parent['compounds'] if c['carbon_count']}
    expected = {h['id'] for h in parent['hypotheses'] if h['compound_id'] in carbon and not h['has_candidate_enzyme_evidence']}
    assert {hid for r in report['rows'] for hid in r['hypothesis_ids']} == expected
    proteins = {p['accession']: p for p in report['reference_proteins']}
    if all((root / lookup['snapshot']).exists() for lookup in report['lookups'] if lookup['status'] == 'retrieved'):
        replay = json.loads(json.dumps(report['rows']))
        lookups = [{**lookup, 'snapshot': str(root / lookup['snapshot'])} if lookup['status'] == 'retrieved' else lookup
                   for lookup in report['lookups']]
        recovered = attach(replay, lookups)
        assert {acc: p['annotated_rhea_ids'] for acc, p in recovered.items()} == {acc: p['annotated_rhea_ids'] for acc, p in proteins.items()}
        assert replay == report['rows']
    for source, digest in report['source_sha256'].items():
        path = root / source
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    for row in report['rows']:
        family_ids = {rid for f in row['rhea_families'].values() for rid in f.values()}
        for match in row['reference_matches']:
            ids = set(proteins[match['accession']]['annotated_rhea_ids'])
            assert set(match['family_annotation_matches']) == ids & family_ids
            assert set(match['exact_source_id_matches']) == ids & set(row['source_reaction_ids'])
            assert match['direction_status'] == 'not-established-for-hypothesis'
