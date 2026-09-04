import json
from cannabis_carbon.phase1_core_evidence import attach


def test_explicit_rhea_directional_join_does_not_cross_namespaces(tmp_path):
    search, network, output = [tmp_path / n for n in ('search.json', 'network.json', 'out.json')]
    search.write_text(json.dumps({'rows': [{'reaction_id': x, 'search_status': 'no-reference-sequence'} for x in ['RHEA:10025', 'MARTS:10025', 'RHEA:10028']]}))
    evidence = [{'enzyme_id': 'p1', 'qualifiers': {'directExperimentalEvidence': False}, 'sources': [{'url': 'https://example.org/evidence'}]}]
    network.write_text(json.dumps({'reactions': [{'id': 'rhea:10024', 'directional_rhea_ids': ['10025'], 'enzyme_ids': ['p1'], 'enzyme_associations': evidence}]}))
    summary = attach(search, network, output)
    assert summary['missing_reference_variants_with_core_evidence'] == 1
    rows = json.loads(output.read_text())['rows']
    assert rows[0]['core_reaction_evidence'][0]['core_reaction_id'] == 'rhea:10024'
    assert not rows[1]['core_reaction_evidence']
    assert not rows[2]['core_reaction_evidence']
    assert rows[0]['core_reaction_evidence'][0]['enzyme_associations'] == evidence
    assert all(row['search_status'] == 'no-reference-sequence' for row in rows)
