import hashlib
import json
from pathlib import Path


def test_carrier_network_replacement_and_boundaries():
    read = lambda n: json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
    report = read('carrier-network'); restored = read('carrier-reconstruction'); old = read('geranial-reduction-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    replaced = {j['model_reaction_id'] for r in restored['restored_equations'] for j in r['source_joins']}
    reactions = {r['id']: r for r in report['reactions']}
    assert len(reactions) == len(report['reactions']) == 18222 - 1316 + 1342
    assert set(report['replaced_projected_reaction_ids']) == replaced
    assert not replaced & reactions.keys()
    assert {r['id'] for r in report['quarantined_original_equations']} == replaced
    for r in restored['restored_equations']:
        for key in ('left', 'right', 'source_joins'):
            assert reactions[r['id']][key] == r[key]
    compounds = {c['id']: c for c in report['compounds']}
    carrier_ids = {c['id'] for c in compounds.values() if c.get('identity_status') == 'source-defined-macromolecular-carrier'}
    expected_carriers = {c['id'] for c in restored['carriers']} & {
        p['compound_id'] for r in restored['restored_equations'] for side in ('left', 'right') for p in r[side]}
    assert carrier_ids == expected_carriers
    assert len(carrier_ids) == 81
    assert not carrier_ids & set(report['external_exchange_compound_ids'])
    assert report['external_exchange_compound_ids'] == old['external_exchange_compound_ids']
    for cid in carrier_ids:
        assert compounds[cid]['carbon_count'] is None
        assert not compounds[cid]['external_uptake_allowed']
    assert len({s['id'] for s in report['directed_steps']}) == len(report['directed_steps'])
    forbidden = set(old['forbidden_step_ids'])
    for s in report['directed_steps']:
        assert s['id'] not in forbidden
        assert not set(s.get('inherited_allowed_step_ids', [])) & forbidden
        r = reactions[s['reaction_id']]
        assert (s['required_inputs'], s['outputs']) in ((r['left'], r['right']), (r['right'], r['left']))
        assert all(p['compound_id'] in compounds for side in ('required_inputs', 'outputs') for p in s[side])
    assert [(t['cannabisdb_id'], t['compound_id']) for t in report['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in old['targets']]
    assert len(report['targets']) == 6220
    assert report['summary']['qualified_co2_pathway_coverage'] is None
