"""Verify the public comparison against exact identities and the paired calculation."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_scope import orientations


def read(path):
    return json.loads(Path(path).read_bytes())


def test_paired_view_preserves_every_identity_result_and_denominator():
    folder = Path('docs/data/pg-named-net-view')
    payload = (folder / 'bundle.json').read_bytes()
    view = json.loads(payload)
    manifest = read(folder / 'index.json')
    assert manifest['file'] == 'bundle.json'
    assert manifest['bytes'] == len(payload)
    assert manifest['sha256'] == hashlib.sha256(payload).hexdigest()
    assert manifest['source_sha256'] == view['source_sha256']
    for path, sha in view['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    net = read('data/reports/phase1-pg-named-net.json')
    historical = read('data/reports/phase1-glycerophospholipid-net.json')
    assert view['historical_inventory_summary_unchanged'] == historical['summary']
    assert view['claim_boundary'] == view['view_boundary'] == net['claim_boundary']
    assert view['view_scenario'] == 'pg-paired-identity'
    pairs = {p['cannabisdb_id']: p for p in net['paired_probes']}
    assert len(pairs) == 23
    keys = ['original_result', 'alternative_baseline_result', 'alternative_extended_result']
    assert [s['id'] for s in view['scenario_options']] == keys
    for scenario in view['scenario_options']:
        key = scenario['id']
        assert len(scenario['targets']) == scenario['summary']['target_records'] == 23
        assert {t['cannabisdb_id'] for t in scenario['targets']} == set(pairs)
        expected_certs = [p[key] for p in net['paired_probes']
                          if p[key]['status'] == 'exact-net-conversion-hypothesis']
        assert scenario['certificates'] == expected_certs
        assert scenario['summary']['target_status_counts'] == dict(Counter(
            t['net_status'] for t in scenario['targets']))
        for target in scenario['targets']:
            pair = pairs[target['cannabisdb_id']]
            result = pair[key]
            assert target['compound_id'] == result['compound_id']
            assert target['net_status'] == result['status']
            assert target['source_url'] == pair['source_url']
            assert target['occurrence_assertion'] == pair['occurrence_assertion']
            assert scenario['title'] in target['label']
            assert target['certificate_compound_id'] == (
                result['compound_id'] if result in expected_certs else None)
            if key == 'original_result':
                assert target['identity_status'] == 'original-encoded-structure; name-conflict-unresolved'
                assert target['compound_id'] != pair['alternative_extended_result']['compound_id']
            else:
                assert target['identity_status'] == pair['identity_status']
    assert [len(s['certificates']) for s in view['scenario_options']] == [0, 0, 23]
    for field in ('targets', 'certificates', 'summary'):
        assert view[field] == view['scenario_options'][-1][field]


def test_public_equations_replay_certificates_and_preserve_hypothesis_provenance():
    view = read('docs/data/pg-named-net-view/bundle.json')
    net = read('data/reports/phase1-pg-named-net.json')
    proposals = read('data/reports/phase1-pg-named-reactions.json')
    source_urls = {r['rule_id']: r['source_url'] for r in proposals['source_records']}
    originals = {r['id']: r for r in net['certificate_reactions']}
    reactions = {r['id']: r for r in view['reactions']}
    assert len(reactions) == len(view['reactions']) == len(originals)
    assert reactions.keys() == originals.keys()
    compounds = {c['id']: c for c in view['compounds']}
    assert len(compounds) == len(view['compounds'])
    canonical = {c['id']: c for c in net['compounds']}
    for cid, compound in compounds.items():
        for key, value in canonical[cid].items():
            if key != 'labels':
                assert compound[key] == value
    for key in ('external_exchange_compound_ids', 'co2_compound_id', 'forbidden_step_ids'):
        assert view[key] == net[key]
    exchange = set(view['external_exchange_compound_ids'])
    assert {cid for cid in exchange if compounds[cid]['carbon_count']} == {view['co2_compound_id']}
    forbidden = set(view['forbidden_step_ids'])
    steps = {s['id']: s for s in orientations(view['reactions']) if s['id'] not in forbidden}
    evidence = {e['id'] for e in view['enzyme_evidence']}
    for rid, reaction in reactions.items():
        for key, value in originals[rid].items():
            if key not in ('sources', 'enzyme_evidence_ids'):
                assert reaction[key] == value
        assert balanced([reaction['left'], reaction['right']], compounds)
        assert reaction['sources'], rid
        assert set(reaction['enzyme_evidence_ids']) <= evidence
        if reaction.get('hypothesis_type'):
            assert reaction['is_route_sensitivity'] is True
            assert reaction['claim_boundary'] in reaction['hypothesis_assumptions']
            if reaction.get('source_reaction_id') in source_urls:
                assert source_urls[reaction['source_reaction_id']] in {
                    url for source in reaction['sources'] for url in source.get('source_urls', [])}
    for scenario in view['scenario_options']:
        for cert in scenario['certificates']:
            validate_certificate(cert, steps, compounds, exchange, view['co2_compound_id'])
            assert not forbidden & {s['step_id'] for s in cert['steps']}
        for target in scenario['targets']:
            cert = pairs_result(net, target['cannabisdb_id'], scenario['id'])
            assert target['missing_candidate_reaction_ids'] == sorted({
                s['reaction_id'] for s in cert.get('steps', [])
                if not reactions[s['reaction_id']]['enzyme_evidence_ids']})


def pairs_result(net, cdb_id, key):
    return next(p[key] for p in net['paired_probes'] if p['cannabisdb_id'] == cdb_id)
