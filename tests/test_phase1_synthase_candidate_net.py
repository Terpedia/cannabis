import hashlib
import json
from fractions import Fraction
from pathlib import Path

from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_synthase_model_preserves_full_target_inventory_and_exact_equations():
    report = read('phase1-synthase-candidate-net')
    parent = read('phase1-replacement-candidate-net')
    network = read('phase1-full-balanced-network')
    ids = set(report['candidate_reaction_evidence_ids'])
    assert len(ids) == 1604
    assert len(ids - parent['candidate_reaction_evidence_ids'].keys()) == 1
    assert parent['candidate_reaction_evidence_ids'].keys() <= ids
    for rid, evidence in parent['candidate_reaction_evidence_ids'].items():
        assert set(evidence) <= set(report['candidate_reaction_evidence_ids'][rid])
    compounds = {c['id']: c for c in network['compounds']}
    assert all(balanced([r['left'], r['right']], compounds) for r in network['reactions'] if r['id'] in ids)
    expected = [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
    for scenario in report['scenarios']:
        assert [(t['cannabisdb_id'], t['compound_id']) for t in scenario['targets']] == expected
        assert len(scenario['targets']) == 6220
    restricted = report['scenarios'][1]
    assert len(restricted['forbidden_step_ids']) == 8
    assert set(parent['scenarios'][1]['forbidden_step_ids']) <= set(restricted['forbidden_step_ids'])
    for link in report['synthase_reference_links']:
        assert link['reaction_id'] in ids and link['core_reaction_id'] not in ids
        assert not link['core_identity_merge_allowed']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_every_synthase_model_certificate_has_exact_net_co2_balance():
    report = read('phase1-synthase-candidate-net')
    network = read('phase1-full-balanced-network')
    compounds = {c['id']: c for c in network['compounds']}
    steps = {s['id']: s for s in orientations([r for r in network['reactions'] if r['id'] in report['candidate_reaction_evidence_ids']])}
    exchanges = set(report['external_exchange_compound_ids'])
    assert {c for c in exchanges if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    for scenario in report['scenarios']:
        certificates = {c['compound_id']: c for c in scenario['certificates']}
        for target in scenario['targets']:
            assert (target['compound_id'] in certificates) == (target['net_status'] == 'exact-net-conversion-hypothesis')
        for cid, cert in certificates.items():
            assert not {s['step_id'] for s in cert['steps']} & set(scenario['forbidden_step_ids'])
            net = exact_net([steps[s['step_id']] for s in cert['steps']], [s['extent'] for s in cert['steps']])
            assert net[cid] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == cert['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == cert['net_exports']
            incoming = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
            outgoing = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
            assert incoming == outgoing > 0
            assert Fraction(cert['external_net_consumption'][report['co2_compound_id']]) == incoming
