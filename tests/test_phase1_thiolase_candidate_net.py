import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_thiolase_inventory_evidence_and_source_lineage():
    report = read('phase1-thiolase-candidate-net')
    parent = read('phase1-synthase-candidate-net')
    network = read('phase1-full-balanced-network')
    search = read('phase1-weighted-gap-search')
    ids = set(report['candidate_reaction_evidence_ids'])
    assert len(ids) == 1605
    assert ids - parent['candidate_reaction_evidence_ids'].keys() == {search['rows'][0]['reaction_id']}
    for rid, eids in parent['candidate_reaction_evidence_ids'].items():
        assert report['candidate_reaction_evidence_ids'][rid] == eids
    assert report['synthase_reference_links'] == parent['synthase_reference_links']
    assert report['constraints'] == parent['constraints']
    assert report['external_exchange_compound_ids'] == parent['external_exchange_compound_ids']
    assert len(report['enzyme_evidence']) == 1
    evidence = report['enzyme_evidence'][0]
    assert evidence['evidence_class'] == search['rows'][0]['evidence_class']
    assert {p['accession'] for p in evidence['screened_proteins']} == set(search['rows'][0]['screened_cannabis_proteins'])
    compounds = {c['id']: c for c in network['compounds']}
    assert all(balanced([r['left'], r['right']], compounds) for r in network['reactions'] if r['id'] in ids)
    expected = [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
    for scenario in report['scenarios']:
        assert len(scenario['targets']) == 6220
        assert [(t['cannabisdb_id'], t['compound_id']) for t in scenario['targets']] == expected
    assert len(report['probes']) == 8
    assert len(report['probe_results']) == 16
    assert {(r['scenario_id'], r['compound_id']) for r in report['probe_results']} == {
        (s['id'], c['id']) for s in report['scenarios'] for c in report['probes']}
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == sha


def test_all_thiolase_target_and_probe_certificates_balance_exactly():
    report = read('phase1-thiolase-candidate-net')
    network = read('phase1-full-balanced-network')
    compounds = {c['id']: c for c in network['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    assert {c for c in exchanges if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    steps = {s['id']: s for s in orientations([r for r in network['reactions'] if r['id'] in report['candidate_reaction_evidence_ids']])}
    for scenario in report['scenarios']:
        certs = {c['compound_id']: c for c in scenario['certificates']}
        for t in scenario['targets']:
            assert (t['compound_id'] in certs) == (t['net_status'] == 'exact-net-conversion-hypothesis')
        probes = [r for r in report['probe_results'] if r['scenario_id'] == scenario['id'] and r['status'] == 'exact-net-conversion-hypothesis']
        for cert in [*certs.values(), *probes]:
            assert not set(scenario['forbidden_step_ids']) & {s['step_id'] for s in cert['steps']}
            net = exact_net([steps[s['step_id']] for s in cert['steps']], [s['extent'] for s in cert['steps']])
            assert net[cert['compound_id']] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == cert['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == cert['net_exports']
            incoming = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
            assert incoming == sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0) > 0
