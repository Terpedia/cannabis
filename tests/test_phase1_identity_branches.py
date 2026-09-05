import hashlib
import json
from fractions import Fraction
from pathlib import Path

from rdkit import Chem
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_identity_branches_preserve_both_sources_and_do_not_inherit_xrefs():
    report = read('phase1-identity-branches')
    review = read('phase1-identity-conflict-review')
    originals = {a['id']: a for a in review['assertions']}
    assert len(report['decisions']) == 4
    assert len(report['assertions']) == 8
    assert len(report['branches']) == 16
    assert len({b['id'] for b in report['branches']}) == 16
    for assertion in report['assertions']:
        assert assertion == originals[assertion['id']]
    for d in report['decisions']:
        assert d['historical_target_replaced'] is False
        assert d['inherited_source_xrefs'] is False
        assert len(d['retained_assertion_ids']) == 2
        for scenario in report['scenario_constraints']:
            assert {b['assertion_id'] for b in report['branches'] if b['cannabisdb_id'] == d['cannabisdb_id'] and b['scenario_id'] == scenario} == set(d['retained_assertion_ids'])
    choices = {d['cannabisdb_id']: d for d in report['decisions']}
    assert choices['CDB006156']['structure_fields']['computed_carbon_count'] == 3
    assert choices['CDB006156']['corroborating_pubchem_cids'] == [753]
    assert choices['CDB006169']['provisional_assertion_id'] is None
    assert choices['CDB000546']['provisional_assertion_id'].endswith(':sdf_derived_assertion')
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_branch_certificates_replay_exactly_and_catalog_producers_are_exhaustive():
    report = read('phase1-identity-branches')
    network = read('phase1-full-balanced-network')
    candidate = read('phase1-replacement-candidate-net')
    compounds = {c['id']: c for c in network['compounds']}
    evidence = candidate['candidate_reaction_evidence_ids']
    steps = {s['id']: s for s in orientations([r for r in network['reactions'] if r['id'] in evidence])}
    exchanges = set(report['external_exchange_compound_ids'])
    for b in report['branches']:
        cid = b['compound_id']
        if cid:
            assert Chem.MolToSmiles(Chem.MolFromSmiles(compounds[cid]['smiles'])) == b['canonical_smiles']
        expected = {}
        for r in network['reactions']:
            net = sum((sign * Fraction(m['coefficient']) for side, sign in [('left', -1), ('right', 1)] for m in r[side] if m['compound_id'] == cid), Fraction())
            if net:
                expected[r['id']] = net
        assert {p['reaction_id'] for p in b['catalog_producers']} == set(expected)
        for p in b['catalog_producers']:
            assert p['producing_orientation'] == ('left-to-right' if expected[p['reaction_id']] > 0 else 'right-to-left')
            assert Fraction(p['net_amount']) == abs(expected[p['reaction_id']])
            assert p['candidate_evidence_ids'] == evidence.get(p['reaction_id'], [])
        result = b['result']
        if result['status'] != 'exact-net-conversion-hypothesis':
            continue
        used = [steps[s['step_id']] for s in result['steps']]
        assert not {s['id'] for s in used} & set(report['scenario_constraints'][b['scenario_id']])
        net = exact_net(used, [s['extent'] for s in result['steps']])
        assert net[cid] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == result['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == result['net_exports']
        carbon_in = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
        carbon_out = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
        assert carbon_in == carbon_out == Fraction(result['net_carbon_in']) == Fraction(result['net_carbon_out'])
        assert {c for c, n in net.items() if n < 0 and compounds[c]['carbon_count']} == {candidate['co2_compound_id']}
