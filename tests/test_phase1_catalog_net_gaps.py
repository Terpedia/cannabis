"""Independently audit chemistry-only certificates, without rerunning the LP."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


@pytest.fixture(scope='module')
def inputs():
    return tuple(read(n) for n in ('phase1-catalog-net-gaps',
        'phase1-full-balanced-network', 'phase1-candidate-net-flux',
        'phase1-candidate-scope', 'phase1-all-reactants-scope'))


def test_catalog_inventory_and_pinned_sources(inputs):
    report, network, baseline, candidate, scope = inputs
    payload = (ROOT / 'data/reports/phase1-catalog-net-gaps.json').read_bytes()
    folder = ROOT / 'docs/data/catalog-net-view'
    assert payload == (folder / 'bundle.json').read_bytes()
    manifest = json.loads((folder / 'index.json').read_text())
    assert manifest['sha256'] == hashlib.sha256(payload).hexdigest()
    assert manifest['bytes'] == len(payload)
    assert manifest['summary'] == report['summary']
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    assert [(t['cannabisdb_id'], t['compound_id']) for t in report['targets']] == [
        (t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
    assert len(report['targets']) == 6220
    assert report['summary']['target_status_counts'] == dict(Counter(t['net_status'] for t in report['targets']))
    assert report['summary']['catalog_equations'] == len(network['reactions'])
    assert report['summary']['baseline_candidate_equations'] == len(candidate['candidate_reaction_evidence_ids'])
    assert report['external_exchange_compound_ids'] == baseline['external_exchange_compound_ids']
    scenario = next(s for s in scope['scenarios'] if s['id'] == 'CO2-plus-all-carbon-free-species')
    assert set(report['external_exchange_compound_ids']) == set(scenario['seed_compound_ids'])
    assert [t['startup_status'] for t in report['targets']] == [t['status'] for t in scenario['targets']]


def test_every_selected_equation_is_exact_and_balanced(inputs):
    report, network, _, candidate, _ = inputs
    original = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    atoms = {}
    for r in report['reactions']:
        assert r['left'] == original[r['id']]['left']
        assert r['right'] == original[r['id']]['right']
        assert r['enzyme_evidence_ids'] == candidate['candidate_reaction_evidence_ids'].get(r['id'], [])
        assert r['missing_candidate_evidence'] == (not r['enzyme_evidence_ids'])
        balance = Counter()
        for side, sign in [('left', -1), ('right', 1)]:
            for m in r[side]:
                cid, coefficient = m['compound_id'], Fraction(m['coefficient'])
                assert coefficient > 0
                if cid not in atoms:
                    mol = Chem.AddHs(Chem.MolFromSmiles(compounds[cid]['smiles']))
                    atoms[cid] = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                    atoms[cid]['charge'] = Chem.GetFormalCharge(mol)
                for key, count in atoms[cid].items():
                    balance[key] += sign * coefficient * count
        assert all(value == 0 for value in balance.values())


def test_all_certificates_close_without_organic_depletion(inputs):
    report, network, baseline, _, _ = inputs
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    assert [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] == ['O=C=O']
    old = {c['compound_id']: c for c in baseline['certificates']}
    for cert in report['certificates']:
        if cert['compound_id'] in old:
            for key, value in old[cert['compound_id']].items():
                assert cert[key] == value
            assert cert['certificate_origin'] == 'unchanged-baseline-witness'
        net = defaultdict(Fraction)
        for step in cert['steps']:
            assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
            sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
            extent = Fraction(step['extent'])
            assert extent > 0
            for side, multiplier in [('left', -sign), ('right', sign)]:
                for m in reactions[step['reaction_id']][side]:
                    net[m['compound_id']] += multiplier * extent * Fraction(m['coefficient'])
        assert net[cert['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == cert['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == cert['net_exports']
        assert sorted(c for c, n in net.items() if not n and c not in exchanges) == cert['zero_net_internal_participants']
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
        assert -net[report['co2_compound_id']] == Fraction(cert['net_carbon_in']) == Fraction(cert['net_carbon_out'])


def test_gap_ranking_is_selected_certificate_membership_not_necessity(inputs):
    report, _, _, candidate, _ = inputs
    supported = candidate['candidate_reaction_evidence_ids']
    certs = {c['compound_id']: c for c in report['certificates']}
    members, structures = defaultdict(list), defaultdict(set)
    for target in report['targets']:
        cert = certs.get(target['certificate_compound_id'])
        assert bool(cert) == (target['net_status'] == 'exact-net-conversion-hypothesis')
        missing = sorted({s['reaction_id'] for s in cert['steps'] if not supported.get(s['reaction_id'])}) if cert else []
        assert target['missing_candidate_reaction_ids'] == missing
        if cert:
            assert cert['missing_candidate_reaction_ids'] == missing
        for rid in missing:
            members[rid].append(target['cannabisdb_id'])
            structures[rid].add(target['compound_id'])
    assert {g['id'] for g in report['gap_priorities']} == set(members)
    for gap in report['gap_priorities']:
        assert gap['selected_certificate_target_ids'] == members[gap['id']]
        assert gap['selected_certificate_target_count'] == len(members[gap['id']])
        assert gap['selected_certificate_structure_count'] == len(structures[gap['id']])
    assert report['summary']['selected_missing_candidate_equations'] == len(members)


def test_no_producer_status_against_full_net_stoichiometry(inputs):
    report, network, _, _, _ = inputs
    producing = set()
    for r in network['reactions']:
        net = defaultdict(Fraction)
        for side, sign in [('left', -1), ('right', 1)]:
            for m in r[side]:
                net[m['compound_id']] += sign * Fraction(m['coefficient'])
        producing.update(c for c, n in net.items() if n)  # both hypothetical directions
    for target in report['targets']:
        if target['net_status'] == 'no-net-producing-catalog-equation':
            assert target['compound_id'] not in producing
        if target['net_status'] == 'solver-reported-infeasible':
            assert target['compound_id'] in producing
            assert target['solver_status'] == 2  # numerical status, not exact impossibility
