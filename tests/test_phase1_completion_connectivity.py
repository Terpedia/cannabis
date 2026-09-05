import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_completion_connectivity import assemble
from cannabis_carbon.phase1_completion_net_view import build as build_view
from cannabis_carbon.phase1_net_flux import NetModel
from cannabis_carbon.phase1_scope import expand

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


@pytest.fixture(scope='module')
def inputs():
    return [read(n) for n in ['phase1-full-balanced-network', 'phase1-candidate-scope',
        'phase1-marts-completions', 'phase1-completion-protein-evidence']]


def test_exact_admission_and_no_source_mutation(inputs):
    before = [json.dumps(r, sort_keys=True) for r in inputs]
    reactions, compounds, added, excluded = assemble(*inputs)
    assert before == [json.dumps(r, sort_keys=True) for r in inputs]
    baseline = inputs[1]['candidate_reaction_evidence_ids']
    completions = {h['id']: h for h in inputs[2]['completions']}
    evidence = {e['id']: e for e in inputs[3]['rows']}
    assert len(reactions) == 1793 and len(added) == 321 and len(excluded) == 380
    assert sum(r['baseline_balanced_equation_exists'] for r in added) == 146
    assert {r['id'] for r in reactions} == set(baseline) | {r['id'] for r in added}
    for r in added:
        h = completions[r['completion_id']]; e = evidence[r['completion_id']]
        assert [r['left'], r['right']] == [h['left'], h['right']]
        assert r['screened_cannabis_proteins'] == e['screened_cannabis_proteins']
        assert r['validation_blockers'] == e['validation_blockers']
        assert r['prior_source_reviews'] == e['prior_source_reviews']
        assert 'enzyme_evidence_ids' not in r  # no promotion to baseline enzyme evidence
        assert r['id'] not in baseline
    broken = copy.deepcopy(inputs[3]); broken['rows'][0]['reaction_id'] = 'wrong'
    with pytest.raises(ValueError, match='identity mismatch'):
        assemble(*inputs[:3], broken)
    altered = copy.deepcopy(inputs[2]); altered['completions'][0]['left'][0]['coefficient'] += 1
    with pytest.raises(ValueError, match='identity mismatch'):
        assemble(inputs[0], inputs[1], altered, inputs[3])


def test_full_inventory_startup_and_every_changed_net_result(inputs):
    report = read('phase1-completion-connectivity'); baseline_net = read('phase1-candidate-net-flux')
    path = ROOT / 'data/reports/phase1-completion-connectivity.json'
    assert path.read_bytes() == (ROOT / 'docs/data/phase1-completion-connectivity.json').read_bytes()
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    reactions, compounds, added, excluded = assemble(*inputs)
    assert added == report['admitted_reactions'] and excluded == report['excluded_completions']
    assert len(report['targets']) == 6220
    assert [t['cannabisdb_id'] for t in report['targets']] == [t['cannabisdb_id'] for t in baseline_net['targets']]
    for old, scenario in zip(inputs[1]['scenarios'], report['startup_scenarios']):
        assert old['seed_compound_ids'] == scenario['seed_compound_ids']
        replay = expand(reactions, set(scenario['seed_compound_ids']))
        assert replay['witnesses'] == scenario['witnesses']
        assert not scenario['newly_available_vs_baseline'] and not scenario['rescued_target_ids']
        for cid, w in scenario['witnesses'].items():
            if w['level']:
                assert all(scenario['witnesses'][m['compound_id']]['level'] < w['level'] for m in w['required_inputs'])
    assert report['external_exchange_compound_ids'] == baseline_net['external_exchange_compound_ids']
    assert [compounds[c]['smiles'] for c in report['external_exchange_compound_ids'] if compounds[c]['carbon_count']] == ['O=C=O']
    model = NetModel(reactions, set(report['external_exchange_compound_ids']))
    cache = {}
    for old, target in zip(baseline_net['targets'], report['targets']):
        cid = target['compound_id']
        assert target['baseline_net_status'] == old['net_status']
        if old['certificate_compound_id']:
            assert target['certificate'] == {'report': 'phase1-candidate-net-flux.json', 'compound_id': cid}
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)['status']
            assert target['sensitivity_net_status'] == cache[cid]
    counts = Counter(t['sensitivity_net_status'] for t in report['targets'])
    assert dict(counts) == report['summary']['sensitivity_net_status_counts']
    assert counts['exact-net-conversion-hypothesis'] == 112


def test_new_certificates_independently_close_elements_isotopes_charge_and_pools(inputs):
    report = read('phase1-completion-connectivity')
    reactions = {r['id']: r for r in report['admitted_reactions'] + report['certificate_baseline_reactions']}
    compounds = {c['id']: c for c in report['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    compositions = {}
    for cid, c in compounds.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(c['smiles']))
        compositions[cid] = (Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms()), Chem.GetFormalCharge(mol))
    for r in reactions.values():
        atoms = defaultdict(int); charge = 0
        for side, sign in [('left', -1), ('right', 1)]:
            for m in r[side]:
                assert isinstance(m['coefficient'], int) and m['coefficient'] > 0
                composition, q = compositions[m['compound_id']]
                for atom, count in composition.items():
                    atoms[atom] += sign * m['coefficient'] * count
                charge += sign * m['coefficient'] * q
        assert not any(atoms.values()) and charge == 0
    assert len(report['additional_net_certificates']) == 11
    added_ids = {r['id'] for r in report['admitted_reactions']}
    for cert in report['additional_net_certificates']:
        net = defaultdict(Fraction)
        assert any(s['reaction_id'] in added_ids for s in cert['steps'])
        for step in cert['steps']:
            reaction = reactions[step['reaction_id']]
            forward = step['direction_mode'] == 'hypothetical-left-to-right'
            assert step['direction_mode'] in ['hypothetical-left-to-right', 'hypothetical-right-to-left']
            extent = Fraction(step['extent']); assert extent > 0
            for side, sign in [('left', -1 if forward else 1), ('right', 1 if forward else -1)]:
                for m in reaction[side]:
                    net[m['compound_id']] += sign * extent * m['coefficient']
        assert net[cert['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == cert['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == cert['net_exports']
        assert sorted(c for c, n in net.items() if not n and c not in exchanges) == cert['zero_net_internal_participants']
        assert -net[report['co2_compound_id']] == Fraction(cert['net_carbon_in']) == Fraction(cert['net_carbon_out'])
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
        assert sum(n * compounds[c]['formal_charge'] for c, n in net.items()) == 0


def test_sensitivity_view_rebuild_preserves_all_baseline_certificates_and_evidence():
    folder = ROOT / 'docs/data/completion-net-view'
    manifest = json.loads((folder / 'index.json').read_text())
    raw = (folder / 'bundle.json').read_bytes(); bundle = json.loads(raw)
    assert len(raw) == manifest['bytes'] < 9_000_000
    assert hashlib.sha256(raw).hexdigest() == manifest['sha256']
    for path, digest in manifest['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    sources = [read(n) for n in ['phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay']]
    baseline = read('phase1-candidate-net-flux')
    assert build_view(read('phase1-completion-connectivity'), baseline, sources,
        read('phase1-completion-protein-evidence'), read('phase1-marts-audit')) == bundle
    assert bundle['certificates'][:len(baseline['certificates'])] == baseline['certificates']
    assert len(bundle['certificates']) == 111  # 112 target records, one exact-structure alias
    evidence = {e['id']: e for e in bundle['enzyme_evidence']}
    for r in bundle['reactions']:
        assert all(eid in evidence for eid in r['enzyme_evidence_ids'])
        if r.get('is_completion_sensitivity'):
            assert all(evidence[eid]['evidence_class'] == 'MARTS-source-homology-for-inferred-stoichiometry' for eid in r['enzyme_evidence_ids'])
