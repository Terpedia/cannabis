import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem
from cannabis_carbon.phase1_purine_candidate_net import build
from cannabis_carbon.phase1_purine_candidate_view import build as view

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_all_targets_and_certificates_independently_replay():
    report = read('phase1-purine-candidate-net')
    network = read('phase1-full-balanced-network')
    baseline = read('phase1-expanded-candidate-net')
    original = {r['id']: r for r in network['reactions']}
    reactions = {r['id']: r for r in report['reactions']}
    compounds = {c['id']: c for c in report['compounds']}
    atoms = {}
    for rid, r in reactions.items():
        assert r['enzyme_evidence_ids'] == report['candidate_reaction_evidence_ids'][rid]
        assert all(r[k] == original[rid][k] for k in ('left', 'right', 'sources'))
        balance = Counter()
        for side, sign in [('left', -1), ('right', 1)]:
            for m in r[side]:
                cid = m['compound_id']
                if cid not in atoms:
                    mol = Chem.AddHs(Chem.MolFromSmiles(compounds[cid]['smiles']))
                    atoms[cid] = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                    atoms[cid]['charge'] = Chem.GetFormalCharge(mol)
                for element, n in atoms[cid].items():
                    balance[element] += sign * Fraction(m['coefficient']) * n
        assert all(n == 0 for n in balance.values())
    assert report['summary']['candidate_equations'] == 1601
    assert len(report['enzyme_evidence']) == 13
    assert sum('unreviewed' in e['evidence_class'] for e in report['enzyme_evidence']) == 3
    assert report['summary']['added_distinct_proteins'] == 33
    exchange = set(report['external_exchange_compound_ids'])
    assert {c for c in exchange if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    old = {c['compound_id']: c for c in baseline['certificates']}
    for scenario in report['scenarios']:
        assert [(t['cannabisdb_id'], t['compound_id']) for t in scenario['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
        forbidden = set(scenario['forbidden_step_ids'])
        certs = {c['compound_id']: c for c in scenario['certificates']}
        assert all(certs[cid] == old[cid] for cid in scenario['preserved_certificate_compound_ids'])
        assert scenario['summary']['target_status_counts'] == dict(Counter(t['net_status'] for t in scenario['targets']))
        for cert in certs.values():
            net = defaultdict(Fraction)
            for step in cert['steps']:
                assert step['reaction_id'] in report['candidate_reaction_evidence_ids']
                assert step['step_id'] not in forbidden
                assert step['step_id'] == step['reaction_id'] + ':' + step['direction_mode']
                assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
                sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
                assert Fraction(step['extent']) > 0
                for side, multiplier in [('left', -sign), ('right', sign)]:
                    for m in reactions[step['reaction_id']][side]:
                        net[m['compound_id']] += multiplier * Fraction(step['extent']) * Fraction(m['coefficient'])
            assert net[cert['compound_id']] == Fraction(cert['target_amount']) >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchange)
            assert cert['external_net_consumption'] == {c: str(-n) for c, n in net.items() if n < 0}
            assert cert['net_exports'] == {c: str(n) for c, n in net.items() if n > 0}
            assert cert['zero_net_internal_participants'] == sorted(c for c, n in net.items() if not n and c not in exchange)
            assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    assert report['scenarios'][0]['summary']['new_target_ids'] == ['CDB004932']
    assert report['scenarios'][1]['summary']['new_target_ids'] == []
    target = next(t for t in report['scenarios'][0]['targets'] if t['new_net_certificate'])
    added_route = next(c for c in report['scenarios'][0]['certificates'] if c['compound_id'] == target['compound_id'])
    new_ids = {e['reaction_id'] for e in report['enzyme_evidence']}
    unreviewed = {e['reaction_id'] for e in report['enzyme_evidence'] if 'unreviewed' in e['evidence_class']}
    assert sum(s['reaction_id'] in new_ids for s in added_route['steps']) == 4
    assert sum(s['reaction_id'] in unreviewed for s in added_route['steps']) == 2
    assert len({s['step_id'] for s in added_route['steps']} & set(report['scenarios'][1]['forbidden_step_ids'])) == 2
    assert [s['summary']['target_status_counts']['exact-net-conversion-hypothesis'] for s in report['scenarios']] == [109, 101]
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_rebuild_preserves_sources_and_rejects_overlap(monkeypatch):
    inputs = [read(n) for n in ('phase1-full-balanced-network', 'phase1-expanded-candidate-net',
        'phase1-direction-sensitivity', 'phase1-plant-purine-search', 'phase1-purine-gap-search')]
    report = read('phase1-purine-candidate-net')
    outcomes = {}
    for s in report['scenarios']:
        rows = {t['compound_id']: {'status': t['net_status']} | {k: t[k] for k in ('solver_status', 'solver_message') if k in t} for t in s['targets']}
        rows.update({c['compound_id']: {k: v for k, v in c.items() if k != 'compound_id'} for c in s['certificates']})
        outcomes[s['summary']['allowed_directed_steps']] = rows
    monkeypatch.setattr('cannabis_carbon.phase1_purine_candidate_net.NetModel.solve', lambda self, cid: outcomes[len(self.steps)][cid])
    before = copy.deepcopy(inputs)
    actual = build(*inputs)
    assert {k: v for k, v in actual.items() if k != 'scipy_version'} == {k: v for k, v in report.items() if k not in ('source_sha256', 'scipy_version')}
    assert inputs == before
    inputs[1]['candidate_reaction_evidence_ids'][report['enzyme_evidence'][0]['reaction_id']] = ['duplicate']
    with pytest.raises(ValueError, match='Supplement overlaps'):
        build(*inputs)


def test_static_view_keeps_both_scenarios_and_reference_confidence():
    names = ('phase1-purine-candidate-net', 'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay',
             'phase1-route-enzyme-overlay', 'phase1-expanded-candidate-net', 'phase1-candidate-direction-review')
    reports = [read(n) for n in names]
    folder = ROOT / 'docs/data/purine-net-view'
    bundle = json.loads((folder / 'bundle.json').read_text())
    manifest = json.loads((folder / 'index.json').read_text())
    assert view(reports[0], reports[1:-1], reports[-1]) == bundle
    assert hashlib.sha256((folder / manifest['file']).read_bytes()).hexdigest() == manifest['sha256']
    for path, sha in manifest['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    evidence = {e['id']: e for e in bundle['enzyme_evidence']}
    assert all(evidence[e['id']] == e for e in reports[0]['enzyme_evidence'])
    assert len(bundle['targets']) == len(bundle['restricted_scenario']['targets']) == 6220
    assert all('not recomputed' in t['startup_status'] for t in bundle['targets'])
    assert len([r for r in bundle['reactions'] if r['new_purine_candidate']]) == 13
