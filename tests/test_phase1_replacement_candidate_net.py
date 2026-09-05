import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem

from cannabis_carbon.phase1_replacement_candidate_net import build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_all_target_records_and_exact_certificates_replay():
    report, parent, network = [read(n) for n in ('phase1-replacement-candidate-net', 'phase1-purine-candidate-net', 'phase1-full-balanced-network')]
    original = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    candidate_ids = report['candidate_reaction_evidence_ids']
    assert len(candidate_ids) == 1603
    assert all(candidate_ids[rid] == ids for rid, ids in parent['candidate_reaction_evidence_ids'].items())
    assert len(report['enzyme_evidence']) == 2
    assert all(e['full_search_report'] == 'phase1-replacement-search.json' for e in report['enzyme_evidence'])
    atoms = {}
    for rid in candidate_ids:
        balance = Counter()
        for side, sign in [('left', -1), ('right', 1)]:
            for m in original[rid][side]:
                cid = m['compound_id']
                if cid not in atoms:
                    mol = Chem.AddHs(Chem.MolFromSmiles(compounds[cid]['smiles']))
                    atoms[cid] = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                    atoms[cid]['charge'] = Chem.GetFormalCharge(mol)
                for element, n in atoms[cid].items():
                    balance[element] += sign * Fraction(m['coefficient']) * n
        assert all(n == 0 for n in balance.values())
    old = {c['compound_id']: c for c in parent['scenarios'][0]['certificates']}
    exchange = set(report['external_exchange_compound_ids'])
    assert {cid for cid in exchange if compounds[cid]['carbon_count']} == {report['co2_compound_id']}
    assert compounds[report['co2_compound_id']]['smiles'] == 'O=C=O'
    assert len(report['constraints']) == 6
    assert report['constraints'][:5] == parent['constraints']
    for scenario, count in zip(report['scenarios'], [109, 101]):
        assert len(scenario['targets']) == 6220
        assert [(t['cannabisdb_id'], t['compound_id']) for t in scenario['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
        assert scenario['summary']['target_status_counts'] == dict(Counter(t['net_status'] for t in scenario['targets']))
        assert scenario['summary']['target_status_counts']['exact-net-conversion-hypothesis'] == count
        assert scenario['summary']['new_target_ids'] == []
        assert not any(t['new_net_certificate'] for t in scenario['targets'])
        forbidden = set(scenario['forbidden_step_ids'])
        assert scenario['summary']['allowed_directed_steps'] == 2 * len(candidate_ids) - len(forbidden)
        assert {c['compound_id'] for c in scenario['certificates']} == set(scenario['preserved_certificate_compound_ids'])
        for cert in scenario['certificates']:
            assert cert == old[cert['compound_id']]
            net = defaultdict(Fraction)
            for step in cert['steps']:
                assert step['step_id'] not in forbidden
                assert step['reaction_id'] in candidate_ids
                assert step['step_id'] == step['reaction_id'] + ':' + step['direction_mode']
                assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
                assert Fraction(step['extent']) > 0
                sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
                for side, factor in [('left', -sign), ('right', sign)]:
                    for m in original[step['reaction_id']][side]:
                        net[m['compound_id']] += factor * Fraction(m['coefficient']) * Fraction(step['extent'])
            assert net[cert['compound_id']] == Fraction(cert['target_amount']) >= 1
            assert all(n >= 0 for cid, n in net.items() if cid not in exchange)
            assert cert['external_net_consumption'] == {c: str(-n) for c, n in net.items() if n < 0}
            assert cert['net_exports'] == {c: str(n) for c, n in net.items() if n > 0}
            assert cert['zero_net_internal_participants'] == sorted(c for c, n in net.items() if not n and c not in exchange)
            assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_assembly_replay_and_fail_closed_direction_and_annotation_checks(monkeypatch):
    inputs = [read(n) for n in ('phase1-full-balanced-network', 'phase1-purine-candidate-net', 'phase1-replacement-search', 'phase1-replacement-references')]
    report = read('phase1-replacement-candidate-net')
    outcomes = {}
    for scenario in report['scenarios']:
        rows = {t['compound_id']: {'status': t['net_status']} | {k: t[k] for k in ('solver_status', 'solver_message') if k in t} for t in scenario['targets']}
        rows.update({c['compound_id']: {k: v for k, v in c.items() if k != 'compound_id'} for c in scenario['certificates']})
        outcomes[scenario['summary']['allowed_directed_steps']] = rows
    monkeypatch.setattr('cannabis_carbon.phase1_purine_candidate_net.NetModel.solve', lambda self, cid: outcomes[len(self.steps)][cid])
    before = copy.deepcopy(inputs)
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    row = next(r for r in inputs[3]['rows'] if 'RHEA:16846' in r['source_reaction_ids'])
    next(s for s in row['sources'] if s['source_reaction_id'] == 'RHEA:16846')['source_left_corresponds_to'] = 'left'
    with pytest.raises(ValueError, match='canonical orientation changed'):
        build(*inputs)
    inputs = before
    inputs[2]['rows'][0]['reference_matches'] = []
    with pytest.raises(ValueError, match='Reference annotation join changed'):
        build(*inputs)
