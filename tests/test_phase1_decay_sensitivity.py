import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem

from cannabis_carbon.phase1_decay_sensitivity import build

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT / path).read_text())


def test_all_replacement_certificates_preserve_balance_boundary_and_direction():
    report = read('data/reports/phase1-decay-sensitivity.json')
    network = read('data/reports/phase1-full-balanced-network.json')
    precursor = read('data/reports/phase1-purine-precursor-audit.json')
    originals = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in report['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    exchange = set(report['external_exchange_compound_ids'])
    forbidden = {c['id'] for c in report['constraints']}
    assert len(forbidden) == 6
    assert {c['id'] for c in precursor['constraints']} < forbidden
    assert report['summary']['allowed_directed_steps'] == 2 * len(originals) - 6
    assert {c for c in exchange if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    assert compounds[report['co2_compound_id']]['smiles'] == 'O=C=O'
    assert len(report['results']) == 39
    assert {r['compound_id'] for r in report['results']} == set(precursor['participant_compound_ids']) | {
        t['compound_id'] for t in precursor['focused_targets']}
    atoms = {}
    for rid, reaction in reactions.items():
        assert reaction == originals[rid]
        balance = Counter()
        for side, sign in [('left', -1), ('right', 1)]:
            for m in reaction[side]:
                cid = m['compound_id']
                if cid not in atoms:
                    mol = Chem.AddHs(Chem.MolFromSmiles(compounds[cid]['smiles']))
                    atoms[cid] = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                    atoms[cid]['charge'] = Chem.GetFormalCharge(mol)
                for element, n in atoms[cid].items():
                    balance[element] += sign * Fraction(m['coefficient']) * n
        assert all(n == 0 for n in balance.values())
    for result in report['results']:
        if result['status'] != 'exact-net-conversion-hypothesis':
            continue
        net = defaultdict(Fraction)
        for step in result['steps']:
            assert step['step_id'] not in forbidden
            assert step['step_id'] == step['reaction_id'] + ':' + step['direction_mode']
            assert Fraction(step['extent']) > 0
            assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
            sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
            for side, factor in [('left', -sign), ('right', sign)]:
                for m in reactions[step['reaction_id']][side]:
                    net[m['compound_id']] += factor * Fraction(m['coefficient']) * Fraction(step['extent'])
        assert net[result['compound_id']] == Fraction(result['target_amount']) >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchange)
        assert result['external_net_consumption'] == {c: str(-n) for c, n in net.items() if n < 0}
        assert result['net_exports'] == {c: str(n) for c, n in net.items() if n > 0}
        assert result['zero_net_internal_participants'] == sorted(c for c, n in net.items() if not n and c not in exchange)
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    candidates = read('data/reports/phase1-purine-candidate-net.json')['candidate_reaction_evidence_ids']
    results = {r['compound_id']: r for r in report['results']}
    assert report['summary']['alternative_certificates'] == 20
    assert report['summary']['selected_missing_candidate_equations'] == len(report['candidate_gaps']) == 66
    for gap in report['candidate_gaps']:
        assert gap['reaction_id'] not in candidates
        assert gap['source_joins'] == originals[gap['reaction_id']]['sources']
        for use in gap['selected_uses']:
            assert {k: v for k, v in use.items() if k != 'probe_compound_id'} in results[use['probe_compound_id']]['steps']


def test_assembly_replay_preserves_inputs_and_rejects_changed_boundary(monkeypatch):
    report = read('data/reports/phase1-decay-sensitivity.json')
    inputs = [read(p) for p in report['source_sha256']]
    captured = {r['compound_id']: {k: v for k, v in r.items()
                if k not in ('compound_id', 'previous_status', 'previous_used_excluded_step')}
                for r in report['results']}
    monkeypatch.setattr('cannabis_carbon.phase1_decay_sensitivity.NetModel.solve', lambda self, cid: captured[cid])
    before = copy.deepcopy(inputs)
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    inputs[2]['selected_canonical_direction'] = 'hypothetical-right-to-left'
    with pytest.raises(ValueError, match='Reviewed orientation changed'):
        build(*inputs)
    inputs = before
    inputs[1]['external_exchange_compound_ids'] = []
    with pytest.raises(ValueError, match='Exchange boundary mismatch'):
        build(*inputs)
