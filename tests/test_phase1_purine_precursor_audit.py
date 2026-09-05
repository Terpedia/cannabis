import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem

from cannabis_carbon.phase1_purine_precursor_audit import build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_every_probe_and_certificate_preserves_exact_chemistry():
    report = read('phase1-purine-precursor-audit')
    network = read('phase1-full-balanced-network')
    compounds = {c['id']: c for c in report['compounds']}
    originals = {r['id']: r for r in network['reactions']}
    reactions = {r['id']: r for r in report['reactions']}
    assert all(r == originals[rid] for rid, r in reactions.items())
    atoms = {}
    for rid, reaction in reactions.items():
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
        assert all(v == 0 for v in balance.values())
    exchange = set(report['external_exchange_compound_ids'])
    assert {c for c in exchange if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    forbidden = {c['id'] for c in report['constraints']}
    probes = set(report['participant_compound_ids']) | {t['compound_id'] for t in report['focused_targets']}
    assert len(report['participant_compound_ids']) == 32
    assert len(report['focused_targets']) == 7
    assert len(report['plant_family_audit']) == 13
    for scenario in report['scenarios']:
        assert {r['compound_id'] for r in scenario['results']} == probes
        assert len(scenario['results']) == len(probes)
        allowed = set(report['baseline_candidate_reaction_ids'])
        if scenario['id'] == 'restricted-plus-plant-hypotheses':
            allowed.update(report['added_candidate_reaction_ids'])
        if scenario['id'] == 'restricted-full-catalog-chemistry':
            allowed = set(originals)
        assert scenario['equation_count'] == len(allowed)
        assert scenario['allowed_directed_steps'] == 2 * len(allowed) - len(forbidden)
        for result in scenario['results']:
            if result['status'] != 'exact-net-conversion-hypothesis':
                continue
            net = defaultdict(Fraction)
            for step in result['steps']:
                assert step['reaction_id'] in allowed
                assert step['step_id'] not in forbidden
                assert step['step_id'] == step['reaction_id'] + ':' + step['direction_mode']
                assert Fraction(step['extent']) > 0
                assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
                sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
                for side, mult in [('left', -sign), ('right', sign)]:
                    for m in reactions[step['reaction_id']][side]:
                        net[m['compound_id']] += mult * Fraction(step['extent']) * Fraction(m['coefficient'])
            assert net[result['compound_id']] == Fraction(result['target_amount']) >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchange)
            assert result['external_net_consumption'] == {c: str(-n) for c, n in net.items() if n < 0}
            assert result['net_exports'] == {c: str(n) for c, n in net.items() if n > 0}
            assert result['zero_net_internal_participants'] == sorted(c for c, n in net.items() if not n and c not in exchange)
            assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    assert all('unreviewed' in e['evidence_class'] for e in report['enzyme_evidence'])
    assert len(report['catalog_candidate_gaps']) == report['summary']['selected_catalog_candidate_gaps'] == 70
    candidate_ids = set(report['baseline_candidate_reaction_ids']) | set(report['added_candidate_reaction_ids'])
    catalog = {r['compound_id']: r for r in report['scenarios'][2]['results']}
    for gap in report['catalog_candidate_gaps']:
        assert gap['reaction_id'] not in candidate_ids
        assert gap['source_joins'] == originals[gap['reaction_id']]['sources']
        assert all(gap[side] == originals[gap['reaction_id']][side] for side in ('left', 'right'))
        for use in gap['selected_uses']:
            assert {k: v for k, v in use.items() if k != 'probe_compound_id'} in catalog[use['probe_compound_id']]['steps']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_replay_nonmutation_and_exchange_guard(monkeypatch):
    inputs = [read(n) for n in ('phase1-full-balanced-network', 'phase1-expanded-candidate-net',
              'phase1-plant-purine-references', 'phase1-plant-purine-search', 'phase1-direction-sensitivity')]
    report = read('phase1-purine-precursor-audit')
    # Captured outcomes avoid asserting a unique minimum-extent solver route.
    # Every successful certificate is independently replayed above.
    by_step_count = {s['allowed_directed_steps']: {r['compound_id']: {k: v for k, v in r.items() if k != 'compound_id'}
                     for r in s['results']} for s in report['scenarios']}
    monkeypatch.setattr('cannabis_carbon.phase1_purine_precursor_audit.NetModel.solve',
                        lambda self, cid: by_step_count[len(self.steps)][cid])
    before = copy.deepcopy(inputs)
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    inputs[1]['external_exchange_compound_ids'].append(next(c['id'] for c in inputs[0]['compounds']
                                                         if c['carbon_count'] > 1))
    with pytest.raises(ValueError, match='Exchange boundary changed'):
        build(*inputs)
