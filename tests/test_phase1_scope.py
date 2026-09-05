import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_scope import expand, extend_network
from cannabis_carbon.balance import _reaction_smiles_balance


def reaction(rid, left, right):
    return {'id': rid, 'left': [{'compound_id': x, 'coefficient': n} for x, n in left],
            'right': [{'compound_id': x, 'coefficient': n} for x, n in right]}


def test_cyclic_missing_cofactor_cannot_bootstrap_from_adjacency():
    reactions = [reaction('r1', [('A', 1), ('B', 2)], [('C', 1)]),
                 reaction('r2', [('C', 1)], [('B', 2)])]
    assert expand(reactions, {'A'})['available'] == {'A'}
    result = expand(reactions + [reaction('r3', [('A', 1)], [('D', 1)]),
                                  reaction('r4', [('D', 1)], [('B', 1)])], {'A'})
    assert result['available'] == {'A', 'B', 'C', 'D'}
    witness = result['witnesses']['C']
    assert witness['level'] == 3
    assert witness['required_inputs'] == reactions[0]['left']
    assert len(result['enabled_orientations']) == len(set(result['enabled_orientations']))


def test_reverse_direction_is_labeled_and_empty_inputs_rejected():
    result = expand([reaction('r', [('A', 1)], [('B', 1)])], {'B'})
    assert result['witnesses']['A']['direction_mode'] == 'hypothetical-right-to-left'
    with pytest.raises(ValueError, match='nonempty'):
        expand([reaction('r', [], [('B', 1)])], set())


def test_full_network_rejects_unbalanced_and_generic_equations_and_deduplicates_reverse():
    source = [dict(rule_id=f'RHEA:{i}', reaction_smarts=smiles, source_url='https://www.rhea-db.org',
                   source_download_url='source', direction_mode='recorded', source_evidence_type='curated')
              for i, smiles in enumerate(['C=C.O>>CCO', 'CCO>>C=C.O', 'C>>CC', '*>>*'])]
    parent = {'compounds': [], 'reactions': [], 'targets': [], 'hypotheses': []}
    network = extend_network(parent, source)
    assert len(network['reactions']) == 1
    assert len(network['reactions'][0]['sources']) == 2
    assert {s['source_left_corresponds_to'] for s in network['reactions'][0]['sources']} == {'left', 'right'}
    assert network['summary']['rhea_balance_status_counts'] == {'balanced': 2, 'imbalanced': 1, 'not-auditable': 1}


def test_published_full_scope_witnesses_require_all_inputs_and_no_hidden_carbon_seeds():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-all-reactants-scope.json'
    report = json.loads(path.read_text())
    network_path = root / report['source_network']
    assert hashlib.sha256(network_path.read_bytes()).hexdigest() == report['source_network_sha256']
    network = json.loads(network_path.read_text())
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    for r in reactions.values():
        equation = '>>'.join('.'.join(compounds[m['compound_id']]['smiles'] for m in r[side] for _ in range(m['coefficient'])) for side in ('left', 'right'))
        element, charge = _reaction_smiles_balance(equation)
        assert element['status'] == charge['status'] == 'balanced'
    for scenario in report['scenarios']:
        seeds = set(scenario['seed_compound_ids'])
        assert [compounds[c]['smiles'] for c in seeds if compounds[c]['carbon_count']] == ['O=C=O']
        result = expand(network['reactions'], seeds)
        assert result['witnesses'] == scenario['witnesses']
        assert [t['cannabisdb_id'] for t in scenario['targets']] == [t['cannabisdb_id'] for t in network['targets']]
        for cid, witness in scenario['witnesses'].items():
            if cid in seeds:
                assert witness['source'] == 'explicit-seed'
                continue
            r = reactions[witness['reaction_id']]
            side = 'left' if witness['direction_mode'] == 'hypothetical-left-to-right' else 'right'
            assert witness['required_inputs'] == r[side]
            assert all(scenario['witnesses'][m['compound_id']]['level'] < witness['level'] for m in witness['required_inputs'])
        for target in scenario['targets']:
            assert (target['status'] != 'blocked') == (target['compound_id'] in result['available'])
            for step in target['blocked_producing_steps']:
                assert step['missing_inputs']
                assert all(m['compound_id'] not in result['available'] for m in step['missing_inputs'])
