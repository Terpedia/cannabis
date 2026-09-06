import hashlib
import json
import re
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from rdkit import Chem
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_odd_chain_elongation import build, REFERENCES


def test_exact_homolog_cycles_and_net_cofactor_bookkeeping():
    read = lambda n: json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
    parent, network, report = map(read, ('alkane-net', 'full-balanced-network', 'odd-chain-elongation'))
    assert build(parent, network) == {k: v for k, v in report.items() if k != 'source_sha256'}
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert report['summary'] == {'cycles': 2, 'balanced_equations': 8,
                                 'new_compound_structures': 6, 'coverage_gain_claimed': 0}
    original = {c['id']: c for c in parent['compounds']}
    compounds = {c['id']: c for c in report['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    refs = {r['id']: r for r in report['reference_reactions']}
    endpoints = []
    for cycle in report['cycles']:
        shift = cycle['precursor_acyl_carbons'] - 18
        mapping = {t['reference_compound_id']: t['compound_id'] for t in cycle['transformations']}
        assert len(mapping) == 5
        for t in cycle['transformations']:
            sm = original[t['reference_compound_id']]['smiles']
            n = len(re.match(r'^C+', sm).group())
            expected = Chem.MolToSmiles(Chem.MolFromSmiles('C' * (n + shift) + sm[n:]), isomericSmiles=True)
            actual = compounds[t['compound_id']]
            assert actual['smiles'] == expected
            assert actual['id'] == stable_id('structure', expected)
            assert actual['carbon_count'] == original[t['reference_compound_id']]['carbon_count'] + shift
            assert actual['formal_charge'] == original[t['reference_compound_id']]['formal_charge']
            assert t['present_in_parent_model'] == (t['compound_id'] in original)
        net = defaultdict(Fraction)
        reference_net = defaultdict(Fraction)
        for i, rid in enumerate(cycle['reaction_ids_in_order']):
            r = reactions[rid]; ref = refs[REFERENCES[i]]
            assert r['reference_reaction_id'] == ref['id']
            assert r['direction_status'] == 'proposed-forward-only'
            assert r['enzyme_evidence_ids'] == []
            assert balanced([r['left'], r['right']], compounds)
            for side, sign in (('left', -1), ('right', 1)):
                assert r[side] == [{**p, 'compound_id': mapping.get(p['compound_id'], p['compound_id'])} for p in ref[side]]
                for p in r[side]:
                    net[p['compound_id']] += sign * Fraction(str(p['coefficient']))
                for p in ref[side]:
                    reference_net[p['compound_id']] += sign * Fraction(str(p['coefficient']))
        net = {c: v for c, v in net.items() if v}
        expected_net = {mapping.get(c, c): v for c, v in reference_net.items() if v}
        assert net == expected_net
        # Only starter and final acyl-CoA survive summation; all three intermediates cancel.
        chain_net = {c: v for c, v in net.items() if c in mapping.values()}
        assert sorted(chain_net.values()) == [-1, 1]
        endpoints.append((next(c for c, v in chain_net.items() if v == -1),
                          next(c for c, v in chain_net.items() if v == 1)))
        cofactor_net = {c: v for c, v in net.items() if c not in mapping.values()}
        assert len(cofactor_net) == 7
        assert sorted(cofactor_net.values()) == [-3, -2, -1, 1, 1, 1, 2]
    assert endpoints[0][1] == endpoints[1][0]
