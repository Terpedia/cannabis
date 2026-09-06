import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from cannabis_carbon.phase1_carrier_reconstruction import RH, reconstruct


def test_restored_carriers_roundtrip_every_source_join():
    report = json.loads(Path('data/reports/phase1-carrier-reconstruction.json').read_bytes())
    audit = json.loads(Path('data/reports/phase1-carrier-catalog-audit.json').read_bytes())
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    originals = {e['id']: e for e in audit['equations']}
    expected = Counter((e['id'], json.dumps(j['source_record'], sort_keys=True))
                       for e in audit['equations'] for j in e['carrier_source_joins'])
    actual = Counter()
    carrier_ids = {c['id'] for c in report['carriers']}
    for equation in report['restored_equations']:
        assert not equation['full_macromolecular_balance_established']
        for join in equation['source_joins']:
            actual[join['model_reaction_id'], json.dumps(join['source_record'], sort_keys=True)] += 1
            projected = {s: Counter({p['compound_id']: p['coefficient'] for p in equation[s]})
                         for s in ('left', 'right')}
            for change in join['restorations']:
                side, cid, amount = change['side'], change['carrier_id'], change['coefficient']
                assert cid in carrier_ids
                assert change['source_occurrence']['carrier'] == cid
                projected[side][cid] -= amount
                assert projected[side][cid] >= 0
                for part in change['reactive_part_projection']:
                    assert part['projected_compound_id'] not in carrier_ids
                    projected[side][part['projected_compound_id']] += amount
            original = originals[join['model_reaction_id']]['model_reaction']
            for side in projected:
                assert +projected[side] == Counter({p['compound_id']: p['coefficient'] for p in original[side]})
    assert actual == expected
    assert sum(actual.values()) == 2739
    assert len(report['restored_equations']) == 1342
    assert len(originals) == 1316
    assert not report['summary']['corrected_pathway_coverage_established']
    # The cytochrome c and cL contexts must remain separate full species.
    used = {r['carrier_id'] for e in report['restored_equations']
            for j in e['source_joins'] for r in j['restorations']}
    assert {RH + 'Compound_' + n for n in ('12863', '12864', '14399', '10350')} <= used


def test_symbolic_coefficients_are_not_silently_made_numeric():
    audit = json.loads(Path('data/reports/phase1-carrier-catalog-audit.json').read_bytes())
    equation = audit['equations'][0]
    join = copy.deepcopy(equation['carrier_source_joins'][0])
    join['carrier_occurrences'][0]['coefficientPredicate'] = RH + 'contains2n'
    with pytest.raises(ValueError, match='symbolic-carrier-coefficient'):
        reconstruct(equation['model_reaction'], join, {}, {})
