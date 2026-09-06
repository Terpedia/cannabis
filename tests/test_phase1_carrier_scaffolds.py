import json
from collections import Counter
from pathlib import Path

from cannabis_carbon.phase1_carrier_scaffolds import RH, audit_equation, scaffold_candidate


def test_entire_reconstruction_has_explicit_scaffold_status():
    source = json.loads(Path('data/reports/phase1-carrier-reconstruction.json').read_bytes())
    report = json.loads(Path('data/reports/phase1-carrier-scaffolds.json').read_bytes())
    candidates = {c['carrier_id']: c for c in report['carrier_scaffold_candidates']}
    assert report['carrier_scaffold_candidates'] == [v for c in source['carriers'] if (v := scaffold_candidate(c))]
    assert len(report['equations']) == len(source['restored_equations']) == 1342
    for actual, original in zip(report['equations'], source['restored_equations']):
        assert actual == audit_equation(original, candidates)
        assert not actual['full_macromolecular_balance_established']
    assert dict(Counter(e['status'] for e in report['equations'])) == report['summary']['equation_status_counts']
    assert candidates[RH + 'Compound_12863']['symbolic_scaffold'] == candidates[RH + 'Compound_12864']['symbolic_scaffold']
    assert candidates[RH + 'Compound_14399']['symbolic_scaffold'] == candidates[RH + 'Compound_10350']['symbolic_scaffold']
    assert candidates[RH + 'Compound_12863']['symbolic_scaffold'] != candidates[RH + 'Compound_14399']['symbolic_scaffold']


def test_mismatched_contexts_and_unknown_carriers_cannot_pass():
    c = {RH + 'Compound_1': {'symbolic_scaffold': '[cytochrome c]'},
         RH + 'Compound_2': {'symbolic_scaffold': '[cytochrome cL]'}}
    e = {'id': 'test', 'left': [{'compound_id': RH + 'Compound_1', 'coefficient': 2}],
         'right': [{'compound_id': RH + 'Compound_2', 'coefficient': 2}], 'source_joins': []}
    assert audit_equation(e, c)['status'] == 'candidate-scaffold-nonconservation'
    assert audit_equation(e, {})['status'] == 'unresolved-carrier-scaffolds'
    c[RH + 'Compound_2']['symbolic_scaffold'] = '[cytochrome c]'
    assert audit_equation(e, c)['status'] == 'conditional-redox-scaffold-conservation'
    e['right'][0]['coefficient'] = 1
    assert audit_equation(e, c)['status'] == 'candidate-scaffold-nonconservation'
