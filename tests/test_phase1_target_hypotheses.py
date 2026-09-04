import hashlib
import json
from pathlib import Path
import pytest

from cannabis_carbon.balance import _reaction_smiles_balance
from cannabis_carbon.phase1_target_coverage import audit_targets
from cannabis_carbon.phase1_target_hypotheses import build


def fixture():
    targets = [{'id': 'ethanol', 'label': 'ethanol', 'smiles': 'CCO'},
               {'id': 'methanol', 'label': 'methanol', 'smiles': 'CO'},
               {'id': 'spectator', 'label': 'sodium', 'smiles': '[Na+]'},
               {'id': 'gap', 'label': 'gap', 'smiles': 'CCCCC'}]
    sources = [{'id': f'test:{i}', 'source_reaction_id': f'TEST:{i}',
                'source_layer': 'test', 'reaction_smiles': smiles}
               for i, smiles in enumerate(['C=C.O>>CCO', 'CCO>>O.C=C',
                   'CO.C=C.O>>CO.CCO', 'CO.C=O.[H][H]>>CO.CO', 'C=C.O.[Na+]>>CCO.[Na+]'])]
    coverage = audit_targets(targets, {'reactions': []}, {'reactions': []}, sources, True)
    return coverage


def test_full_equation_dedup_net_production_and_no_silent_seeds():
    result = build(fixture(), {'reactions': []}, {'reactions': []})
    assert len(result['reactions']) == 4  # reverse encodings merge, spectator variants do not
    assert sum(len(r['sources']) for r in result['reactions']) == 5
    paired = next(r for r in result['reactions'] if len(r['sources']) == 2)
    assert {s['source_left_corresponds_to'] for s in paired['sources']} == {'left', 'right'}
    targets = {t['cannabisdb_id']: t for t in result['targets']}
    assert targets['spectator']['status'] == 'balanced-participation-only-no-net-production'
    assert targets['spectator']['hypothesis_ids'] == []
    assert targets['gap']['status'] == 'no-exact-encoded-reaction-match'
    assert targets['spectator']['carbon_count'] == 0
    assert 'uptake' in targets['spectator']['next_step']
    assert result['summary']['carbon_bearing_target_records'] == 3
    compounds = {c['id']: c for c in result['compounds']}
    assert all(t['compound_id'] in compounds for t in result['targets'])
    methanol = [h for h in result['hypotheses'] if h['cannabisdb_id'] == 'methanol']
    assert len(methanol) == 1
    assert methanol[0]['net_target_coefficient'] == 1
    assert any('requires-target-bootstrap' in b for b in methanol[0]['blockers'])
    assert all(h['status'] == 'blocked' for h in result['hypotheses'])
    assert all(i['availability'] == 'unestablished' for h in result['hypotheses'] for i in h['required_inputs'])
    ethanol = next(h for h in result['hypotheses'] if h['reaction_id'] == paired['id'])
    assert {compounds[i['compound_id']]['smiles'] for i in ethanol['required_inputs']} == {'C=C', 'O'}


def test_source_family_evidence_is_candidate_only():
    coverage = fixture()
    coverage['reaction_ledger'][0]['source_reaction_id'] = 'RHEA:2'
    core = {'id': 'rhea:1', 'directional_rhea_ids': ['2'], 'enzyme_ids': ['uniprot:example'],
            'enzyme_associations': [{'enzyme_id': 'uniprot:example', 'qualifiers': {'directExperimentalEvidence': False}}]}
    result = build(coverage, {'reactions': [core]}, {'reactions': []})
    evidence = result['enzyme_evidence'][0]
    assert evidence['enzyme_associations'] == core['enzyme_associations']
    assert any(h['has_candidate_enzyme_evidence'] for h in result['hypotheses'])
    assert all(h['status'] == 'blocked' and 'Cannabis-enzyme-activity-unconfirmed' in h['blockers'] for h in result['hypotheses'])


def test_rejects_false_balance_claim():
    coverage = fixture()
    coverage['reaction_ledger'][0]['reaction_smiles'] = 'CC>>C'
    with pytest.raises(ValueError, match='independent audit'):
        build(coverage, {'reactions': []}, {'reactions': []})


def test_published_hypotheses_full_inputs_net_coefficients_and_sources():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-target-hypotheses.json'
    report = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-target-hypotheses.json').read_bytes()
    for source, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    coverage = json.loads((root / 'data/reports/phase1-target-rhea-coverage.json').read_text())
    assert [t['cannabisdb_id'] for t in report['targets']] == [t['cannabisdb_id'] for t in coverage['targets']]
    compounds = {c['id']: c for c in report['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    sources = [s['coverage_record_id'] for r in reactions.values() for s in r['sources']]
    assert len(set(sources)) == len(sources)
    assert set(sources) == {r['id'] for r in coverage['reaction_ledger'] if r['computed_balance_status'] == 'balanced'}
    for reaction in reactions.values():
        equation = '>>'.join('.'.join(compounds[m['compound_id']]['smiles'] for m in reaction[side]
                                     for _ in range(m['coefficient'])) for side in ('left', 'right'))
        element, charge = _reaction_smiles_balance(equation)
        assert element['status'] == charge['status'] == 'balanced'
    hypotheses = {h['id']: h for h in report['hypotheses']}
    for target in report['targets']:
        assert target['compound_id'] in compounds
        for hid in target['hypothesis_ids']:
            assert hypotheses[hid]['cannabisdb_id'] == target['cannabisdb_id']
    for h in hypotheses.values():
        reaction = reactions[h['reaction_id']]
        input_side, output_side = ('left', 'right') if h['direction_mode'] == 'hypothetical-left-to-right' else ('right', 'left')
        assert [{k: i[k] for k in ('compound_id', 'coefficient')} for i in h['required_inputs']] == reaction[input_side]
        assert h['outputs'] == reaction[output_side]
        delta = sum(m['coefficient'] for m in h['outputs'] if m['compound_id'] == h['compound_id']) - sum(m['coefficient'] for m in h['required_inputs'] if m['compound_id'] == h['compound_id'])
        assert delta == h['net_target_coefficient'] > 0
        assert h['status'] == 'blocked'
