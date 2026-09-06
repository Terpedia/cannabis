import copy
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_fnsii_redox_hypothesis import APIGENIN, NARINGENIN, audit, build, describe


def test_exact_redox_hypothesis_preserves_protonation_and_explicit_assumptions(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-fnsii-redox-hypothesis.json').read_text())
    assert report == build()
    parts = report['participants']
    assert parts[NARINGENIN]['formula'] == 'C15H12O5'
    assert parts[APIGENIN]['formula'] == 'C15H9O5-'
    assert parts['NADPH']['formal_charge'] == -4
    assert parts['NADP']['formal_charge'] == -3
    assert parts['NADPH']['elements']['C'] == parts['NADP']['elements']['C'] == 21
    assert parts['RHEA-COMP:11964']['formal_charge'] == -2
    assert parts['RHEA-COMP:11965']['formal_charge'] == -3
    for step in report['steps']:
        assert step['balance'] == audit(step['stoichiometry'], parts)
        assert step['balance']['balanced'] is True
        assert step['model_eligible'] is False
    assert report['net']['balance'] == audit(report['net']['stoichiometry'], parts)
    assert report['net']['balance']['balanced'] is True
    assert set(report['canceled_participants']) == {'RHEA-COMP:11964', 'RHEA-COMP:11965', 'proton'}
    assert not {'RHEA-COMP:14627', 'RHEA-COMP:14628'} & parts.keys()
    assert report['compatible_fnsii_partners'] == report['external_organic_seeds'] == []
    assert report['model_eligible'] is report['summary']['candidate_model_changed'] is False
    assert report['summary']['new_exact_enzyme_assignments'] == 0


def test_audit_detects_generic_proton_count_and_missing_carbon_carrier(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = build()
    wrong_proton = copy.deepcopy(report['steps'][0]['stoichiometry'])
    wrong_proton['proton'] = 1
    result = audit(wrong_proton, report['participants'])
    assert not result['balanced']
    assert result['element_delta']['H'] == result['charge_delta'] == -1
    missing_nadp = copy.deepcopy(report['net']['stoichiometry'])
    del missing_nadp['NADP']
    result = audit(missing_nadp, report['participants'])
    assert not result['balanced'] and result['element_delta']['C'] == -21
    missing_water = copy.deepcopy(report['net']['stoichiometry'])
    missing_water['water'] = 1
    assert not audit(missing_water, report['participants'])['balanced']
    with pytest.raises(ValueError, match='exact'):
        describe('*CC')
