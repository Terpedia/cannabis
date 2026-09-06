import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_fnsii_addition_sensitivity import FNSII, CHI, added_reaction
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations


def test_factorial_sensitivity_retains_full_baseline_and_exact_carbon_accounting(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-fnsii-addition-sensitivity.json').read_text())
    parent = json.loads(Path('data/reports/phase1-remaining-candidate-net.json').read_text())
    redox = json.loads(Path('data/reports/phase1-fnsii-redox-hypothesis.json').read_text())
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert report['model_eligible'] is False
    assert report['baseline_candidate_reaction_evidence_ids'] == parent['candidate_reaction_evidence_ids']
    baseline_ids = set(parent['candidate_reaction_evidence_ids'])
    assert len(baseline_ids) == 1609
    assert {r['id'] for r in report['reactions']} == baseline_ids | {CHI, FNSII}
    fns = next(r for r in report['reactions'] if r['id'] == FNSII)
    assert fns == added_reaction(redox, report['compounds'])
    assert fns['model_eligible'] is False and fns['enzyme_evidence_ids'] == []
    assert report['external_exchange_compound_ids'] == parent['external_exchange_compound_ids']
    compounds = {c['id']: c for c in report['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    carbon_exchanges = [c for c in exchanges if compounds[c]['carbon_count']]
    assert len(carbon_exchanges) == 1
    assert compounds[carbon_exchanges[0]]['smiles'] == 'O=C=O'
    baseline = next(s for s in parent['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    steps = {s['id']: s for s in orientations(report['reactions'])}
    assert [s['reaction_count'] for s in report['scenarios']] == [1609, 1610, 1610, 1611]
    for scenario in report['scenarios']:
        assert scenario['forbidden_step_ids'] == baseline['forbidden_step_ids'] + [r + ':hypothetical-right-to-left' for r in scenario['added_reaction_ids']]
        assert {r['cannabisdb_id'] for r in scenario['rows']} == {'CDB005071', 'CDB005072'}
        for row in scenario['rows']:
            if scenario['id'] != 'CHI-and-FNSII':
                assert row['status'] == 'solver-reported-infeasible'
                continue
            assert row['status'] == 'exact-net-conversion-hypothesis'
            used = {s['reaction_id'] for s in row['steps']}
            assert {CHI, FNSII} <= used <= baseline_ids | {CHI, FNSII}
            assert not set(scenario['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
            net = exact_net([steps[s['step_id']] for s in row['steps']], [s['extent'] for s in row['steps']])
            assert net[row['compound_id']] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == row['net_exports']
            carbon_in = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
            carbon_out = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
            assert carbon_in == carbon_out > 0
            assert -net[carbon_exchanges[0]] == carbon_in
            assert net[fns['source_participant_resolution']['NADPH']] >= 0


def test_exact_structure_resolution_rejects_missing_or_duplicate_matches(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    redox = json.loads(Path('data/reports/phase1-fnsii-redox-hypothesis.json').read_text())
    with pytest.raises(ValueError, match='unique structure'):
        added_reaction(redox, [])
    compounds = json.loads(Path('data/reports/phase1-full-balanced-network.json').read_text())['compounds']
    fns = added_reaction(redox, compounds)
    c = next(c for c in compounds if c['id'] == fns['source_participant_resolution']['NADPH'])
    with pytest.raises(ValueError, match='unique structure'):
        added_reaction(redox, compounds + [c])
