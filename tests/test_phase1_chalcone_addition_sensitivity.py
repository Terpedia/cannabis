import hashlib
import json
from fractions import Fraction
from pathlib import Path
from cannabis_carbon.phase1_chalcone_reference_search import RID
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations


def test_chalcone_counterfactual_uses_complete_baseline_and_exact_net_certificates(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-addition-sensitivity.json').read_text())
    parent = json.loads(Path('data/reports/phase1-remaining-candidate-net.json').read_text())
    priority = json.loads(Path('data/reports/phase1-current-gap-priority.json').read_text())
    gap = next(r for r in priority['rows'] if r['reaction_id'] == RID)
    old = next(s for s in parent['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert report['model_eligible'] is False
    assert {r['id'] for r in report['reactions']} == set(parent['candidate_reaction_evidence_ids']) | {RID}
    baseline, addition = report['scenarios']
    assert baseline['reaction_count'] == 1609 and addition['reaction_count'] == 1610
    assert baseline['forbidden_step_ids'] == old['forbidden_step_ids']
    assert addition['forbidden_step_ids'] == old['forbidden_step_ids'] + [RID + ':hypothetical-right-to-left']
    assert report['external_exchange_compound_ids'] == parent['external_exchange_compound_ids']
    compounds = {c['id']: c for c in report['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    assert [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] == ['O=C=O']
    steps = {s['id']: s for s in orientations(report['reactions'])}
    for scenario in report['scenarios']:
        assert {r['cannabisdb_id'] for r in scenario['rows']} == set(gap['remaining_target_ids'])
        for row in scenario['rows']:
            if row['status'] != 'exact-net-conversion-hypothesis':
                continue
            assert not set(scenario['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
            assert RID in {s['reaction_id'] for s in row['steps']}
            net = exact_net([steps[s['step_id']] for s in row['steps']], [s['extent'] for s in row['steps']])
            assert net[row['compound_id']] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == row['net_exports']
            assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    assert baseline['summary'] == {'solver-reported-infeasible': 7}
    assert addition['summary'] == {'exact-net-conversion-hypothesis': 5, 'solver-reported-infeasible': 2}
