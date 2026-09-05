import hashlib
import json
from fractions import Fraction
from pathlib import Path
from cannabis_carbon.phase1_allantoate_sensitivity import RID as ALLANTOATE
from cannabis_carbon.phase1_ureide_sensitivity import RID as UREIDOGLYCINE
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_allantoate_alternative_gaps import build


def test_ureide_routes_preserve_exact_net_and_all_carbon():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / 'data/reports/phase1-ureide-sensitivity.json').read_text())
    original = json.loads((root / 'data/reports/phase1-allantoate-sensitivity.json').read_text())
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    assert len(report['rows']) == 14
    assert {r['cannabisdb_id'] for r in report['rows']} == {r['cannabisdb_id'] for r in original['rows']}
    assert set(report['forbidden_step_ids']) == set(original['forbidden_step_ids']) | {UREIDOGLYCINE + ':hypothetical-left-to-right'}
    assert len(report['forbidden_step_ids']) == 10
    assert ALLANTOATE + ':hypothetical-left-to-right' in report['forbidden_step_ids']
    assert report['direction_review']['model_eligible'] is False
    compounds = {c['id']: c for c in report['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    assert exchanges == set(original['external_exchange_compound_ids'])
    assert [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] == ['O=C=O']
    steps = {s['id']: s for s in orientations(report['reactions'])}
    for row in report['rows']:
        assert row['status'] == 'exact-net-conversion-hypothesis'
        assert not set(report['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
        net = exact_net([steps[s['step_id']] for s in row['steps']], [Fraction(s['extent']) for s in row['steps']])
        assert net[row['compound_id']] == Fraction(row['target_amount']) >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == row['net_exports']
        carbon_in = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
        carbon_out = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
        assert carbon_in == carbon_out > 0


def test_new_gap_inventory_replays_every_used_missing_equation(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    path = 'data/reports/phase1-ureide-sensitivity.json'
    report = json.loads(Path('data/reports/phase1-ureide-alternative-gaps.json').read_text())
    rebuilt = build([p for p in report['source_sha256'] if p.endswith('search.json')], route_path=path)
    rebuilt['schema'] = 'cannabis-ureide-alternative-gaps-v1'
    assert rebuilt == report
    routes = json.loads(Path(path).read_text())
    model = json.loads(Path('data/reports/phase1-remaining-candidate-net.json').read_text())
    expected = {s['reaction_id'] for r in routes['rows'] for s in r['steps']} - model['candidate_reaction_evidence_ids'].keys()
    assert {g['reaction_id'] for g in report['candidate_gaps']} == expected
    for gap in report['candidate_gaps']:
        assert set(gap['target_ids']) == {r['cannabisdb_id'] for r in routes['rows'] if any(s['reaction_id'] == gap['reaction_id'] for s in r['steps'])}
