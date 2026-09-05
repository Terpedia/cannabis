import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations


def test_remaining_flavonoid_witnesses_retain_all_cofactors_and_gap_directions(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-remaining-gaps.json').read_text())
    parent = json.loads(Path('data/reports/phase1-chalcone-addition-sensitivity.json').read_text())['scenarios'][1]
    candidate = json.loads(Path('data/reports/phase1-remaining-candidate-net.json').read_text())
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert report['model_eligible'] is False
    assert {r['cannabisdb_id'] for r in report['rows']} == {r['cannabisdb_id'] for r in parent['rows'] if r['status'] != 'exact-net-conversion-hypothesis'} == {'CDB005071', 'CDB005072'}
    assert report['forbidden_step_ids'] == parent['forbidden_step_ids']
    steps = {s['id']: s for s in orientations(report['reactions'])}
    compounds = {c['id']: c for c in report['compounds']}
    exchanges = set(report['external_exchange_compound_ids'])
    assert [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] == ['O=C=O']
    missing = set()
    for row in report['rows']:
        assert row['status'] == 'exact-net-conversion-hypothesis'
        assert not set(report['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
        net = exact_net([steps[s['step_id']] for s in row['steps']], [s['extent'] for s in row['steps']])
        assert net[row['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchanges)
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
        assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
        missing.update(s['reaction_id'] for s in row['steps'] if s['reaction_id'] not in candidate['candidate_reaction_evidence_ids'] and s['reaction_id'] != report['hypothetical_discounted_reaction_id'])
    assert missing == {g['reaction_id'] for g in report['candidate_gaps']}
    assert len(missing) == 1
    gap = report['candidate_gaps'][0]
    assert len(gap['reaction']['left']) == 5 and len(gap['reaction']['right']) == 3
    assert {u['direction_mode'] for u in gap['selected_uses']} == {'hypothetical-right-to-left'}
