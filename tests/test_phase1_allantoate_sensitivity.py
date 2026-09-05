import hashlib
import json
from fractions import Fraction
from pathlib import Path
from cannabis_carbon.phase1_allantoate_sensitivity import RID
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations


def test_all_affected_routes_replay_without_condensation():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / 'data/reports/phase1-allantoate-sensitivity.json').read_text())
    for path, digest in report['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    priority = json.loads((root / 'data/reports/phase1-current-gap-priority.json').read_text())
    gap = next(r for r in priority['rows'] if r['reaction_id'] == RID)
    assert {r['cannabisdb_id'] for r in report['rows']} == set(gap['remaining_target_ids'])
    assert len(report['rows']) == 14
    assert len(report['forbidden_step_ids']) == 9
    assert RID + ':hypothetical-left-to-right' in report['forbidden_step_ids']
    steps = {s['id']: s for s in orientations(report['reactions'])}
    exchanges = set(report['external_exchange_compound_ids'])
    for row in report['rows']:
        assert row['status'] == 'exact-net-conversion-hypothesis'
        assert not set(report['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
        net = exact_net([steps[s['step_id']] for s in row['steps']], [Fraction(s['extent']) for s in row['steps']])
        assert net[row['compound_id']] >= 1
        assert all(n >= 0 for cid, n in net.items() if cid not in exchanges)
        assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
