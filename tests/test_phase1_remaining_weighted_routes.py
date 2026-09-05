import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def test_remaining_diagnostic_covers_every_unresolved_target_and_replays_certificates():
    report = json.loads((ROOT / 'data/reports/phase1-remaining-weighted-routes.json').read_text())
    candidate = json.loads((ROOT / 'data/reports/phase1-thiolase-candidate-net.json').read_text())
    baseline = next(s for s in candidate['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    assert len(report['targets']) == 6220
    assert report['forbidden_step_ids'] == baseline['forbidden_step_ids']
    assert report['external_exchange_compound_ids'] == candidate['external_exchange_compound_ids']
    results = {r['compound_id']: r for r in report['results']}
    needed = set()
    for before, after in zip(baseline['targets'], report['targets']):
        assert all(after[k] == v for k, v in before.items())
        if before['net_status'] not in ('exact-net-conversion-hypothesis', 'explicit-exchange-species; not a synthesis target'):
            needed.add(before['compound_id'])
            assert after['diagnostic_status'] == results[before['compound_id']]['status']
        else:
            assert after['diagnostic_status'] == 'retained-baseline-status; not rescreened'
    assert set(results) == needed
    assert dict(Counter(t['diagnostic_status'] for t in report['targets'])) == report['summary']['target_diagnostic_status_counts']
    steps = {s['id']: s for s in orientations(report['reactions'])}
    compounds = {c['id']: c for c in report['compounds']}
    exchange = set(report['external_exchange_compound_ids'])
    assert {c for c in exchange if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    gaps = set()
    for result in results.values():
        selected = result.get('steps', [])
        missing = {s['reaction_id'] for s in selected} - candidate['candidate_reaction_evidence_ids'].keys()
        assert set(result['selected_missing_reaction_ids']) == missing
        gaps.update(missing)
        if result['status'] != 'exact-net-conversion-hypothesis':
            assert not selected
            continue
        assert not set(report['forbidden_step_ids']) & {s['step_id'] for s in selected}
        net = exact_net([steps[s['step_id']] for s in selected], [s['extent'] for s in selected])
        assert net[result['compound_id']] >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchange)
        assert {c: str(-n) for c, n in net.items() if n < 0} == result['external_net_consumption']
        assert {c: str(n) for c, n in net.items() if n > 0} == result['net_exports']
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
    assert gaps == {g['reaction_id'] for g in report['candidate_gaps']}
    sources = {}
    for path, sha in report['source_sha256'].items():
        payload = (ROOT / path).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == sha
        if path.endswith('search.json'):
            sources[path] = json.loads(payload)
    for gap in report['candidate_gaps']:
        assert gap['prior_searches'] == [{'report': p, 'row': row} for p, source in sorted(sources.items()) for row in source['rows'] if row['reaction_id'] == gap['reaction_id']]
