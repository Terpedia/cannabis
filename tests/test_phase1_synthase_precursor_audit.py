import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_all_exact_probes_and_certificates_preserve_boundaries():
    report = read('phase1-synthase-precursor-audit')
    model = read('phase1-synthase-candidate-net')
    network = read('phase1-full-balanced-network')
    compounds = {c['id']: c for c in network['compounds']}
    steps = {s['id']: s for s in orientations(network['reactions'])}
    probes = {p['id'] for p in report['probes']}
    assert len(probes) == 8 and len(report['results']) == 32
    exchanges = set(report['external_exchange_compound_ids'])
    assert {c for c in exchanges if compounds[c]['carbon_count']} == {report['co2_compound_id']}
    for scenario in report['scenarios']:
        rows = [r for r in report['results'] if r['scenario_id'] == scenario['id']]
        assert {r['compound_id'] for r in rows} == probes and len(rows) == 8
        for row in rows:
            if row['status'] != 'exact-net-conversion-hypothesis':
                continue
            assert not set(scenario['forbidden_step_ids']) & {s['step_id'] for s in row['steps']}
            if scenario['id'].startswith('candidate:'):
                assert {s['reaction_id'] for s in row['steps']} <= model['candidate_reaction_evidence_ids'].keys()
            net = exact_net([steps[s['step_id']] for s in row['steps']], [s['extent'] for s in row['steps']])
            assert net[row['compound_id']] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == row['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == row['net_exports']
            incoming = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
            outgoing = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
            assert incoming == outgoing > 0
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_selected_gaps_preserve_every_prior_search_row():
    report = read('phase1-synthase-precursor-audit')
    model = read('phase1-synthase-candidate-net')
    searches = {p: json.loads((ROOT / p).read_text()) for p in report['source_sha256'] if p.endswith('search.json')}
    used = {s['reaction_id'] for r in report['results'] for s in r.get('steps', [])}
    assert {g['reaction_id'] for g in report['candidate_gaps']} == used - model['candidate_reaction_evidence_ids'].keys()
    for gap in report['candidate_gaps']:
        expected = [{'report': path, 'row': row} for path, source in sorted(searches.items()) for row in source['rows'] if row['reaction_id'] == gap['reaction_id']]
        assert gap['prior_searches'] == expected
        assert gap['selected_uses'] == [{'probe_id': r['id'], **s} for r in report['results'] for s in r.get('steps', []) if s['reaction_id'] == gap['reaction_id']]
