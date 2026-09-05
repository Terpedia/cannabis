import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_net_flux import exact_net
from cannabis_carbon.phase1_scope import orientations

ROOT = Path(__file__).resolve().parents[1]


def test_every_cbga_input_and_net_witness_is_preserved():
    r = json.loads((ROOT / 'data/reports/phase1-cbga-input-audit.json').read_text())
    model = json.loads((ROOT / 'data/reports/phase1-synthase-candidate-net.json').read_text())
    probes = {m['compound_id'] for reaction in r['selected_reactions'] for side in ('left', 'right') for m in reaction[side]}
    assert probes == set(r['probe_compound_ids']) and len(probes) == 10
    assert len(r['selected_reactions']) == 3 and len(r['results']) == 20
    for reaction in r['selected_reactions']:
        assert reaction['candidate_evidence_ids'] == model['candidate_reaction_evidence_ids'][reaction['id']]
    compounds = {c['id']: c for c in r['compounds']}
    steps = {s['id']: s for s in orientations(r['reactions'])}
    exchanges = set(r['external_exchange_compound_ids'])
    assert {c for c in exchanges if compounds[c]['carbon_count']} == {r['co2_compound_id']}
    for sid, forbidden in r['scenario_constraints'].items():
        results = [x for x in r['results'] if x['scenario_id'] == sid]
        assert {x['compound_id'] for x in results} == probes and len(results) == 10
        successful = [x for x in results if x['status'] == 'exact-net-conversion-hypothesis']
        assert len(successful) == 1
        assert compounds[successful[0]['compound_id']]['carbon_count'] == 10
        for result in successful:
            assert not set(forbidden) & {s['step_id'] for s in result['steps']}
            net = exact_net([steps[s['step_id']] for s in result['steps']], [s['extent'] for s in result['steps']])
            assert net[result['compound_id']] >= 1
            assert all(n >= 0 for c, n in net.items() if c not in exchanges)
            assert {c: str(-n) for c, n in net.items() if n < 0} == result['external_net_consumption']
            assert {c: str(n) for c, n in net.items() if n > 0} == result['net_exports']
            incoming = sum(-n * compounds[c]['carbon_count'] for c, n in net.items() if n < 0)
            outgoing = sum(n * compounds[c]['carbon_count'] for c, n in net.items() if n > 0)
            assert incoming == outgoing > 0
    for p, sha in r['source_sha256'].items():
        assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == sha
