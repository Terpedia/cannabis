import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import input_closure
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_closure_keeps_all_inputs_and_handles_cycles_without_claiming_supply():
    rows = [{'compound_id':'a','reaction_id':'r1','required_precursor_ids':['b','c']},
            {'compound_id':'b','reaction_id':'r2','required_precursor_ids':['a','d']}]
    assert input_closure('a',rows) == (['c','d'],['r1','r2'])


def test_full_frontier_audit_preserves_boundaries_and_replays_certificates():
    r = json.loads(Path('data/reports/phase1-glycerophospholipid-input-audit.json').read_bytes())
    current = json.loads(Path('data/reports/phase1-glycerophospholipid-net.json').read_bytes())
    proposals = json.loads(Path('data/reports/phase1-glycerophospholipid-synthesis.json').read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    for k in ('external_exchange_compound_ids','co2_compound_id','forbidden_step_ids'):
        assert r[k] == current[k]
    assert r['summary']['full_model_balanced_equations'] == current['summary']['balanced_equations'] == 18082
    assert r['summary']['new_CO2_route_claims'] == 0
    assert {i['compound_id'] for i in r['input_results']} == {i['compound_id'] for i in proposals['frontier']}
    assert len(r['input_results']) == 21
    assert len(r['targets']) == 157
    cs = {c['id']:c for c in r['compounds']}
    forbidden = set(r['forbidden_step_ids'])
    steps = {s['id']:s for s in orientations(r['certificate_reactions']) if s['id'] not in forbidden}
    for reaction in r['certificate_reactions']:
        assert balanced([reaction['left'],reaction['right']],cs)
    for result in r['input_results']:
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate(result,steps,cs,set(r['external_exchange_compound_ids']),r['co2_compound_id'])
    source_targets = {t['cannabisdb_id']:t for t in proposals['targets']}
    statuses = {t['cannabisdb_id']:t['net_status'] for t in current['targets']}
    results = {i['compound_id']:i for i in r['input_results']}
    for t in r['targets']:
        original = source_targets[t['cannabisdb_id']]
        assert t['compound_id'] == original['compound_id']
        assert t['net_status'] == statuses[t['cannabisdb_id']]
        assert (t['terminal_input_ids'],t['template_reaction_ids']) == input_closure(original['source_form_compound_id'],proposals['precursor_candidates'])
        assert t['unresolved_terminal_input_ids'] == [c for c in t['terminal_input_ids'] if results[c]['status'] not in (
            'exact-net-conversion-hypothesis','explicit-exchange-species; not a synthesis target')]
