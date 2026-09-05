import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_net_flux import NetModel
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_direction_sensitivity import build

ROOT = Path(__file__).resolve().parents[1]


def reaction(rid, a, b):
    return {'id':rid,'left':[{'compound_id':a,'coefficient':1}], 'right':[{'compound_id':b,'coefficient':1}]}


def test_direction_restriction_removes_only_requested_orientation():
    r = reaction('r','a','b')
    forward = 'r:hypothetical-left-to-right'
    model = NetModel([r],{'a'},forbidden_step_ids=[forward])
    assert [s['direction_mode'] for s in model.steps] == ['hypothetical-right-to-left']
    assert model.solve('b')['status'] == 'no-net-producing-candidate-equation'
    assert NetModel([r],{'a'}).solve('b')['status'] == 'exact-net-conversion-hypothesis'
    alternate = NetModel([r,reaction('s','a','c'),reaction('t','c','b')],{'a'},forbidden_step_ids=[forward])
    result = alternate.solve('b')
    assert result['status'] == 'exact-net-conversion-hypothesis'
    assert {s['reaction_id'] for s in result['steps']} == {'s','t'}
    with pytest.raises(ValueError,match='Unknown forbidden'):
        NetModel([r],{'a'},forbidden_step_ids=['typo'])


def test_full_sensitivity_inventory_preserved_witnesses_and_replay(monkeypatch):
    names = ['phase1-full-balanced-network','phase1-expanded-candidate-net','phase1-candidate-direction-review','phase1-direction-sensitivity']
    network,expanded,review,report = [json.loads((ROOT/'data/reports'/(n+'.json')).read_text()) for n in names]
    before = copy.deepcopy((network,expanded,review))
    forbidden = {r['id'] for r in report['constraints']}
    admitted = [r for r in network['reactions'] if r['id'] in expanded['candidate_reaction_evidence_ids']]
    all_steps = {s['id']:s for s in orientations(admitted)}
    restricted = NetModel(admitted,expanded['external_exchange_compound_ids'],forbidden_step_ids=forbidden)
    assert {s['id'] for s in restricted.steps} == all_steps.keys()-forbidden
    assert len(forbidden)==5 and len(restricted.steps)==3171
    assert report['summary']['unique_allowed_reactant_compounds'] == len({p['compound_id'] for s in restricted.steps for p in s['required_inputs']}) == 1667
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']] == [(t['cannabisdb_id'],t['compound_id']) for t in expanded['targets']]
    assert len(report['targets'])==6220
    old = {c['compound_id']:c for c in expanded['certificates']}
    expected = {cid for cid,c in old.items() if not forbidden & {s['step_id'] for s in c['steps']}}
    assert set(report['preserved_certificate_compound_ids'])==expected
    assert len(expected)==100
    assert not report['alternative_certificates']
    assert {t['restricted_net_status'] for t in report['targets'] if t['was_new_expanded_target']}=={'solver-reported-infeasible'}
    for c in report['constraints']:
        r=next(r for r in review['reviews'] if r['id']==c['review_id'])
        allowed='hypothetical-left-to-right' if r['source_left_corresponds_to']=='left' else 'hypothetical-right-to-left'
        assert c['forbidden_direction_mode'] != allowed
        assert r['reaction_id']+':'+allowed in {s['id'] for s in restricted.steps}
    for path,sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    # Assembly replay avoids depending on identical optimizer paths across
    # SciPy versions. Tiny network above independently tests real LP restriction.
    outcomes={t['compound_id']:{'status':t['restricted_net_status']} | {k:t[k] for k in ('solver_status','solver_message') if k in t} for t in report['targets']}
    monkeypatch.setattr('cannabis_carbon.phase1_direction_sensitivity.NetModel.solve',lambda self,cid:outcomes[cid])
    rebuilt=build(network,expanded,review)
    assert {k:v for k,v in rebuilt.items() if k!='scipy_version'}=={k:v for k,v in report.items() if k not in ('source_sha256','scipy_version')}
    assert (network,expanded,review)==before
