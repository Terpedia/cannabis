import pytest
import json
import hashlib
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_net_flux import NetModel
from cannabis_carbon.phase1_net_obstruction import solve,validate


def reaction(rid,left,right):
    return {'id':rid,'left':[{'compound_id':c,'coefficient':1} for c in left],
            'right':[{'compound_id':c,'coefficient':1} for c in right]}


def test_exact_obstruction_does_not_confuse_pool_recycling_with_net_production():
    model=NetModel([reaction('r',['precursor'],['target'])],{'CO2'})
    result=solve(model,'target')
    assert result['status']=='exact-stoichiometric-obstruction'
    assert result['weights']=={'precursor':'1','target':'1'}
    assert len(validate(model,'target',result['weights']))==2
    with pytest.raises(ValueError,match='increases'):
        validate(model,'target',{'target':'1'})
    with pytest.raises(ValueError,match='weights'):
        validate(model,'target',{'target':'1','CO2':'1'})
    with pytest.raises(ValueError,match='weights'):
        validate(model,'target',{'target':'-1'})


def test_available_input_and_regenerated_pool_preclude_false_obstruction():
    model=NetModel([reaction('r',['CO2','carrier'],['target','spent']),reaction('regen',['spent'],['carrier'])],{'CO2'})
    assert model.solve('target')['status']=='exact-net-conversion-hypothesis'
    assert solve(model,'target')['status']=='no-obstruction-certificate-produced'
    with pytest.raises(ValueError,match='increases'):
        validate(model,'target',{'target':'1','carrier':'1','spent':'1'})


@pytest.mark.parametrize('scenario,equations', [('alkane',18203), ('odd-chain',18211)])
def test_full_alkane_model_obstructions_replay_every_allowed_direction(scenario,equations):
    root=Path('data/reports'); read=lambda n:json.loads((root/n).read_bytes())
    report=read('phase1-'+scenario+'-obstructions.json'); current=read('phase1-'+scenario+'-net.json')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    layers=[read(n) for n in current['baseline_certificate_reports']]
    reactions,compounds,model=assemble_current(read('phase1-full-balanced-network.json'),read('phase1-marts-completions.json'),current,layers)
    assert report['forbidden_step_ids']==current['forbidden_step_ids']
    assert report['external_exchange_compound_ids']==current['external_exchange_compound_ids']
    assert report['summary']['balanced_equations']==len(reactions)==equations
    assert report['summary']['allowed_steps']==len(model.steps)
    assert {t['cannabisdb_id'] for t in report['targets']}=={'CDB000155','CDB000157'}
    exact=[r for r in report['targets'] if r['status']=='exact-stoichiometric-obstruction']
    assert report['summary']['exact_obstructions']==len(exact)==2
    for r in exact:
        assert r['weighted_steps']==validate(model,r['compound_id'],r['weights'])
        assert r['checked_allowed_steps']==len(model.steps)
    assert {c['id']:c for c in report['compounds']}=={cid:compounds[cid] for r in exact for cid in r['weights']}
