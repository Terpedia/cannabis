import pytest
import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_medium_boundary import UptakeLimitedModel
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_medium_boundary import blocked_species
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def reaction(id,left,right):
    return {'id':id,'left':[{'compound_id':c,'coefficient':1} for c in left],
            'right':[{'compound_id':c,'coefficient':1} for c in right]}


def test_blocked_uptake_does_not_block_secretion():
    r=reaction('r',['carbon'],['target','waste'])
    m=UptakeLimitedModel([r],{'carbon','waste'},{'carbon'})
    c=m.solve('target')
    assert c['status']=='exact-net-conversion-hypothesis'
    assert c['external_net_consumption']=={'carbon':'1'}
    assert c['net_exports']=={'target':'1','waste':'1'}
    assert m.solve('waste')=={'status':'external-inventory-species-not-assessed-for-synthesis','uptake_allowed':False}


def test_missing_input_cannot_be_imported_but_can_be_regenerated():
    r=reaction('r',['carbon','carrier'],['target','spent'])
    m=UptakeLimitedModel([r],{'carbon','carrier','spent'},{'carbon'})
    assert m.solve('target')['status']=='solver-reported-infeasible'
    regeneration=reaction('regen',['spent'],['carrier'])
    c=UptakeLimitedModel([r,regeneration],{'carbon','carrier','spent'},{'carbon'}).solve('target')
    assert c['status']=='exact-net-conversion-hypothesis'
    assert c['external_net_consumption']=={'carbon':'1'}
    assert not c['zero_net_internal_participants']
    with pytest.raises(ValueError,match='subset'):
        UptakeLimitedModel([r],{'carbon'},{'carbon','unknown'})


def test_complete_medium_restriction_report_and_exact_witnesses():
    root=Path('data/reports'); read=lambda n:json.loads((root/n).read_bytes())
    report=read('phase1-medium-boundary.json'); parent=read('phase1-selenium-forward-net.json')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    layers=[read(n) for n in parent['baseline_certificate_reports']]
    reactions,compounds,_=assemble_current(read('phase1-full-balanced-network.json'),read('phase1-marts-completions.json'),parent,layers)
    exchange=set(parent['external_exchange_compound_ids'])
    assert report['blocked_inputs']==blocked_species(compounds,exchange)
    blocked={r['compound_id'] for r in report['blocked_inputs']}
    uptake=exchange-blocked
    assert set(report['allowed_uptake_compound_ids'])==uptake
    assert set(report['allowed_external_output_compound_ids'])==exchange
    assert report['forbidden_step_ids']==parent['forbidden_step_ids']
    model=UptakeLimitedModel(list(reactions.values()),exchange,uptake,report['forbidden_step_ids'])
    steps={s['id']:s for s in model.steps}
    certs={c['compound_id']:c for c in report['certificates']}
    assert len(certs)==len(report['certificates'])==report['summary']['certificate_structures']
    for c in certs.values():
        validate_certificate(c,steps,compounds,uptake,report['co2_compound_id'])
        assert not set(c['external_net_consumption']) & blocked
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']]==[(t['cannabisdb_id'],t['compound_id']) for t in parent['targets']]
    assert report['summary']['target_records']==len(report['targets'])==6220
    assert report['summary']['balanced_equations']==len(reactions)==18177
    assert report['summary']['net_status_counts']==dict(Counter(t['net_status'] for t in report['targets']))
    for t,p in zip(report['targets'],parent['targets']):
        assert t['parent_net_status']==p['net_status']
        assert (t['net_status']=='exact-net-conversion-hypothesis')==(t['compound_id'] in certs)
        if t['compound_id'] in exchange:
            assert t['net_status']=='external-inventory-species-not-assessed-for-synthesis'
            assert t['uptake_allowed']==(t['compound_id'] in uptake)
    assert report['summary']['lost_covered_records']==sum(t['parent_net_status']=='exact-net-conversion-hypothesis' and t['net_status']!='exact-net-conversion-hypothesis' for t in report['targets'])
