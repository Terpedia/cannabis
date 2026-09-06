import hashlib
import json
from pathlib import Path
from rdkit import Chem
import pytest
from cannabis_carbon.phase1_selenium_connectivity_audit import selenium_cut
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current


def step(name, inputs, outputs):
    return {'id':name, 'required_inputs':[{'compound_id':c,'coefficient':n} for c,n in inputs],
            'outputs':[{'compound_id':c,'coefficient':n} for c,n in outputs]}


def test_relaxed_connectivity_does_not_claim_all_input_supply_and_detects_closed_cycle():
    counts={'external':1,'a':1,'b':1,'target':2}
    closed=[step('ab',[('a',1)],[('b',1)]),step('ba',[('b',1)],[('a',1)])]
    cut=selenium_cut(closed,counts,{'external'})
    assert cut['unreachable_compound_ids']==['a','b','target']
    assert all(d['unreachable_selenium_weight_delta']=='0' for d in cut['step_weight_deltas'])
    relaxed=selenium_cut(closed+[step('mix',[('external',1),('a',1)],[('target',1)])],counts,{'external'})
    assert 'target' in relaxed['relaxed_reachable_compound_ids']
    assert 'a' in relaxed['unreachable_compound_ids']  # Connectivity is deliberately only an upper bound.
    with pytest.raises(ValueError,match='without selenium input'):
        selenium_cut([step('birth',[],[('a',1)])],counts,{'external'})
    with pytest.raises(ValueError,match='increases its weight'):
        selenium_cut([step('unbalanced',[('a',1)],[('b',2)])],counts,{'external'})


def test_selenium_cut_replays_all_allowed_equations_and_every_target():
    root=Path('data/reports');read=lambda p:json.loads(p.read_bytes())
    report=read(root/'phase1-selenium-connectivity-audit.json')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    current=read(root/'phase1-local-speciation-net.json')
    rs,cs,model=assemble_current(read(root/'phase1-full-balanced-network.json'),read(root/'phase1-marts-completions.json'),
                                current,[read(root/n) for n in current['baseline_certificate_reports']])
    counts={cid:sum(a.GetAtomicNum()==34 for a in Chem.MolFromSmiles(c['smiles']).GetAtoms()) for cid,c in cs.items()}
    replay=selenium_cut(model.steps,counts,current['external_exchange_compound_ids'])
    assert all(report[k]==v for k,v in replay.items())
    assert report['summary']['whole_model_equations']==len(rs)==18174
    assert report['summary']['whole_model_directed_steps']==len(model.steps)
    assert report['summary']['coverage_gain_claimed']==0
    assert report['external_exchange_compound_ids']==current['external_exchange_compound_ids']
    assert report['forbidden_step_ids']==current['forbidden_step_ids']
    unreachable=set(report['unreachable_compound_ids'])
    assert unreachable.isdisjoint(current['external_exchange_compound_ids'])
    targets={t['cannabisdb_id']:t for t in current['targets'] if counts[t['compound_id']]}
    assert len(targets)==len(report['targets'])==3
    for t in report['targets']:
        assert all(t[k]==v for k,v in targets[t['cannabisdb_id']].items())
        assert t['compound_id'] in unreachable and t['selenium_count']==counts[t['compound_id']]>0
        assert t['net_status']=='solver-reported-infeasible'
    se={cid for cid,n in counts.items() if n}
    expected={r['id']:r for r in rs.values() if any(p['compound_id'] in se for side in ('left','right') for p in r[side])}
    assert {r['id']:r for r in report['reactions']}==expected
