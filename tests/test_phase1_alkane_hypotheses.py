import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_alkane_hypotheses import build,REFERENCE
from cannabis_carbon.phase1_marts_completions import balanced


def test_exact_chain_window_and_complete_deformylation_bookkeeping():
    read=lambda name:json.loads(Path('data/reports/phase1-'+name+'.json').read_bytes())
    report=read('alkane-hypotheses'); current=read('ketone-stereo-net'); network=read('full-balanced-network')
    reference=next(r for r in network['reactions'] if r['id']==REFERENCE)
    assert {k:v for k,v in report.items() if k!='source_sha256'}==build(current,reference)
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    assert report['summary']=={'full_inventory_records':6220,'straight_chain_no_producer_records':25,'balanced_proposals':4,'coverage_gain_claimed':0}
    cs={c['id']:c for c in report['compounds']}; reactions={r['id']:r for r in report['reactions']}
    for t in report['targets']:
        assert t['aldehyde_carbon_count']==t['alkane_carbon_count']+1
        if 'hypothesis_id' not in t:
            assert not 13<=t['aldehyde_carbon_count']<=22
            continue
        r=reactions[t['hypothesis_id']]
        assert not t['aldehyde_present_in_parent_model']
        assert balanced([r['left'],r['right']],cs)
        assert len(r['left'])==len(r['right'])==4
        assert not r['enzyme_evidence_ids']
        assert any(cs[p['compound_id']]['smiles']=='O=C[O-]' for p in r['right'])
        assert any(cs[p['compound_id']]['smiles']=='O=O' for p in r['left'])
        assert sorted(p['coefficient'] for p in r['left'])==[1,1,1,2]
        assert sorted(p['coefficient'] for p in r['right'])==[1,1,1,2]
