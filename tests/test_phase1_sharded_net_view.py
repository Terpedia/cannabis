import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_sharded_net_view import write_view


def fixture():
    return {'schema':'fixture', 'summary':{'target_records':3}, 'claim_boundary':'hypotheses only',
        'targets':[{'cannabisdb_id':x, 'compound_id':'b', 'certificate_compound_id':'b',
                    'missing_candidate_reaction_ids':['r']} for x in ('one','alias')] + [
                    {'cannabisdb_id':'gap','compound_id':'c','certificate_compound_id':None,
                     'missing_candidate_reaction_ids':[]}],
        'certificates':[{'compound_id':'b','steps':[{'reaction_id':'r','extent':'1/2'}]}],
        'reactions':[{'id':'r','left':[{'compound_id':'a','coefficient':2}],
                      'right':[{'compound_id':'b','coefficient':1}], 'enzyme_evidence_ids':[]}],
        'compounds':[{'id':'a'},{'id':'b'},{'id':'c'}], 'enzyme_evidence':[]}


def checked(root, ref):
    data = (root/ref['file']).read_bytes()
    assert len(data)==ref['bytes']
    assert hashlib.sha256(data).hexdigest()==ref['sha256']
    return json.loads(data)


def test_shards_reconstruct_every_record_equation_and_certificate(tmp_path):
    bundle = fixture(); before = copy.deepcopy(bundle)
    stats = write_view(bundle,tmp_path)
    manifest = json.loads((tmp_path/'index.json').read_bytes())
    index = checked(tmp_path,manifest)
    shared = checked(tmp_path,index['shared_chemistry'])
    assert shared == {k:bundle[k] for k in ('reactions','compounds','enzyme_evidence')}
    assert stats['target_records']==3 and stats['certificates']==1
    assert index['summary']==bundle['summary']
    assert index['claim_boundary']==bundle['claim_boundary']
    assert [t['missing_candidate_reaction_count'] for t in index['targets']]==[1,1,0]
    assert all('missing_candidate_reaction_ids' not in t for t in index['targets'])
    assert checked(tmp_path,index['certificate_files']['b'])==bundle['certificates'][0]
    assert len(list((tmp_path/'certificates').iterdir()))==1
    assert bundle==before


@pytest.mark.parametrize('bad', ['participant','reaction','evidence','gap-count','identity','duplicate'])
def test_incomplete_or_inconsistent_shards_fail_closed(tmp_path,bad):
    bundle=fixture()
    if bad=='participant':bundle['compounds']=bundle['compounds'][1:]
    if bad=='reaction':bundle['reactions']=[]
    if bad=='evidence':bundle['reactions'][0]['enzyme_evidence_ids']=['missing']
    if bad=='gap-count':bundle['targets'][0]['missing_candidate_reaction_ids']=[]
    if bad=='identity':bundle['targets'][0]['compound_id']='a'
    if bad=='duplicate':bundle['certificates']*=2
    with pytest.raises(ValueError):write_view(bundle,tmp_path)
