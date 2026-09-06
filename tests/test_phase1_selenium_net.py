import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_selenium_net import proposals

ROOT=Path('data/reports')
def read(name):
    return json.loads((ROOT/('phase1-'+name+'.json')).read_bytes())


def test_proposals_preserve_exact_source_and_explicit_protons():
    parent=read('local-speciation-net'); source=read('selenium-source-identity')
    compounds={c['id']:c for c in parent['compounds']}
    reactions, merged=proposals(source,compounds)
    assert len(reactions)==3
    assert reactions[0]['left']==source['reaction']['left']
    assert reactions[0]['right']==source['reaction']['right']
    assert all(balanced([r['left'],r['right']],merged) for r in reactions)
    for r in reactions[1:]:
        assert {'compound_id':stable_id('structure','[H+]'),'coefficient':1} in r['right']
        assert 'computed' in r['source_evidence_type']
    assert set(compounds)<=set(merged)
    for cid in compounds:
        assert compounds[cid]==merged[cid]
    without_proton={k:v for k,v in compounds.items() if k!=stable_id('structure','[H+]')}
    with pytest.raises(ValueError,match='Missing exact'):
        proposals(source,without_proton)


@pytest.mark.parametrize('mode',['forward','reversible'])
def test_full_inventory_and_all_exact_certificates(mode):
    report=read('selenium-'+mode+'-net'); parent=read('local-speciation-net')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    layers=[json.loads((ROOT/n).read_bytes()) for n in report['baseline_certificate_reports']]
    reactions,compounds,model=assemble_current(read('full-balanced-network'),read('marts-completions'),report,layers)
    steps={s['id']:s for s in model.steps}
    inherited={c['compound_id']:c for doc in layers for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    assert len(inherited)==2720
    new={c['compound_id']:c for c in report['new_certificates']}
    assert not inherited.keys() & new.keys()
    assert report['external_exchange_compound_ids']==parent['external_exchange_compound_ids']
    assert report['co2_compound_id']==parent['co2_compound_id']
    reverse='selenium-source:KEGG-R03601:hypothetical-right-to-left'
    assert set(report['forbidden_step_ids'])==set(parent['forbidden_step_ids'])|({reverse} if mode=='forward' else set())
    assert (reverse in steps)==(mode=='reversible')
    assert 'selenium-source:KEGG-R03601:hypothetical-left-to-right' in steps
    for r in report['added_reactions'][1:]:
        assert all(r['id']+':hypothetical-'+d in steps for d in ('left-to-right','right-to-left'))
    for c in [*inherited.values(),*new.values()]:
        validate_certificate(c,steps,compounds,set(report['external_exchange_compound_ids']),report['co2_compound_id'])
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']]==[(t['cannabisdb_id'],t['compound_id']) for t in parent['targets']]
    assert report['summary']['target_records']==len(report['targets'])==6220
    assert report['summary']['balanced_equations']==len(reactions)==18177
    assert report['summary']['new_certificate_structures']==len(new)
    assert report['summary']['net_status_counts']==dict(Counter(t['net_status'] for t in report['targets']))
    for t,p in zip(report['targets'],parent['targets']):
        assert t['parent_net_status']==p['net_status']
        assert t['new_certificate']==(t['compound_id'] in new)
        assert (t['net_status']=='exact-net-conversion-hypothesis')==(t['compound_id'] in inherited or t['compound_id'] in new)
    for c in new.values():
        assert set(c['added_lipid_reaction_ids'])=={s['reaction_id'] for s in c['steps']} & {r['id'] for r in report['added_reactions']}
        assert 'selenium-source:KEGG-R03601' in c['added_lipid_reaction_ids']
