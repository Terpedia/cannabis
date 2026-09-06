import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_lipid_acylation_net import equation_key
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_full_inventory_forward_only_redox_hypotheses_and_preserved_baseline():
    root=Path('data/reports'); read=lambda n:json.loads((root/('phase1-'+n+'.json')).read_bytes())
    report=read('ketone-stereo-net'); parent=read('selenium-forward-net'); proposals=read('ketone-stereo-hypotheses')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    layers=[json.loads((root/n).read_bytes()) for n in report['baseline_certificate_reports']]
    reactions,compounds,model=assemble_current(read('full-balanced-network'),read('marts-completions'),report,layers)
    steps={s['id']:s for s in model.steps}
    prior={c['compound_id']:c for layer in layers for c in layer.get('certificates',[])+layer.get('new_certificates',[])}
    assert len(prior)==2721
    new={c['compound_id']:c for c in report['new_certificates']}
    assert not prior.keys() & new.keys()
    for c in [*prior.values(),*new.values()]:
        validate_certificate(c,steps,compounds,set(parent['external_exchange_compound_ids']),parent['co2_compound_id'])
    added={r['id']:r for r in report['added_reactions']}; proposed={r['id']:r for r in proposals['reactions']}
    assert set(added)|{j['hypothesis_id'] for j in report['existing_equation_joins']}==set(proposed)
    for rid,r in added.items():
        assert r==proposed[rid]
        assert rid+':hypothetical-left-to-right' in steps
        assert rid+':hypothetical-right-to-left' not in steps
    for join in report['existing_equation_joins']:
        assert equation_key(proposed[join['hypothesis_id']])==equation_key(reactions[join['existing_reaction_id']])
    assert set(report['forbidden_step_ids'])==set(parent['forbidden_step_ids'])|{rid+':hypothetical-right-to-left' for rid in added}
    assert report['external_exchange_compound_ids']==parent['external_exchange_compound_ids']
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']]==[(t['cannabisdb_id'],t['compound_id']) for t in parent['targets']]
    assert len(report['targets'])==6220
    assert report['summary']['balanced_equations']==18177+len(added)
    assert report['summary']['net_status_counts']==dict(Counter(t['net_status'] for t in report['targets']))
    for t in report['targets']:
        assert (t['net_status']=='exact-net-conversion-hypothesis')==(t['compound_id'] in prior or t['compound_id'] in new)
        assert t['new_certificate']==(t['compound_id'] in new)
    for c in new.values():
        assert set(c['added_lipid_reaction_ids'])=={s['reaction_id'] for s in c['steps']} & added.keys()
        assert c['added_lipid_reaction_ids']
