"""Paired original/name-derived PG probes; historical target inventory is unchanged."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_lipid_acylation_net import equation_key
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_reaction_completion_net import validate_certificate


def run():
    root=Path('data/reports')
    current=json.loads((root/'phase1-glycerophospholipid-net.json').read_bytes())
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',
        root/'phase1-glycerophospholipid-net.json',root/'phase1-pg-named-reactions.json']+[
            root/n for n in current['baseline_certificate_reports']]
    docs=[json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Changed source snapshot')
    network,completions,current,proposals,*layers=docs
    reactions,compounds,baseline=assemble_current(network,completions,current,layers)
    keys={equation_key(r):r['id'] for r in reactions.values()}
    added=[];joins=[];forbidden=set(current['forbidden_step_ids'])
    allowed={'glycerophospholipid-speciation','pgp-hydrolysis','pgp-synthesis','cdp-dag-synthesis','sn1-acylation','sn2-acylation'}
    for c in proposals['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles']!=c['smiles']:
            raise ValueError('Identity conflict')
        compounds.setdefault(c['id'],c)
    for r in proposals['reactions']:
        if r['hypothesis_type'] not in allowed or not balanced([r['left'],r['right']],compounds):
            raise ValueError('Unexpected or unbalanced proposal')
        key=equation_key(r)
        if key in keys:
            joins.append({'hypothesis_id':r['id'],'existing_reaction_id':keys[key]})
            continue
        reactions[r['id']]=r;keys[key]=r['id'];added.append(r)
        if r['hypothesis_type']!='glycerophospholipid-speciation':
            forbidden.add(r['id']+':hypothetical-right-to-left')
    model=NetModel(list(reactions.values()),current['external_exchange_compound_ids'],forbidden_step_ids=forbidden)
    steps={s['id']:s for s in model.steps}; baseline_steps={s['id']:s for s in baseline.steps}
    pairs=[]
    def solve(m,s,cid):
        result={'compound_id':cid,**m.solve(cid)}
        if result['status']=='exact-net-conversion-hypothesis':
            validate_certificate(result,s,compounds,set(current['external_exchange_compound_ids']),current['co2_compound_id'])
        return result
    for i,t in enumerate(proposals['identity_alternatives'],1):
        original=solve(model,steps,t['original_compound_id'])
        before=solve(baseline,baseline_steps,t['alternative_compound_id'])
        after=solve(model,steps,t['alternative_compound_id'])
        pairs.append({'cannabisdb_id':t['cannabisdb_id'],'source_name':t['source_name'],
            'original_result':original,'alternative_baseline_result':before,'alternative_extended_result':after,
            'identity_status':t['identity_status'],'source_url':t['source_url'],
            'occurrence_assertion':t['occurrence_assertion']})
        print(f'Named PG pair {i}/23: original={original["status"]}; alternative={after["status"]}',flush=True)
    used={s['reaction_id'] for p in pairs for k in ('original_result','alternative_baseline_result','alternative_extended_result') for s in p[k].get('steps',[])}
    report={'schema':'cannabis-carbon.phase1-pg-named-net.v1','paired_probes':pairs,
        'added_reactions':added,'existing_equation_joins':joins,'certificate_reactions':[reactions[r] for r in sorted(used)],
        'compounds':list(compounds.values()),'forbidden_step_ids':sorted(forbidden),
        'external_exchange_compound_ids':current['external_exchange_compound_ids'],'co2_compound_id':current['co2_compound_id'],
        'historical_inventory_summary_unchanged':current['summary'],
        'summary':{'paired_records':len(pairs),'balanced_equations':len(reactions),'added_equations':len(added),'existing_equation_joins':len(joins),
            **{k+'_counts':dict(Counter(p[k]['status'] for p in pairs)) for k in ('original_result','alternative_baseline_result','alternative_extended_result')},
            'changed_historical_identities':0,'historical_coverage_gain_claimed':0},
        'claim_boundary':'Separate 23-pair identity sensitivity, not a full-inventory coverage update or source correction. Names specify geometry omitted from original structures. No reaction links the two identity assertions. CO2 remains the sole external carbon input; inherited direction bounds are unchanged and new chemistry is source-forward except explicit reversible proton exchanges. Exact certificates permit regenerated pre-existing pools, not startup, energetics, compartments, detected occurrence or demonstrated Cannabis activity.',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'phase1-pg-named-net.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__=='__main__':
    run()
