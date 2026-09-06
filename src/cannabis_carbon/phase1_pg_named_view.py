"""Static paired-identity view; separate 23-record denominator and original assertions."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_chemistry_route_view import reaction_sources
from .phase1_net_view import build as attach_evidence


def run():
    paths=[Path('data/reports/phase1-pg-named-net.json'),Path('data/reports/phase1-pg-named-reactions.json'),
        Path('docs/data/glycerophospholipid-net-view/bundle.json')]+[Path('data/reports/'+n+'.json') for n in
        ('phase1-target-hypotheses','phase1-screened-enzyme-overlay','phase1-route-enzyme-overlay')]
    net,proposals,previous,*evidence=[json.loads(p.read_bytes()) for p in paths]
    for doc in (net,proposals,previous):
        for p,sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Changed view source')
    prior={r['id']:r for r in previous['reactions']}
    sources={r['rule_id']:r for r in proposals['source_records']}
    reactions=[]
    for r in net['certificate_reactions']:
        if r['id'] in prior:
            old=prior[r['id']]
            if any(old[s]!=r[s] for s in ('left','right')):
                raise ValueError('Prior equation identity conflict')
            item=dict(old)
        else:
            hypothesis=bool(r.get('hypothesis_type'))
            if hypothesis and r.get('source_reaction_id') not in sources and not r.get('sources'):
                raise ValueError('Missing inherited hypothesis source snapshot')
            item={**r,'enzyme_evidence_ids':r.get('enzyme_evidence_ids',[]),
                'sources':reaction_sources(r,sources) if hypothesis else r.get('sources',[])}
        if r.get('hypothesis_type'):
            item.update({'is_route_sensitivity':True,'hypothesis_assumptions':[r['claim_boundary'],net['claim_boundary']]})
        if not item.get('sources'):
            raise ValueError('Missing reaction display provenance: '+r['id'])
        item['missing_candidate_evidence']=not bool(item['enzyme_evidence_ids'])
        reactions.append(item)
    byid={r['id']:r for r in reactions}
    options=[]
    for key,title in [('original_result','Original encoded structures, extended chemistry'),
                      ('alternative_baseline_result','Name-derived alternatives, baseline chemistry'),
                      ('alternative_extended_result','Name-derived alternatives, extended chemistry')]:
        targets=[];certs=[]
        for pair in net['paired_probes']:
            result=pair[key];good=result['status']=='exact-net-conversion-hypothesis'
            if good:certs.append(result)
            targets.append({'cannabisdb_id':pair['cannabisdb_id'],'label':pair['source_name']+' · '+title,
                'compound_id':result['compound_id'],'certificate_compound_id':result['compound_id'] if good else None,
                'net_status':result['status'],'startup_status':'not established; identity alternatives are separate',
                'identity_status':'original-encoded-structure; name-conflict-unresolved' if key=='original_result' else pair['identity_status'],'source_url':pair['source_url'],
                'occurrence_assertion':pair['occurrence_assertion'],
                'missing_candidate_reaction_ids':sorted({s['reaction_id'] for s in result.get('steps',[]) if not byid[s['reaction_id']]['enzyme_evidence_ids']})})
        options.append({'id':key,'title':title,'targets':targets,'certificates':certs,
            'summary':{'target_records':23,'target_status_counts':dict(Counter(t['net_status'] for t in targets))}})
    required={p['compound_id'] for r in reactions for s in ('left','right') for p in r[s]} | {
        t['compound_id'] for option in options for t in option['targets']} | set(net['external_exchange_compound_ids'])
    report={'schema':'cannabis-carbon.phase1-pg-named-view.v1','view_scenario':'pg-paired-identity',
        **options[-1],'scenario_options':options,'reactions':reactions,
        'compounds':[c for c in net['compounds'] if c['id'] in required],
        'claim_boundary':net['claim_boundary'],'view_boundary':net['claim_boundary'],
        'external_exchange_compound_ids':net['external_exchange_compound_ids'],
        'co2_compound_id':net['co2_compound_id'],'forbidden_step_ids':net['forbidden_step_ids'],
        'historical_inventory_summary_unchanged':net['historical_inventory_summary_unchanged'],
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    bundle=attach_evidence(report,evidence)
    folder=Path('docs/data/pg-named-net-view');folder.mkdir(parents=True,exist_ok=True)
    payload=(json.dumps(bundle,separators=(',',':'))+'\n').encode()
    (folder/'bundle.json').write_bytes(payload)
    manifest={'schema':report['schema'],'file':'bundle.json','bytes':len(payload),
        'sha256':hashlib.sha256(payload).hexdigest(),'source_sha256':report['source_sha256']}
    (folder/'index.json').write_text(json.dumps(manifest,separators=(',',':'))+'\n')
    print(json.dumps({'bytes':len(payload),'reactions':len(reactions),'paired_targets':23}),flush=True)


if __name__=='__main__':
    run()
