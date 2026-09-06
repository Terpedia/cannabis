"""Extend the on-demand map with exact selenium chemistry and medium annotations."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_sharded_net_view import write_view
from .phase1_net_view import build as attach_evidence


def run():
    folder=Path('docs/data/local-speciation-net-view')
    read=lambda p:json.loads(p.read_bytes())
    manifest=read(folder/'index.json')
    def verified(ref):
        payload=(folder/ref['file']).read_bytes()
        if len(payload)!=ref['bytes'] or hashlib.sha256(payload).hexdigest()!=ref['sha256']:
            raise ValueError('Prior map hash mismatch')
        return json.loads(payload)
    base=verified(manifest); shared=verified(base['shared_chemistry'])
    certs=[verified(ref) for ref in base['certificate_files'].values()]
    paths=[Path('data/reports/phase1-selenium-forward-net.json'),
           Path('data/reports/phase1-medium-inventory.json'),Path('data/curation/light-reaction-requirements.json'),
           folder/'index.json',folder/'bundle.json']+[Path('data/reports/'+n+'.json') for n in
           ('phase1-target-hypotheses','phase1-screened-enzyme-overlay','phase1-route-enzyme-overlay')]
    current,medium,light=[read(p) for p in paths[:3]]
    for doc in (current,medium,base):
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source')
    certs+=current['new_certificates']; bycert={c['compound_id']:c for c in certs}
    if len(bycert)!=len(certs):
        raise ValueError('Duplicate exact witness')
    reactions={r['id']:r for r in shared['reactions']}
    for r in current['certificate_reactions']:
        if r['id'] in reactions:
            if any(r[s]!=reactions[r['id']][s] for s in ('left','right')):
                raise ValueError('Changed prior reaction')
            continue
        sources=r.get('sources',[])
        if r.get('source_url'):
            sources=sources+[{'source_urls':[r['source_url']],'evidence_type':r['source_evidence_type'],
                              'claim_boundary':r['claim_boundary']}]
        if not sources:
            raise ValueError('Missing new reaction provenance')
        reactions[r['id']]={**r,'enzyme_evidence_ids':r.get('enzyme_evidence_ids',[]),
            'is_route_sensitivity':bool(r.get('hypothesis_type')),
            'hypothesis_assumptions':[r['claim_boundary']] if r.get('hypothesis_type') else [],
            'missing_candidate_evidence':not bool(r.get('enzyme_evidence_ids')),
            'sources':sources}
    for annotation in light['reactions']:
        rid=annotation['model_reaction_id']
        if rid in reactions:
            reactions[rid]={**reactions[rid],'light_requirement_annotation':annotation}
    labels={c['id']:c.get('labels',[]) for c in shared['compounds']}
    for t in current['targets']:
        if t['label'] not in labels.setdefault(t['compound_id'],[]):
            labels[t['compound_id']].append(t['label'])
    inputs={c['id']:c for c in medium['exchange_species']}
    used={p['compound_id'] for r in reactions.values() for side in ('left','right') for p in r[side]}
    compounds=[]
    for c in current['compounds']:
        if c['id'] not in used:
            continue
        item={**c,'labels':labels.get(c['id'],[])}
        if c['id'] in inputs:
            row=inputs[c['id']]
            item['medium_annotation']={k:row[k] for k in ('review_flags','consuming_structure_count','nutrient_essentiality')}
        compounds.append(item)
    targets=[{**t,'certificate_compound_id':t['compound_id'] if t['compound_id'] in bycert else None,
        'missing_candidate_reaction_ids':sorted({s['reaction_id'] for s in bycert.get(t['compound_id'],{}).get('steps',[])
                                               if not reactions[s['reaction_id']]['enzyme_evidence_ids']})}
        for t in current['targets']]
    report={'schema':'cannabis-carbon.selenium-view.v1','view_scenario':'reaction-first-chemistry',
        'targets':targets,'certificates':certs,'reactions':list(reactions.values()),'compounds':compounds,
        'enzyme_evidence':shared['enzyme_evidence'],
        'summary':{'target_records':len(targets),'target_status_counts':dict(Counter(t['net_status'] for t in targets))},
        'model_summary':current['summary'],'medium_summary':medium['summary'],
        'light_reaction_requirements':light,
        'forbidden_step_ids':current['forbidden_step_ids'],
        'external_exchange_compound_ids':current['external_exchange_compound_ids'],'co2_compound_id':current['co2_compound_id'],
        'view_boundary':'Selenium incorporation is source-forward only. Node details annotate external-input dependencies; edge details retain source light requirements where reviewed. Those light annotations are not enforced energy constraints. The permissive 102-species boundary is not a minimum defined medium.',
        'claim_boundary':current['claim_boundary'],
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    report=attach_evidence(report,[read(p) for p in paths[5:]])
    Path('docs/data/light-reaction-requirements.json').write_bytes(paths[2].read_bytes())
    print(json.dumps(write_view(report,'docs/data/selenium-net-view')),flush=True)


if __name__=='__main__':
    run()
