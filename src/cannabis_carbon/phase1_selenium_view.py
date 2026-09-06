"""Extend the on-demand map with exact selenium chemistry and medium annotations."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_sharded_net_view import write_view
from .phase1_net_view import build as attach_evidence


def run(*, ketone=False, c17=False):
    if ketone and c17:
        raise ValueError('Select one scenario')
    folder=Path('docs/data/ketone-stereo-net-view' if c17 else 'docs/data/selenium-net-view' if ketone else 'docs/data/local-speciation-net-view')
    read=lambda p:json.loads(p.read_bytes())
    manifest=read(folder/'index.json')
    def verified(ref):
        payload=(folder/ref['file']).read_bytes()
        if len(payload)!=ref['bytes'] or hashlib.sha256(payload).hexdigest()!=ref['sha256']:
            raise ValueError('Prior map hash mismatch')
        return json.loads(payload)
    base=verified(manifest); shared=verified(base['shared_chemistry'])
    certs=[verified(ref) for ref in base['certificate_files'].values()]
    paths=[Path('data/reports/phase1-c17-elongation-net.json' if c17 else 'data/reports/phase1-ketone-stereo-net.json' if ketone else 'data/reports/phase1-selenium-forward-net.json'),
           Path('data/reports/phase1-medium-inventory.json'),Path('data/curation/light-reaction-requirements.json'),
           folder/'index.json',folder/'bundle.json']+[Path('data/reports/'+n+'.json') for n in
           ('phase1-target-hypotheses','phase1-screened-enzyme-overlay','phase1-route-enzyme-overlay')]
    current,medium,light=[read(p) for p in paths[:3]]
    extra_layers=[]
    audit=None
    if c17:
        extra_paths=[Path('data/reports/phase1-alkane-net.json'),Path('data/reports/phase1-odd-chain-net.json'),
                     Path('data/reports/phase1-c17-evidence-audit.json')]
        extra_layers=[read(p) for p in extra_paths[:2]]
        audit=read(extra_paths[2])
        paths+=extra_paths
    completion_path=Path('data/reports/phase1-marts-completions.json')
    completion_report=read(completion_path)
    completions={c['id']:c for c in completion_report['completions']}
    variants={v['id']:v for v in completion_report['variants']}
    for doc in (current,medium,base,*extra_layers,*([audit] if audit else [])):
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source')
    for layer in extra_layers:
        certs+=layer['new_certificates']
    certs+=current['new_certificates']; bycert={c['compound_id']:c for c in certs}
    if len(bycert)!=len(certs):
        raise ValueError('Duplicate exact witness')
    reactions={r['id']:r for r in shared['reactions']}
    for r in [r for layer in [*extra_layers,current] for r in layer['certificate_reactions']]:
        if r['id'] in reactions:
            if any(r[s]!=reactions[r['id']][s] for s in ('left','right')):
                raise ValueError('Changed prior reaction')
            continue
        sources=r.get('sources',[])
        if r.get('source_url'):
            source_type=r.get('source_evidence_type')
            if source_type is None and r.get('hypothesis_type')=='source-mapped-protonation':
                source_type='source-mapped-acid-base-hypothesis-not-curated-biochemical-reaction'
            if source_type is None:
                raise ValueError('Missing source evidence classification')
            sources=sources+[{'source_urls':[r['source_url']],'evidence_type':source_type,
                              'claim_boundary':r['claim_boundary']}]
        if not sources and r.get('completion_ids'):
            records=[completions[cid] for cid in r['completion_ids']]
            for record in records:
                if record['balanced_equation_id']!=r['id'] or any(record[s]!=r[s] for s in ('left','right')):
                    raise ValueError('Completion provenance equation mismatch')
            sources=[{'evidence_type':'inferred-stoichiometric-completion',
                      'completion':record,'source_variant':variants[record['variant_id']],
                      'claim_boundary':record['claim_boundary']} for record in records]
            r={**r,'hypothesis_type':'inferred-stoichiometric-completion',
               'claim_boundary':completion_report['claim_boundary']}
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
    if audit:
        for cert in audit['certificates']:
            if cert['role']!='new-inventory-target-certificate':
                continue
            for s in cert['steps']:
                rid=s['reaction_id']
                if rid not in reactions:
                    raise ValueError('Audited reaction missing from map')
                note={k:s[k] for k in ('evidence_class','review_flags',
                    'cannabis_physiological_direction_established_by_this_audit','enzyme_assignment_established_by_this_audit')}
                notes=reactions[rid].setdefault('certificate_direction_annotations',{})
                if s['direction_mode'] in notes and notes[s['direction_mode']]!=note:
                    raise ValueError('Conflicting direction audit')
                notes[s['direction_mode']]=note
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
    report['source_sha256'][str(completion_path)]=hashlib.sha256(completion_path.read_bytes()).hexdigest()
    report=attach_evidence(report,[read(p) for p in paths[5:8]])
    if ketone:
        report['schema']='cannabis-carbon.ketone-stereo-view.v1'
        report['view_boundary']=('Separate ketone-mediated redox hypothesis scenario. Newly added reactions are '
            'forward-only reaction-class analogies, not target-specific enzyme evidence. Medium annotations '
            'and consumer counts describe the earlier selenium-forward certificate set, not this scenario. '
            'The unchanged permissive exchange boundary is not a minimum defined medium.')
        report['medium_annotation_scope']='selenium-forward baseline; not recomputed for ketone hypotheses'
    if c17:
        report['schema']='cannabis-carbon.c17-view.v1'
        report['view_boundary']=('Separate alkane and odd-chain elongation hypothesis scenario. '
            'New steps are proposed forward-only chain-length analogies, not confirmed Cannabis activity. '
            'Direction review flags describe the four newly audited target certificates, not a full-network physiology audit. '
            'Medium annotations and consumer counts refer to the earlier selenium-forward certificate set. '
            'The unchanged permissive 102-species boundary is not a minimum defined medium.')
        report['medium_annotation_scope']='selenium-forward baseline; not recomputed for C17 hypotheses'
        report['certificate_evidence_audit_summary']=audit['summary']
    Path('docs/data/light-reaction-requirements.json').write_bytes(paths[2].read_bytes())
    print(json.dumps(write_view(report,'docs/data/c17-net-view' if c17 else 'docs/data/ketone-stereo-net-view' if ketone else 'docs/data/selenium-net-view')),flush=True)


if __name__=='__main__':
    run()
