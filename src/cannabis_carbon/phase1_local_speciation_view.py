"""Evidence-preserving full-inventory view with certificates loaded on demand."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_chemistry_route_view import reaction_sources
from .phase1_net_view import build as attach_evidence
from .phase1_sharded_net_view import write_view


def run():
    paths = [Path('data/reports/phase1-local-speciation-net.json'),
             Path('docs/data/glycerophospholipid-net-view/bundle.json')] + [
                 Path('data/reports/'+n+'.json') for n in
                 ('phase1-target-hypotheses','phase1-screened-enzyme-overlay','phase1-route-enzyme-overlay')]
    current, previous, *evidence = [json.loads(p.read_bytes()) for p in paths]
    for doc in (current, previous):
        for p, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale view input')
    certificates = previous['certificates']+current['new_certificates']
    bycert = {c['compound_id']:c for c in certificates}
    if len(bycert)!=len(certificates):
        raise ValueError('Duplicate exact certificates')
    reactions = {r['id']:r for r in previous['reactions']}
    for r in current['certificate_reactions']:
        if r['id'] in reactions:
            if any(reactions[r['id']][s]!=r[s] for s in ('left','right')):
                raise ValueError('Changed prior equation')
            continue
        item = {**r, 'enzyme_evidence_ids':r.get('enzyme_evidence_ids', [])}
        if r.get('hypothesis_type'):
            item.update({'is_route_sensitivity':True, 'hypothesis_assumptions':[r['claim_boundary']]})
            if r['hypothesis_type'] in ('source-mapped-protonation','amino-phospholipid-speciation','glycerophospholipid-speciation'):
                item['sources'] = reaction_sources(r,{})
            elif r.get('source_url'):
                item['sources'] = r.get('sources', []) + [{
                    'source_urls':[r['source_url']], 'evidence_type':r.get('source_evidence_type','computed-reaction-hypothesis-not-demonstrated-Cannabis-activity'),
                    'claim_boundary':r['claim_boundary']}]
        if not item.get('sources'):
            raise ValueError('Missing display provenance: '+r['id'])
        item['missing_candidate_evidence'] = not bool(item['enzyme_evidence_ids'])
        reactions[r['id']] = item
    targets = [{**t, 'certificate_compound_id':t['compound_id'] if t['compound_id'] in bycert else None,
        'startup_status':'not established by this net certificate',
        'missing_candidate_reaction_ids':sorted({s['reaction_id'] for s in bycert.get(t['compound_id'],{}).get('steps',[])
                                               if not reactions[s['reaction_id']]['enzyme_evidence_ids']})}
        for t in current['targets']]
    used = {p['compound_id'] for r in reactions.values() for side in ('left','right') for p in r[side]}
    report = {'schema':'cannabis-carbon.local-speciation-view.v1', 'view_scenario':'reaction-first-chemistry',
        'targets':targets, 'certificates':certificates, 'reactions':list(reactions.values()),
        'compounds':[c for c in current['compounds'] if c['id'] in used],
        'summary':{'target_records':len(targets), 'target_status_counts':dict(Counter(t['net_status'] for t in targets))},
        'model_summary':current['summary'], 'forbidden_step_ids':current['forbidden_step_ids'],
        'external_exchange_compound_ids':current['external_exchange_compound_ids'], 'co2_compound_id':current['co2_compound_id'],
        'view_boundary':'Reaction-first sensitivity across all 6,220 historical records. Local proton relocation and H+ exchanges are explicit reversible assumptions, not curated reactions or established tissue speciation.',
        'claim_boundary':current['claim_boundary'],
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    bundle = attach_evidence(report,evidence)
    print(json.dumps(write_view(bundle,'docs/data/local-speciation-net-view')),flush=True)


if __name__=='__main__':
    run()
