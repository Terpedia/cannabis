"""Full-inventory selenium incorporation sensitivity, with explicit acid-base bridges."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_catalog import stable_id
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_lipid_acylation_net import equation_key
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_reaction_completion_net import validate_certificate

BOUNDARY = ('Conditional stoichiometric CO2 conversion, not demonstrated Cannabis biology. '
    'KEGG R03601 supports the exact neutral source equation but its NAME mentions sulfide '
    'in conflict with its selenium DEFINITION/EQUATION. Explicit acid-base bridges are '
    'computed hypotheses, not source-curated reactions. New source-forward-only and '
    'source-reversible scenarios are kept separate; neither establishes physiological '
    'direction, enzyme capability, energetics, compartments, startup or flux. Existing '
    'direction exclusions and all external exchanges are unchanged. Regenerated '
    'pre-existing pools are allowed. CO2 remains the sole external carbon source.')


def proposals(source, compounds):
    compounds = dict(compounds)
    for c in source['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Source identity conflict')
        compounds.setdefault(c['id'], c)
    reaction = {'id': 'selenium-source:KEGG-R03601',
        **{k:source['reaction'][k] for k in ('left','right','balance_status')},
        'hypothesis_type':'selenium-incorporation', 'source_url':source['source_url'],
        'source_ec':source['reaction']['source_ec'],
        'source_evidence_type':'KEGG-exact-structure-equation-not-Cannabis-assay',
        'source_name_definition_conflict':source['source_name_definition_conflict'],
        'claim_boundary':BOUNDARY, 'enzyme_evidence_ids':[]}
    reactions = [reaction]
    # No generalized normalization: two named, fully explicit proton-only equations.
    for label, a, b in [('hydrogen-selenide','[SeH2]','[SeH-]'),
                        ('acetic-acid','CC(=O)O','CC(=O)[O-]')]:
        left = [{'compound_id':stable_id('structure',a),'coefficient':1}]
        right = [{'compound_id':stable_id('structure',s),'coefficient':1} for s in (b,'[H+]')]
        if any(p['compound_id'] not in compounds for p in left+right):
            raise ValueError('Missing exact acid-base endpoint')
        reactions.append({'id':'selenium-speciation:'+label, 'left':left,'right':right,
            'hypothesis_type':'explicit-acid-base', 'source_url':source['source_url'],
            'source_evidence_type':'computed-proton-balance-not-KEGG-curated-speciation',
            'balance_status':'independently-element-isotope-charge-balanced',
            'claim_boundary':BOUNDARY,'enzyme_evidence_ids':[]})
    for r in reactions:
        if not balanced([r['left'],r['right']], compounds):
            raise ValueError('Selenium proposal is not balanced')
    return reactions, compounds


def run():
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    parent_path=root/'phase1-local-speciation-net.json'; parent=read(parent_path)
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',
           parent_path,root/'phase1-selenium-source-identity.json']+[root/n for n in parent['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale input snapshot')
    reactions,compounds,_=assemble_current(docs[0],docs[1],parent,docs[4:])
    added,compounds=proposals(docs[3],compounds)
    keys={equation_key(r) for r in reactions.values()}
    if any(equation_key(r) in keys for r in added):
        raise ValueError('Existing equation needs direction-preserving join review')
    reactions.update({r['id']:r for r in added})
    inherited={c['compound_id']:c for doc in [*docs[4:],parent]
               for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    if set(inherited)!={t['compound_id'] for t in parent['targets'] if t['net_status']=='exact-net-conversion-hypothesis'}:
        raise ValueError('Incomplete baseline certificates')
    exchange=set(parent['external_exchange_compound_ids']); co2=parent['co2_compound_id']
    if {c for c in exchange if compounds[c]['carbon_count']}!={co2}:
        raise ValueError('External carbon boundary changed')
    participants={p['compound_id'] for r in reactions.values() for side in ('left','right') for p in r[side]}
    for mode in ('forward','reversible'):
        forbidden=sorted(set(parent['forbidden_step_ids']) | (
            {'selenium-source:KEGG-R03601:hypothetical-right-to-left'} if mode=='forward' else set()))
        model=NetModel(list(reactions.values()),exchange,forbidden_step_ids=forbidden)
        steps={s['id']:s for s in model.steps}
        for cert in inherited.values():
            validate_certificate(cert,steps,compounds,exchange,co2)
        cache={}; new={}; targets=[]
        for i,t in enumerate(parent['targets'],1):
            cid=t['compound_id']
            if cid in inherited:
                result=inherited[cid]
            else:
                if cid not in cache:
                    cache[cid]={'compound_id':cid,**model.solve(cid)}
                result=cache[cid]
                if result['status']=='exact-net-conversion-hypothesis' and cid not in new:
                    validate_certificate(result,steps,compounds,exchange,co2)
                    used=sorted({s['reaction_id'] for s in result['steps']} & {r['id'] for r in added})
                    if not used:
                        raise ValueError('New certificate lacks added chemistry')
                    new[cid]={**result,'added_lipid_reaction_ids':used}
            targets.append({k:t[k] for k in ('cannabisdb_id','label','compound_id')} | {
                'parent_net_status':t['net_status'],
                'net_status':result['status'].replace('no-net-producing-candidate-equation','no-net-producing-equation'),
                'balanced_participant':cid in participants,'new_certificate':cid in new,
                'added_lipid_reaction_ids':new.get(cid,{}).get('added_lipid_reaction_ids',[])})
            if i%500==0:
                print(f'Selenium {mode}: {i}/6220; {len(new)} new certificates',flush=True)
        used={s['reaction_id'] for c in new.values() for s in c['steps']}
        report={'schema':'cannabis-carbon.phase1-selenium-net.v1','scenario':mode,
            'summary':{'target_records':len(targets),'balanced_equations':len(reactions),
                'added_lipid_equations':len(added),'existing_equation_joins':0,
                'net_status_counts':dict(Counter(t['net_status'] for t in targets)),
                'balanced_participant_records':sum(t['balanced_participant'] for t in targets),
                'new_certificate_structures':len(new),'new_certificate_records':sum(t['new_certificate'] for t in targets)},
            'targets':targets,'new_certificates':list(new.values()),'added_reactions':added,
            'existing_equation_joins':[], 'certificate_reactions':[reactions[r] for r in sorted(used)],
            'compounds':list(compounds.values()),'external_exchange_compound_ids':sorted(exchange),
            'co2_compound_id':co2,'forbidden_step_ids':forbidden,'claim_boundary':BOUNDARY,
            'baseline_certificate_reports':parent['baseline_certificate_reports']+[parent_path.name],
            'lipid_evidence_reports':['phase1-selenium-source-identity.json'],
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
        (root/f'phase1-selenium-{mode}-net.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
        print(mode+' '+json.dumps(report['summary']),flush=True)


if __name__=='__main__':
    run()
