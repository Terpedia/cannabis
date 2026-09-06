"""Separate net uptake permissions from disposal; no physiological-medium claim."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_net_flux import NetModel
from .phase1_reaction_completion_net import validate_certificate


class UptakeLimitedModel:
    def __init__(self, reactions, exchange_ids, uptake_ids, forbidden_step_ids=()):
        self.exchange_ids=set(exchange_ids)
        self.uptake_ids=set(uptake_ids)
        if not self.uptake_ids<=self.exchange_ids:
            raise ValueError('Uptake must be a subset of declared external species')
        # NetModel constrains every non-input species to nonnegative net production,
        # not zero. Thus blocked inputs may still be produced and disposed of.
        self.model=NetModel(reactions,self.uptake_ids,forbidden_step_ids=forbidden_step_ids)
        self.steps=self.model.steps

    def solve(self,target):
        if target in self.exchange_ids:
            return {'status':'external-inventory-species-not-assessed-for-synthesis',
                    'uptake_allowed':target in self.uptake_ids}
        result=self.model.solve(target)
        if result['status']=='exact-net-conversion-hypothesis':
            if not set(result['external_net_consumption'])<=self.uptake_ids:
                raise ValueError('Forbidden uptake in exact witness')
            result['zero_net_internal_participants']=[c for c in result['zero_net_internal_participants']
                                                     if c not in self.exchange_ids]
        return result


def blocked_species(compounds,exchange):
    rows=[]
    for cid in sorted(exchange):
        c=compounds[cid]; mol=Chem.MolFromSmiles(c['smiles'])
        atoms=Counter(a.GetSymbol() for a in mol.GetAtoms())
        if atoms['Fe']>1 and atoms['S']>1:
            rows.append({'compound_id':cid,'smiles':c['smiles'],'reason':'iron-sulfur-cluster-input'})
        elif c['smiles']=='OO':
            rows.append({'compound_id':cid,'smiles':c['smiles'],'reason':'hydrogen-peroxide-input'})
    return rows


def run():
    RDLogger.DisableLog('rdApp.warning')
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    current_path=root/'phase1-selenium-forward-net.json'; parent=read(current_path)
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',current_path]+[root/n for n in parent['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale model input')
    reactions,compounds,_=assemble_current(docs[0],docs[1],parent,docs[3:])
    exchange=set(parent['external_exchange_compound_ids']); blocked=blocked_species(compounds,exchange)
    uptake=exchange-{r['compound_id'] for r in blocked}
    model=UptakeLimitedModel(list(reactions.values()),exchange,uptake,parent['forbidden_step_ids'])
    steps={s['id']:s for s in model.steps}; co2=parent['co2_compound_id']
    prior={c['compound_id']:c for doc in [*docs[3:],parent] for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    if set(prior)!={t['compound_id'] for t in parent['targets'] if t['net_status']=='exact-net-conversion-hypothesis'}:
        raise ValueError('Missing baseline witnesses')
    cache={}; certs={}; targets=[]; reused=0
    for i,t in enumerate(parent['targets'],1):
        cid=t['compound_id']
        if cid not in cache:
            if cid in prior and set(prior[cid]['external_net_consumption'])<=uptake:
                result=prior[cid]; reused+=1
            else:
                result={'compound_id':cid,**model.solve(cid)}
            if result['status']=='exact-net-conversion-hypothesis':
                validate_certificate(result,steps,compounds,uptake,co2)
                certs[cid]=result
            cache[cid]=result
        result=cache[cid]
        targets.append({k:t[k] for k in ('cannabisdb_id','label','compound_id')} | {
            'parent_net_status':t['net_status'],'net_status':result['status'],
            'prior_witness_reusable':cid in prior and set(prior[cid]['external_net_consumption'])<=uptake,
            'uptake_allowed':cid in uptake if cid in exchange else None})
        if i%250==0:
            print(f'Medium boundary {i}/6220: {len(certs)} exact witnesses; {reused} reused',flush=True)
    report={'schema':'cannabis-carbon.phase1-medium-boundary.v1',
        'scenario':'no-net-uptake-of-iron-sulfur-clusters-or-hydrogen-peroxide',
        'targets':targets,'certificates':list(certs.values()),'blocked_inputs':blocked,
        'allowed_uptake_compound_ids':sorted(uptake),'allowed_external_output_compound_ids':sorted(exchange),
        'forbidden_step_ids':parent['forbidden_step_ids'],'co2_compound_id':co2,
        'summary':{'target_records':len(targets),'balanced_equations':len(reactions),
            'original_external_species':len(exchange),'allowed_uptake_species':len(uptake),
            'blocked_uptake_species':len(blocked),'certificate_structures':len(certs),
            'reused_certificate_structures':reused,
            'net_status_counts':dict(Counter(t['net_status'] for t in targets)),
            'lost_covered_records':sum(t['parent_net_status']=='exact-net-conversion-hypothesis' and t['net_status']!='exact-net-conversion-hypothesis' for t in targets)},
        'claim_boundary':'Diagnostic uptake restriction, not a curated plant medium or minimum medium. Only net import of explicitly listed iron-sulfur clusters and hydrogen peroxide is blocked; their production and disposal remain allowed. No reactions, source identities or reaction directions are changed. Other permissive external inputs remain. Full inventory retained, with external inventory entries separately unassessed as synthesis targets. Exact witnesses allow recycled pre-existing pools, not zero-pool startup, growth, transport, compartments or energetically feasible metabolism. Lost coverage diagnoses this model boundary, not biological absence or nutrient essentiality.',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'phase1-medium-boundary.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__=='__main__':
    run()
