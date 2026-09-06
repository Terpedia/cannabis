"""Exact nonnegative stoichiometric obstruction certificates, not biological absence."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from .phase1_glycerophospholipid_input_audit import assemble_current


def validate(model,target,weights):
    w={c:Fraction(v) for c,v in weights.items()}
    if not set(w)<=set(model.internal_ids) or any(v<=0 for v in w.values()) or w.get(target,0)<1:
        raise ValueError('Invalid nonnegative internal obstruction weights')
    touched=[]
    for step in model.steps:
        delta=sum((sign*Fraction(str(p['coefficient']))*w.get(p['compound_id'],0)
                   for side,sign in [('required_inputs',-1),('outputs',1)] for p in step[side]),Fraction())
        if delta>0:
            raise ValueError('Allowed step increases obstruction weight')
        if any(p['compound_id'] in w for side in ('required_inputs','outputs') for p in step[side]):
            touched.append({'step_id':step['id'],'weighted_net_change':str(delta)})
    return touched


def solve(model,target):
    if target not in model.index:
        return {'compound_id':target,'status':'target-not-in-internal-matrix'}
    bounds=[(0,None)]*len(model.internal_ids); bounds[model.index[target]]=(1,1)
    # model.matrix = -S_internal; require S_internal.T @ weights <= 0.
    result=linprog(np.ones(len(model.internal_ids)),A_ub=-model.matrix.T,
        b_ub=np.zeros(len(model.steps)),bounds=bounds,method='highs',
        options={'time_limit':30,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    base={'compound_id':target,'solver_status':int(result.status),'solver_message':result.message}
    if not result.success:
        return {**base,'status':'no-obstruction-certificate-produced'}
    weights={c:str(Fraction(float(v)).limit_denominator(1000000)) for c,v in zip(model.internal_ids,result.x)
             if Fraction(float(v)).limit_denominator(1000000)}
    try:
        touched=validate(model,target,weights)
    except ValueError as e:
        return {**base,'status':'numerical-obstruction-failed-exact-validation','validation_error':str(e)}
    return {**base,'status':'exact-stoichiometric-obstruction','weights':weights,
            'checked_allowed_steps':len(model.steps),'weighted_steps':touched}


def run():
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    current_path=root/'phase1-alkane-net.json'; current=read(current_path)
    paths=[current_path,root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json']+[root/n for n in current['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale model input')
    reactions,compounds,model=assemble_current(docs[1],docs[2],current,docs[3:])
    targets=[t for t in current['targets'] if t['cannabisdb_id'] in ('CDB000155','CDB000157')]
    results=[]
    for target in targets:
        result={**solve(model,target['compound_id']),'cannabisdb_id':target['cannabisdb_id'],'label':target['label']}
        results.append(result); print(json.dumps({k:v for k,v in result.items() if k not in ('weights','weighted_steps')}),flush=True)
    used={c for r in results for c in r.get('weights',{})}
    report={'schema':'cannabis-carbon.phase1-net-obstruction.v1','targets':results,
        'compounds':[compounds[c] for c in sorted(used)],'forbidden_step_ids':current['forbidden_step_ids'],
        'external_exchange_compound_ids':current['external_exchange_compound_ids'],
        'summary':{'assessed_target_records':len(results),'exact_obstructions':sum(r['status']=='exact-stoichiometric-obstruction' for r in results),'balanced_equations':len(reactions),'allowed_steps':len(model.steps)},
        'proof':'Every weight is nonnegative and restricted to internal compounds; target weight is at least one. Every allowed directed reaction has weighted net change <= 0, checked with exact fractions. Nonnegative reaction extents therefore cannot produce a positive net target while leaving every other internal net amount nonnegative. External species carry zero weight.',
        'claim_boundary':'Proof applies only to the pinned reaction set, exact identities, directions and exchange boundary. It permits recycled pre-existing pools, but not their net depletion. It is not biological absence, an enzyme assignment, a minimal cut, a complete search for missing chemistry, or a minimum-medium result. Failure to produce a certificate proves nothing about feasibility.',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'phase1-alkane-obstructions.json').write_text(json.dumps(report,separators=(',',':'))+'\n')


if __name__=='__main__':
    run()
