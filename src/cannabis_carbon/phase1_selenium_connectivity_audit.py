"""Necessary selenium-connectivity cut; relaxed paths are not full pathways."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_glycerophospholipid_input_audit import assemble_current


def selenium_cut(steps, counts, exchange):
    available = {cid for cid in exchange if counts[cid]}
    initial = sorted(available)
    relevant = []
    for s in steps:
        inputs = {p['compound_id'] for p in s['required_inputs'] if counts[p['compound_id']]}
        outputs = {p['compound_id'] for p in s['outputs'] if counts[p['compound_id']]}
        if outputs and not inputs:
            raise ValueError('Selenium created without selenium input')
        if inputs or outputs:
            relevant.append((s, inputs, outputs))
    while True:
        new = set().union(*(outputs for _,inputs,outputs in relevant if inputs & available)) if relevant else set()
        if new <= available:
            break
        available.update(new)
    unreachable = {cid for cid,n in counts.items() if n and cid not in available}
    deltas = []
    for s, _, _ in relevant:
        delta = sum(Fraction(p['coefficient'])*counts[p['compound_id']] for p in s['outputs'] if p['compound_id'] in unreachable) - sum(
            Fraction(p['coefficient'])*counts[p['compound_id']] for p in s['required_inputs'] if p['compound_id'] in unreachable)
        if delta > 0:
            raise ValueError('Invalid necessary cut: an allowed step increases its weight')
        deltas.append({'step_id':s['id'],'unreachable_selenium_weight_delta':str(delta)})
    return {'external_selenium_compound_ids':initial, 'relaxed_reachable_compound_ids':sorted(available),
            'unreachable_compound_ids':sorted(unreachable), 'step_weight_deltas':deltas}


def run():
    RDLogger.DisableLog('rdApp.warning')
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    current=read(root/'phase1-local-speciation-net.json')
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',
           root/'phase1-local-speciation-net.json']+[root/n for n in current['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source snapshot')
    reactions,compounds,model=assemble_current(docs[0],docs[1],current,docs[3:])
    counts={cid:sum(a.GetAtomicNum()==34 for a in Chem.MolFromSmiles(c['smiles']).GetAtoms()) for cid,c in compounds.items()}
    cut=selenium_cut(model.steps,counts,current['external_exchange_compound_ids'])
    unreachable=set(cut['unreachable_compound_ids'])
    targets=[{**t,'selenium_count':counts[t['compound_id']],
              'selenium_connectivity_status':'unreachable-even-in-relaxed-selenium-graph' if t['compound_id'] in unreachable else 'relaxed-connectivity-only-not-a-pathway'}
             for t in current['targets'] if counts[t['compound_id']]]
    if any(t['net_status']=='exact-net-conversion-hypothesis' and t['compound_id'] in unreachable for t in targets):
        raise ValueError('Cut contradicts a saved certificate')
    se_ids={cid for cid,n in counts.items() if n}
    selected=[r for r in reactions.values() if any(p['compound_id'] in se_ids for side in ('left','right') for p in r[side])]
    used={p['compound_id'] for r in selected for side in ('left','right') for p in r[side]}|se_ids
    report={'schema':'cannabis-carbon.phase1-selenium-connectivity-audit.v1', **cut, 'targets':targets,
        'compounds':[compounds[cid] for cid in sorted(used)], 'reactions':selected,
        'selenium_atom_counts':{cid:counts[cid] for cid in sorted(used)},
        'forbidden_step_ids':current['forbidden_step_ids'],
        'external_exchange_compound_ids':current['external_exchange_compound_ids'],
        'summary':{'whole_model_equations':len(reactions),'whole_model_directed_steps':len(model.steps),
            'selenium_structures':len(se_ids),'selenium_target_records':len(targets),
            'target_status_counts':dict(Counter(t['selenium_connectivity_status'] for t in targets)),
            'coverage_gain_claimed':0},
        'claim_boundary':'Necessary model obstruction only, not biological absence. Relaxed selenium connectivity ignores all non-selenium prerequisites and uses any selenium input to connect all selenium outputs. Reachability is therefore an upper bound, never a complete reaction pathway. For unreachable structures, the number of selenium atoms in the unreachable set cannot increase in any allowed step; every external exchange has zero such weight. Nondepleting internal pools cannot supply a positive export of these structures. No direction or exchange is changed. Review source-backed selenium-incorporation chemistry rather than assume mineral supply is absent.',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'phase1-selenium-connectivity-audit.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)
    print(json.dumps([{'id':t['cannabisdb_id'],'label':t['label'],'status':t['selenium_connectivity_status']} for t in targets]),flush=True)


if __name__=='__main__':
    run()
