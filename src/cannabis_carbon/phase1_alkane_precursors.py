"""Explicit acyl-CoA reduction hypotheses upstream of alkane deformylation."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

REFERENCE='balanced-equation:69c53ec6cf79e634b00f942c2f7084205466d9655ed5ccc003cd4520e0ee78b1'
BOUNDARY=('Homolog-extension hypothesis for acyl-CoA reduction to a free aldehyde. Exact NADPH, '
    'proton and free CoA bookkeeping is retained from the Terpedia hexadecanoyl-CoA equation. '
    'Activity at C19-C22, Cannabis enzyme identity and physiological direction are unestablished. '
    'The current IUBMB EC 1.2.1.50 describes an acyl-protein complex and separately notes acceptance '
    'of acyl-CoA; the EC is a reaction-class reference, not an exact target enzyme assignment. '
    'No protein/ACP intermediate is merged with CoA. Acyl-CoA presence is not a supply certificate.')


def build(current,alkanes,reference):
    compounds={c['id']:c for c in current['compounds']}; original=set(compounds)
    compounds.update({c['id']:c for c in alkanes['compounds']})
    bysm={c['smiles']:c['id'] for c in compounds.values()}
    aid=bysm['C'*16+'=O']
    acyl=next(p['compound_id'] for p in reference['left'] if compounds[p['compound_id']]['smiles'].startswith('C'*16+'(=O)SCC'))
    tail=compounds[acyl]['smiles'][16:]
    if aid not in {p['compound_id'] for p in reference['right']} or not balanced([reference['left'],reference['right']],compounds):
        raise ValueError('Reference reduction orientation or balance changed')
    reactions=list(alkanes['reactions']); pairs=[]
    for target in alkanes['targets']:
        if 'hypothesis_id' not in target:
            continue
        n=target['aldehyde_carbon_count']; sm=Chem.MolToSmiles(Chem.MolFromSmiles('C'*n+tail),isomericSmiles=True)
        cid=stable_id('structure',sm); mol=Chem.MolFromSmiles(sm)
        compounds.setdefault(cid,{'id':cid,'smiles':sm,'formula':rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge':Chem.GetFormalCharge(mol),'carbon_count':sum(a.GetAtomicNum()==6 for a in mol.GetAtoms())})
        left=[{**p,'compound_id':cid if p['compound_id']==acyl else p['compound_id']} for p in reference['left']]
        right=[{**p,'compound_id':target['aldehyde_compound_id'] if p['compound_id']==aid else p['compound_id']} for p in reference['right']]
        if not balanced([left,right],compounds):
            raise ValueError('Acyl-CoA reduction failed exact balance')
        rid=stable_id('alkane-acyl-reduction-hypothesis',[left,right])
        reactions.append({'id':rid,'left':left,'right':right,'hypothesis_type':'aldehyde-forming-acyl-CoA-reduction',
            'source_url':'https://iubmb.qmul.ac.uk/enzyme/EC1/2/1/50.html',
            'source_evidence_type':'reaction-class-analogy-not-target-specific-assay',
            'reference_reaction_id':REFERENCE,'enzyme_evidence_ids':[],
            'direction_status':'proposed-forward-only-not-Cannabis-physiological-evidence',
            'balance_status':'independently-element-isotope-charge-balanced','claim_boundary':BOUNDARY})
        pairs.append({**target,'acyl_CoA_compound_id':cid,'acyl_CoA_present_in_parent_model':cid in original,
            'aldehyde_reduction_hypothesis_id':rid,'supply_status':'requires-full-input-net-test'})
    used={p['compound_id'] for r in reactions for side in ('left','right') for p in r[side]}
    return {'schema':'cannabis-carbon.phase1-alkane-precursors.v1','pairs':pairs,'reactions':reactions,
        'compounds':[compounds[c] for c in sorted(used)],'reference_reaction':reference,
        'claim_boundary':BOUNDARY+' '+alkanes['claim_boundary'],
        'summary':{'target_records':len(pairs),'balanced_equations':len(reactions),
            'new_acyl_CoA_structures':len({p['acyl_CoA_compound_id'] for p in pairs if not p['acyl_CoA_present_in_parent_model']}),
            'coverage_gain_claimed':0}}


def run():
    paths=[Path('data/reports/phase1-'+n+'.json') for n in ('ketone-stereo-net','alkane-hypotheses','full-balanced-network')]
    current,alkanes,network=[json.loads(p.read_bytes()) for p in paths]
    for doc in (current,alkanes,network):
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale input')
    report=build(current,alkanes,next(r for r in network['reactions'] if r['id']==REFERENCE))
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-alkane-precursors.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']))


if __name__=='__main__':
    run()
