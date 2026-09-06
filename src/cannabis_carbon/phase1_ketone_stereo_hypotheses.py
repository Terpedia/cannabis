"""Explicit secondary-alcohol redox hypotheses; never a free stereoisomer edge."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced
from .phase1_stereochemistry_queue import pair_status

REFERENCE='balanced-equation:e4e2372fb33c0b3331a19d1d78dc0a5ce12980cf55b85fb347d8313196e28699'
BORNEOL='structure:38306a8d3b2da885fba17957b7b48f0df27dd9f8a811a7298930d93ab8cf41d1'
CAMPHOR='structure:f97f1ad1c50d7bb0015ed251094c5f5e38c8108e7e56d0ad2f0d309335f0584f'
BOUNDARY=('Computed secondary-alcohol/ketone redox hypothesis. Reference borneol chemistry '
    'supports the reaction class and exact cofactor bookkeeping, not substrate scope, '
    'stereoselectivity, reversibility or activity of a Cannabis enzyme. Oxidation and '
    'reduction are two separate required reactions with an exact shared ketone; no '
    'direct stereoisomerization edge is created. All other encoded stereochemistry, '
    'isotopes and connectivity are retained. Missing target-specific enzyme capability '
    'is a testable hypothesis. No CO2 route or coverage gain is claimed before full-input testing.')


def ketones(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid alcohol structure')
    result={}
    for c in mol.GetAtoms():
        if c.GetAtomicNum()!=6 or c.GetFormalCharge() or c.GetIsAromatic() or c.GetTotalNumHs()!=1 or c.GetDegree()!=3 or c.GetChiralTag()==Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        neighbors=list(c.GetNeighbors()); oxygen=[o for o in neighbors if o.GetAtomicNum()==8 and not o.GetFormalCharge() and o.GetDegree()==1 and o.GetTotalNumHs()==1]
        if len(oxygen)!=1 or sum(n.GetAtomicNum()==6 for n in neighbors)!=2:
            continue
        o=oxygen[0]
        if any(b.GetBondType()!=Chem.BondType.SINGLE for b in c.GetBonds()):
            continue
        changed=Chem.RWMol(mol)
        for index in (c.GetIdx(),o.GetIdx()):
            a=changed.GetAtomWithIdx(index); a.SetNumExplicitHs(0); a.SetNoImplicit(False)
        changed.GetAtomWithIdx(c.GetIdx()).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        changed.GetBondBetweenAtoms(c.GetIdx(),o.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        product=changed.GetMol(); Chem.SanitizeMol(product)
        key=Chem.MolToSmiles(product,isomericSmiles=True)
        result[key]={'alcohol_carbon_atom_index':c.GetIdx(),'oxygen_atom_index':o.GetIdx()}
    return result


def build(queue,current,reference):
    compounds={c['id']:c for c in current['compounds']}
    original_ids=set(compounds)
    if CAMPHOR not in {p['compound_id'] for p in reference['left']} or BORNEOL not in {p['compound_id'] for p in reference['right']}:
        raise ValueError('Reference sides changed')
    if not balanced([reference['left'],reference['right']],compounds):
        raise ValueError('Unbalanced reference')
    if compounds[CAMPHOR]['smiles'] not in ketones(compounds[BORNEOL]['smiles']):
        raise ValueError('Reference does not match exact local oxidation')
    proposals={}; pairs=[]
    for t in queue['targets']:
        if t['source_identity_conflict']:
            continue
        cid=t['compound_id']; target_products=ketones(compounds[cid]['smiles'])
        if not target_products:
            continue
        for alternative in t['alternatives']:
            other=alternative['compound_id']
            if pair_status(compounds[cid]['smiles'],compounds[other]['smiles'])!='specified-stereoisomer-difference':
                continue
            partner_products=ketones(compounds[other]['smiles'])
            for smiles in sorted(target_products.keys() & partner_products.keys()):
                kid=stable_id('structure',smiles); mol=Chem.MolFromSmiles(smiles)
                existing=kid in original_ids
                compounds.setdefault(kid,{'id':kid,'smiles':smiles,'formula':rdMolDescriptors.CalcMolFormula(mol),
                    'formal_charge':Chem.GetFormalCharge(mol),'carbon_count':sum(a.GetAtomicNum()==6 for a in mol.GetAtoms())})
                ids=[]
                for alcohol,reverse in ((other,True),(cid,False)):
                    sides=[[{**p,'compound_id':kid if p['compound_id']==CAMPHOR else p['compound_id']} for p in reference['left']],
                           [{**p,'compound_id':alcohol if p['compound_id']==BORNEOL else p['compound_id']} for p in reference['right']]]
                    if reverse:
                        sides.reverse()
                    if not balanced(sides,compounds):
                        raise ValueError('Local redox hypothesis failed exact balance')
                    rid=stable_id('ketone-stereo-hypothesis',sides); ids.append(rid)
                    proposals[rid]={'id':rid,'left':sides[0],'right':sides[1],
                        'hypothesis_type':'secondary-alcohol-oxidation' if reverse else 'stereoselective-ketone-reduction',
                        'source_url':'https://enzyme.expasy.org/EC/1.1.1.198',
                        'reference_reaction_id':reference['id'],
                        'source_evidence_type':'reaction-class-analogy-not-target-substrate-or-Cannabis-evidence',
                        'direction_status':'proposed-forward-sensitivity-only-not-physiological-direction',
                        'balance_status':'independently-element-isotope-charge-balanced','enzyme_evidence_ids':[],
                        'claim_boundary':BOUNDARY}
                pairs.append({'cannabisdb_id':t['cannabisdb_id'],'label':t['label'],'compound_id':cid,
                    'partner_compound_id':other,'ketone_compound_id':kid,'ketone_already_in_model':existing,
                    'partner_has_saved_conditional_CO2_certificate':alternative['partner_has_saved_conditional_CO2_certificate'],
                    'partner_oxidation_hypothesis_id':ids[0],'target_reduction_hypothesis_id':ids[1],
                    'target_local_change':target_products[smiles],'partner_local_change':partner_products[smiles]})
    used={p['compound_id'] for r in proposals.values() for side in ('left','right') for p in r[side]}
    return {'schema':'cannabis-carbon.phase1-ketone-stereo-hypotheses.v1','pairs':pairs,
        'reactions':list(proposals.values()),'compounds':[compounds[c] for c in sorted(used)],
        'reference_reaction':reference,'claim_boundary':BOUNDARY,
        'summary':{'target_records':len({p['cannabisdb_id'] for p in pairs}),'pair_hypotheses':len(pairs),
            'balanced_equation_proposals':len(proposals),'reaction_type_counts':dict(Counter(r['hypothesis_type'] for r in proposals.values())),
            'coverage_gain_claimed':0}}


def run():
    RDLogger.DisableLog('rdApp.warning'); root=Path('data/reports')
    paths=[root/'phase1-stereochemistry-queue.json',root/'phase1-selenium-forward-net.json',root/'phase1-full-balanced-network.json']
    queue,current,network=[json.loads(p.read_bytes()) for p in paths]
    for doc in (queue,current,network):
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source')
    reference=next(r for r in network['reactions'] if r['id']==REFERENCE)
    report=build(queue,current,reference)
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root/'phase1-ketone-stereo-hypotheses.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)
    print(json.dumps(report['pairs']),flush=True)


if __name__=='__main__':
    run()
