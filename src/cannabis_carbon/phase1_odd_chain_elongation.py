"""Exact four-step homolog extensions; not established Cannabis substrate scope."""
import hashlib
import json
import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

REFERENCES=[
 'balanced-equation:f968c8cf85d6e5a896c509c818e2797237dafb30ee7a25f3bef240eb5a628715',
 'balanced-equation:87fee30ea8941caec912cf4939585d22d8e11ebcb55b4cdf2a494f0a34fc64fb',
 'balanced-equation:465c3fa43d84977c88e4a869aa70ffb26aa9afe781d3bbbdbab2f99b2efdecd7',
 'balanced-equation:0a8ee6caf6715e687dba406b6921cf609744a557b7131fc01b7b2d686aff5d9c']
BOUNDARY=('Odd-chain homolog sensitivity derived from a four-reaction C18-to-C20 acyl-CoA '
 'elongation sequence in Terpedia. Only the terminal saturated chain length is changed; '
 'all encoded intermediate and CoA stereochemistry, bond orders, cofactors and stoichiometry '
 'are retained. C17-to-C19 and C19-to-C21 substrate scope and Cannabis activity are unestablished. '
 'Every added step is proposed forward-only, not evidence of physiological direction. '
 'This report alone establishes neither precursor supply nor a CO2 pathway. No coverage gain is claimed.')


def build(current,network,shifts=(-1,1)):
    compounds={c['id']:c for c in current['compounds']}; original=set(compounds)
    refs=[next(r for r in network['reactions'] if r['id']==rid) for rid in REFERENCES]
    reactions=[]; cycles=[]
    for shift in shifts:
        changed={}; transformations=[]
        for r in refs:
            for side in ('left','right'):
                for p in r[side]:
                    cid=p['compound_id']; sm=compounds[cid]['smiles']; prefix=re.match(r'^C{17,}',sm)
                    if cid in changed or not prefix or '(=O)SCC' not in sm:
                        continue
                    n=len(prefix.group()); candidate='C'*(n+shift)+sm[n:]
                    mol=Chem.MolFromSmiles(candidate); canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
                    new_id=stable_id('structure',canonical); changed[cid]=new_id
                    compounds.setdefault(new_id,{'id':new_id,'smiles':canonical,
                        'formula':rdMolDescriptors.CalcMolFormula(mol),'formal_charge':Chem.GetFormalCharge(mol),
                        'carbon_count':sum(a.GetAtomicNum()==6 for a in mol.GetAtoms())})
                    transformations.append({'reference_compound_id':cid,'compound_id':new_id,
                        'terminal_chain_carbon_delta':shift,'present_in_parent_model':new_id in original})
        if len(changed)!=5:
            raise ValueError('Expected five exact acyl-chain structures in reference cycle')
        ids=[]
        for i,r in enumerate(refs):
            sides=[[{**p,'compound_id':changed.get(p['compound_id'],p['compound_id'])} for p in r[s]] for s in ('left','right')]
            if not balanced(sides,compounds):
                raise ValueError('Unbalanced chain-length extension')
            rid=stable_id('odd-chain-elongation-hypothesis',sides); ids.append(rid)
            reactions.append({'id':rid,'left':sides[0],'right':sides[1],
                'hypothesis_type':'odd-chain-elongation','cycle_step':i+1,
                'reference_reaction_id':r['id'],'source_url':r['sources'][0]['source_urls'][0],
                'source_evidence_type':'chain-length-analogy-not-target-substrate-or-Cannabis-evidence',
                'enzyme_evidence_ids':[],'balance_status':'independently-element-isotope-charge-balanced',
                'direction_status':'proposed-forward-only','claim_boundary':BOUNDARY})
        cycles.append({'precursor_acyl_carbons':18+shift,'product_acyl_carbons':20+shift,
                       'reaction_ids_in_order':ids,'transformations':transformations})
    used={p['compound_id'] for r in reactions for side in ('left','right') for p in r[side]}
    return {'schema':'cannabis-carbon.phase1-odd-chain-elongation.v1','cycles':cycles,'reactions':reactions,
        'reference_reactions':refs,'compounds':[compounds[c] for c in sorted(used)],'claim_boundary':BOUNDARY,
        'summary':{'cycles':len(cycles),'balanced_equations':len(reactions),'new_compound_structures':len(used-original),'coverage_gain_claimed':0}}


def run():
    paths=[Path('data/reports/phase1-'+n+'.json') for n in ('alkane-net','full-balanced-network')]
    docs=[json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source')
    report=build(*docs)
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-odd-chain-elongation.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']))


if __name__=='__main__':
    run()
