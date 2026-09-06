"""Conservative straight-chain aldehyde deformylation hypotheses, not plant CER1."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

REFERENCE='balanced-equation:0f4737961dbef786ba6922fec1ae7b562c73646bba4fee5c2a06da9f6b91497d'
BOUNDARY=('Chain-length reaction-class analogy to the exact Terpedia octadecanal deformylation '
    'equation, not a demonstrated Cannabis reaction, target-specific substrate assay or plant CER1 '
    'assignment. The full NADPH, oxygen, proton, formate and water bookkeeping is retained. '
    'New aldehydes have no supply route established by this report. No coverage gain is claimed. '
    'Plant decarbonylation and cyanobacterial deformylation must not be conflated.')


def build(current,reference):
    compounds={c['id']:c for c in current['compounds']}; original=set(compounds)
    smiles={c['smiles']:c['id'] for c in compounds.values()}
    aldehyde=smiles['C'*18+'=O']; alkane=smiles['C'*17]
    if aldehyde not in {p['compound_id'] for p in reference['left']} or alkane not in {p['compound_id'] for p in reference['right']}:
        raise ValueError('Reference orientation changed')
    if not balanced([reference['left'],reference['right']],compounds):
        raise ValueError('Unbalanced reference')
    reactions={}; targets=[]
    for t in current['targets']:
        sm=compounds[t['compound_id']]['smiles']
        if t['net_status']!='no-net-producing-equation' or not sm or set(sm)!={'C'}:
            continue
        n=len(sm); precursor_smiles=sm+'C=O'; pid=stable_id('structure',precursor_smiles)
        row={'cannabisdb_id':t['cannabisdb_id'],'compound_id':t['compound_id'],'label':t['label'],
             'alkane_carbon_count':n,'aldehyde_smiles':precursor_smiles,'aldehyde_carbon_count':n+1,
             'aldehyde_compound_id':pid,'aldehyde_present_in_parent_model':pid in original,
             'status':'outside-conservative-chain-length-window-no-proposal'}
        # Deliberately bounded homolog test. This does not assert activity at every length.
        if 13<=n+1<=22:
            mol=Chem.MolFromSmiles(precursor_smiles)
            compounds.setdefault(pid,{'id':pid,'smiles':precursor_smiles,
                'formula':rdMolDescriptors.CalcMolFormula(mol),'formal_charge':0,'carbon_count':n+1})
            left=[{**p,'compound_id':pid if p['compound_id']==aldehyde else p['compound_id']} for p in reference['left']]
            right=[{**p,'compound_id':t['compound_id'] if p['compound_id']==alkane else p['compound_id']} for p in reference['right']]
            if not balanced([left,right],compounds):
                raise ValueError('Unbalanced homolog proposal')
            rid=stable_id('alkane-deformylation-hypothesis',[left,right])
            reactions[rid]={'id':rid,'left':left,'right':right,
                'hypothesis_type':'aldehyde-deformylation',
                'direction_status':'proposed-forward-only-not-Cannabis-physiological-evidence',
                'source_url':'https://iubmb.qmul.ac.uk/enzyme/EC4/1/99/5.html',
                'source_evidence_type':'reaction-class-analogy-not-target-specific-assay',
                'reference_reaction_id':REFERENCE,'enzyme_evidence_ids':[],
                'balance_status':'independently-element-isotope-charge-balanced','claim_boundary':BOUNDARY}
            row.update(status='balanced-hypothesis-with-precursor-supply-unestablished',hypothesis_id=rid)
        targets.append(row)
    used={p['compound_id'] for r in reactions.values() for side in ('left','right') for p in r[side]}
    return {'schema':'cannabis-carbon.phase1-alkane-hypotheses.v1','targets':targets,
        'reactions':list(reactions.values()),'compounds':[compounds[c] for c in sorted(used)],
        'reference_reaction':reference,'claim_boundary':BOUNDARY,
        'summary':{'full_inventory_records':len(current['targets']),'straight_chain_no_producer_records':len(targets),
                   'balanced_proposals':len(reactions),'coverage_gain_claimed':0},
        'source_scope':{'url':'https://iubmb.qmul.ac.uk/enzyme/EC4/1/99/5.html','accessed':'2026-09-06',
            'claim':'The nomenclature defines long-chain aldehyde deformylation with oxygen and reducing equivalents, producing alkane and formate. The glossary describes C13–C22 aldehydes. This work uses that as a conservative proposal window, not measured activity for all homologs.',
            'bookkeeping_note':'Exact charged structures and proton coefficient come from the independently balanced Terpedia reference, not an uncharged shorthand equation.'}}


def run():
    paths=[Path('data/reports/phase1-ketone-stereo-net.json'),Path('data/reports/phase1-full-balanced-network.json')]
    current,network=[json.loads(p.read_bytes()) for p in paths]
    reference=next(r for r in network['reactions'] if r['id']==REFERENCE)
    report=build(current,reference)
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-alkane-hypotheses.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']))


if __name__=='__main__':
    run()
