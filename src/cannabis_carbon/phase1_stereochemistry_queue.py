"""Existing stereo reaction evidence and missing-producer searches, never auto-edges."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from .phase1_current_reaction_gaps import producing_steps
from .phase1_glycerophospholipid_input_audit import assemble_current


@lru_cache(maxsize=None)
def identity(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):
        raise ValueError('Concrete structure required')
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    stereo=Chem.FindPotentialStereo(mol,cleanIt=True,flagPossible=True)
    canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
    counts=Counter(str(s.specified) for s in stereo)
    stripped=Chem.Mol(mol); Chem.RemoveStereochemistry(stripped)
    return {'canonical_smiles':canonical,'constitution_key':Chem.MolToSmiles(stripped,isomericSmiles=True),
        'potential_stereo_count':len(stereo),'stereo_status_counts':dict(counts),
        'fully_specified':all(str(s.specified)=='Specified' for s in stereo)}


def pair_status(a,b):
    x,y=identity(a),identity(b)
    if x['constitution_key']!=y['constitution_key'] or x['canonical_smiles']==y['canonical_smiles']:
        return None
    return 'specified-stereoisomer-difference' if x['fully_specified'] and y['fully_specified'] else 'unspecified-stereo-identity-review'


def build(current,compounds,reactions,conflicts):
    allowed=producing_steps(list(reactions.values()),current['forbidden_step_ids'])
    existing=[]; used=set(); selected=set()
    for r in reactions.values():
        net=defaultdict(Fraction)
        for side,sign in [('left',-1),('right',1)]:
            for p in r[side]:
                net[p['compound_id']]+=sign*Fraction(p['coefficient'])
        left=[c for c,n in net.items() if n<0 and compounds[c]['carbon_count']]
        right=[c for c,n in net.items() if n>0 and compounds[c]['carbon_count']]
        for a in left:
            for b in right:
                if -net[a]!=net[b]:
                    continue
                status=pair_status(compounds[a]['smiles'],compounds[b]['smiles'])
                if status:
                    existing.append({'reaction_id':r['id'],'left_compound_id':a,'right_compound_id':b,
                        'pair_status':status,'pair_coefficient':str(net[b]),
                        'allowed_producing_directions':{c:sorted(s for s in allowed[c] if s.rsplit(':',1)[0]==r['id']) for c in (a,b)},
                        'claim_boundary':'Pair embedded in an existing full equation. Other participants remain required. Structural classification alone does not prove enzyme capability or Cannabis occurrence.'})
                    used.update((a,b));selected.add(r['id'])
    indexes=defaultdict(set)
    for cid in allowed:
        if allowed[cid] and compounds[cid]['carbon_count']:
            indexes[identity(compounds[cid]['smiles'])['constitution_key']].add(cid)
    coverage={t['compound_id'] for t in current['targets'] if t['net_status']=='exact-net-conversion-hypothesis'}
    targets=[]
    for t in current['targets']:
        if t['net_status']!='no-net-producing-equation':
            continue
        cid=t['compound_id']; candidates=[]
        for other in sorted(indexes.get(identity(compounds[cid]['smiles'])['constitution_key'],set())-{cid}):
            status=pair_status(compounds[cid]['smiles'],compounds[other]['smiles'])
            if status:
                candidates.append({'compound_id':other,'pair_status':status,
                    'partner_has_saved_conditional_CO2_certificate':other in coverage,
                    'partner_producing_step_ids':sorted(allowed[other])})
                used.add(other); selected.update(s.rsplit(':',1)[0] for s in allowed[other])
        if candidates:
            targets.append({**t,'source_identity_conflict':t['cannabisdb_id'] in conflicts,
                'identity':identity(compounds[cid]['smiles']),'alternatives':candidates,
                'next_action':'source-identity-resolution-first' if t['cannabisdb_id'] in conflicts else
                    'seek-exact-stereoisomerization-source' if any(p['pair_status']=='specified-stereoisomer-difference' for p in candidates) else
                    'resolve-unspecified-stereochemistry-before-proposing-reaction'})
            used.add(cid)
    for rid in selected:
        used.update(p['compound_id'] for side in ('left','right') for p in reactions[rid][side])
    return {'schema':'cannabis-carbon.phase1-stereochemistry-queue.v1','existing_reaction_pairs':existing,
        'targets':targets,'compounds':[compounds[c] for c in sorted(used)],
        'reactions':[reactions[r] for r in sorted(selected)],
        'summary':{'inventory_records':len(current['targets']),'balanced_model_equations':len(reactions),
            'existing_reaction_pair_counts':dict(Counter(r['pair_status'] for r in existing)),
            'no_producer_records_with_stereo_leads':len(targets),
            'next_action_counts':dict(Counter(t['next_action'] for t in targets)),
            'records_with_specified_partner_and_saved_route':sum(not t['source_identity_conflict'] and any(a['pair_status']=='specified-stereoisomer-difference' and a['partner_has_saved_conditional_CO2_certificate'] for a in t['alternatives']) for t in targets),
            'new_reactions':0,'coverage_gain_claimed':0},
        'rdkit_version':rdBase.rdkitVersion,
        'claim_boundary':'Diagnostic only. Exact stereoisomers remain separate identities. Removing stereo is a search key, not an equation. Pairs with any unassigned potential stereochemistry require identity review, not automatic stereoisomerization. Known full equations preserve sources and all substrates; paired-structure similarity alone establishes neither a reaction nor Cannabis activity. Protonation, isotopes, connectivity and bond orders are not normalized. This does not exhaust all coupled stereo transformations, enzyme substrate scope or enhanced-stereochemistry representations.'}


def run():
    RDLogger.DisableLog('rdApp.warning'); root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    parent_path=root/'phase1-selenium-forward-net.json'; current=read(parent_path)
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',parent_path,
           root/'phase1-no-producer-audit.json']+[root/n for n in current['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale model input')
    reactions,compounds,_=assemble_current(docs[0],docs[1],current,docs[4:])
    conflicts={t['cannabisdb_id'] for t in docs[3]['source_identity_conflicts']}
    report=build(current,compounds,reactions,conflicts)
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root/'phase1-stereochemistry-queue.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)
    print(json.dumps([{k:t[k] for k in ('cannabisdb_id','label','next_action')} for t in report['targets'] if t['next_action']=='seek-exact-stereoisomerization-source']),flush=True)


if __name__=='__main__':
    run()
