"""Exact PG/PGP speciation plus recursively instantiated precursor reactions."""
import hashlib
import json
from collections import Counter, deque
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical
from .phase1_marts_completions import balanced
from .phase1_cardiolipin_precursors import instantiate as phospho, RULES
from .phase1_glycerolipid_precursors import instantiate as acyl

BOUNDARY = ('Generic-source chemistry hypothesis, not confirmed Cannabis biology. '
    'All substrates and coproducts are explicit, with exact identities and encoded '
    'stereochemistry retained. Synthesis/hydrolysis is source-forward only; separate '
    'proton exchanges are explicitly assumed reversible, not curated reactions. '
    'No organic input is seeded. Inventory presence is not proof of supply; this '
    'proposal report establishes neither CO2 conversion nor pool startup or flux.')


def build(parent, audit, catalog):
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    known = set(compounds)
    for c in audit['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Compound identity conflict')
        compounds.setdefault(c['id'], dict(c))
    rules = RULES | {'sn2-acylation': 'RHEA:19710', 'sn1-acylation': 'RHEA:15326'}
    sources = {k:next(s for s in catalog if s['rule_id'] == rid) for k,rid in rules.items()}
    reactions, rows, targets = {}, [], []
    used = set()
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id':cid, 'smiles':smiles,
            'formula':rdMolDescriptors.CalcMolFormula(mol), 'formal_charge':Chem.GetFormalCharge(mol),
            'carbon_count':sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    for t in audit['targets']:
        if not balanced([t['left'], t['right']], compounds):
            raise ValueError('Unbalanced speciation')
        rid = stable_id('glycerophospholipid-speciation', [t['left'], t['right']])
        reactions.setdefault(rid, {'id':rid, 'left':t['left'], 'right':t['right'],
            'hypothesis_type':'glycerophospholipid-speciation', 'speciation_type':t['speciation_type'],
            'source_reaction_id':t['source_reaction_id'], 'source_url':t['source_url'],
            'direction_status':'explicit-reversible-speciation-assumption-not-source-reaction',
            'enzyme_evidence_ids':[], 'balance_status':'independently-element-isotope-charge-balanced',
            'claim_boundary':t['claim_boundary'] + ' Reversibility is an explicit sensitivity assumption.'})
        used.update(p['compound_id'] for side in (t['left'], t['right']) for p in side)
        targets.append({k:t[k] for k in ('cannabisdb_id','compound_id','label','lipid_class','net_status')} | {
            'source_form_compound_id':t['source_form_compound_id'], 'speciation_reaction_id':rid})
    roots = {t['source_form_compound_id'] for t in targets}
    queue = deque(sorted(roots)); visited = set()
    while queue:
        cid = queue.popleft()
        if cid in visited:
            continue
        if len(visited) >= 2000:
            raise ValueError('Explicit 2000-structure expansion bound exceeded; no partial report')
        visited.add(cid)
        for kind, source in sources.items():
            mol = Chem.MolFromSmiles(compounds[cid]['smiles'])
            for candidate in (phospho if kind in RULES else acyl)(mol, source):
                sides = [[{'compound_id':p, 'coefficient':n} for p,n in
                          sorted(Counter(participant(s) for s in candidate[side]).items())]
                         for side in ('reactant_smiles','product_smiles')]
                if not balanced(sides, compounds):
                    raise ValueError('Unbalanced synthesis')
                rid = stable_id('glycerophospholipid-synthesis', [kind,sides])
                reactions.setdefault(rid, {'id':rid, 'left':sides[0], 'right':sides[1],
                    'hypothesis_type':kind, 'source_reaction_id':source['rule_id'], 'source_url':source['source_url'],
                    'direction_status':'forward-only-generic-source-hypothesis', 'enzyme_evidence_ids':[],
                    'balance_status':'independently-element-isotope-charge-balanced', 'claim_boundary':BOUNDARY})
                required = sorted(p['compound_id'] for p in sides[0])
                rows.append({'compound_id':cid, 'reaction_id':rid, 'required_precursor_ids':required})
                queue.extend(p for p in required if p not in visited)
    proposed = {r['compound_id'] for r in rows}
    frontier = visited - proposed
    return {'schema':'cannabis-carbon.phase1-glycerophospholipid-synthesis.v1',
        'source_records':list(sources.values()), 'targets':targets, 'reactions':list(reactions.values()),
        'compounds':[compounds[c] for c in sorted(used | visited)], 'precursor_candidates':rows,
        'frontier':[{'compound_id':c, 'in_parent_inventory':c in known,
            'status':'no-producer-from-this-template-set; parent-pathway-supply-not-assumed'} for c in sorted(frontier)],
        'claim_boundary':BOUNDARY,
        'summary':{'target_records':len(targets), 'root_structures':len(roots),
            'roots_without_synthesis_hypotheses':len(roots-proposed), 'structures_examined':len(visited),
            'balanced_equations':len(reactions), 'equations_by_type':dict(Counter(r['hypothesis_type'] for r in reactions.values())),
            'frontier_structures':len(frontier), 'frontier_absent_from_parent_inventory':len(frontier-known),
            'new_CO2_route_claims':0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-amino-phospholipid-net.json'),
             Path('data/reports/phase1-glycerophospholipid-speciation.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs[:2]:
        for p,sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-glycerophospholipid-synthesis.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
