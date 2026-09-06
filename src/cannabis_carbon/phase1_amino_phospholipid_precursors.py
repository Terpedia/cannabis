"""Backward expansion of every required amino-phospholipid synthesis input."""
import hashlib
import json
from collections import Counter, deque
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical
from .phase1_marts_completions import balanced
from .phase1_amino_phospholipid_synthesis import instantiate as amino, KINDS, BOUNDARY
from .phase1_amino_phospholipid_speciation import RULES as AMINO_RULES
from .phase1_cardiolipin_precursors import instantiate as cdp
from .phase1_glycerolipid_precursors import instantiate as acyl


def build(parent, synthesis, catalog):
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    known = set(compounds)
    for c in synthesis['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Input compound identity conflict')
        compounds.setdefault(c['id'], dict(c))
    rules = {KINDS[k]: v for k, v in AMINO_RULES.items()} | {
        'cdp-dag-synthesis': 'RHEA:16230', 'sn2-acylation': 'RHEA:19710', 'sn1-acylation': 'RHEA:15326'}
    sources = {k: next(s for s in catalog if s['rule_id'] == rid) for k, rid in rules.items()}
    roots = {p for t in synthesis['targets'] for p in t['required_precursor_ids']}
    queue = deque(sorted(roots)); visited = set(); used = set(); reactions = {}; rows = []
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol); cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles, 'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol), 'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    while queue:
        cid = queue.popleft()
        if cid in visited:
            continue
        if len(visited) >= 2000:
            raise ValueError('Expansion exceeded explicit 2000-structure audit bound; no partial report')
        visited.add(cid)
        mol = Chem.MolFromSmiles(compounds[cid]['smiles'])
        for kind, source in sources.items():
            candidates = (amino(mol, source, kind) if kind in KINDS.values() else
                          cdp(mol, source) if kind == 'cdp-dag-synthesis' else acyl(mol, source))
            for candidate in candidates:
                sides = [[{'compound_id': p, 'coefficient': n} for p, n in
                          sorted(Counter(participant(s) for s in candidate[side]).items())]
                         for side in ('reactant_smiles', 'product_smiles')]
                if not balanced(sides, compounds):
                    raise ValueError('Precursor equation fails full balance')
                rid = stable_id('amino-phospholipid-precursor', [kind, sides])
                reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                    'hypothesis_type': kind, 'source_reaction_id': source['rule_id'], 'source_url': source['source_url'],
                    'direction_status': 'forward-only-generic-source-hypothesis', 'enzyme_evidence_ids': [],
                    'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
                required = sorted(p['compound_id'] for p in sides[0])
                queue.extend(p for p in required if p not in visited)
                rows.append({'compound_id': cid, 'reaction_id': rid, 'required_precursor_ids': required})
    proposed = {r['compound_id'] for r in rows}
    absent_roots = roots - known
    # Inventory presence is recorded, never treated as proof of producing reactions.
    frontier = visited - proposed
    return {'schema': 'cannabis-carbon.phase1-amino-phospholipid-precursors.v1',
        'source_records': list(sources.values()), 'reactions': list(reactions.values()),
        'compounds': [compounds[c] for c in sorted(used | visited)], 'precursor_candidates': rows,
        'frontier': [{'compound_id': c, 'in_parent_inventory': c in known,
                      'status': 'no-producer-from-this-template-set; parent-pathway-supply-not-assumed'} for c in sorted(frontier)],
        'claim_boundary': BOUNDARY + ' Backward expansion stops where this explicit template set has no match; parent inventory presence is not a supply certificate.',
        'summary': {'root_precursor_structures': len(roots), 'structures_examined': len(visited),
            'balanced_equations': len(reactions), 'equations_by_type': dict(Counter(r['hypothesis_type'] for r in reactions.values())),
            'missing_root_precursors_with_producing_hypotheses': len(absent_roots & proposed),
            'missing_root_precursors_without_producing_hypotheses': len(absent_roots - proposed),
            'frontier_structures': len(frontier), 'frontier_absent_from_parent_inventory': len(frontier - known),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-source-mapped-protonation-net.json'),
             Path('data/reports/phase1-amino-phospholipid-synthesis.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs[:2]:
        for p, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-amino-phospholipid-precursors.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
