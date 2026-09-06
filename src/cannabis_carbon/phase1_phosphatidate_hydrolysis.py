"""Source-forward phosphatidate hydrolysis hypotheses for exact DAG structures."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches
from .phase1_marts_completions import balanced

BOUNDARY = ('Exact structural instantiation of generic RHEA:27430, not demonstrated '
    'Cannabis activity or exact-substrate specificity. Source-forward hydrolysis '
    'requires explicit phosphatidate dianion and water and produces exact DAG '
    'and hydrogen phosphate dianion. Acyl-chain identity, isotopes and encoded '
    'stereochemistry are preserved; unassigned glycerol stereo is not completed. '
    'Neutral phosphatidic acid is not merged with the dianion. No precursor '
    'availability, CO2 conversion, energy feasibility or enzyme assignment is '
    'established by a balanced equation alone.')


def source_parts(source):
    if source['rule_id'] != 'RHEA:27430':
        raise ValueError('Require source-forward hydrolysis record')
    left, right = [[Chem.MolFromSmiles(s) for s in side.split('.')]
                   for side in source['reaction_smarts'].split('>>')]
    def generic(m):
        return any(a.GetAtomicNum() == 0 for a in m.GetAtoms())
    pa = next(m for m in left if generic(m))
    dag = next(m for m in right if generic(m))
    water = next(m for m in left if not generic(m))
    phosphate = next(m for m in right if not generic(m))
    if canonical(water) != 'O' or canonical(phosphate) != 'O=P([O-])([O-])O':
        raise ValueError('Unexpected source water/phosphate identity')
    if Chem.GetFormalCharge(pa) != -2 or Chem.GetFormalCharge(dag) != 0:
        raise ValueError('Unexpected source lipid charges')
    return pa, dag, water, phosphate


def instantiate(dag, source):
    pa_template, dag_template, water, phosphate = source_parts(source)
    if not exact_scaffold_matches(dag, dag_template):
        return []
    oh = [a.GetIdx() for a in dag.GetAtoms()
          if a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 1]
    phosphorus = [a.GetIdx() for a in phosphate.GetAtoms() if a.GetAtomicNum() == 15]
    leaving = [a.GetIdx() for a in phosphate.GetAtoms()
               if a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 1]
    if len(oh) != 1 or len(phosphorus) != 1 or len(leaving) != 1:
        raise ValueError('Require unique DAG hydroxyl and phosphate hydroxyl')
    # Reverse construction proposes a precursor; the admitted reaction is hydrolysis.
    offset = dag.GetNumAtoms()
    edit = Chem.RWMol(Chem.CombineMols(dag, phosphate))
    edit.RemoveBond(offset + phosphorus[0], offset + leaving[0])
    edit.AddBond(oh[0], offset + phosphorus[0], Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    fragments = Chem.GetMolFrags(edit, asMols=True)
    pa = next(m for m in fragments if m.GetNumAtoms() > 1)
    released = next(m for m in fragments if m.GetNumAtoms() == 1)
    if canonical(released) != canonical(water) or not exact_scaffold_matches(pa, pa_template):
        raise ValueError('Generated precursor fails exact source scaffold/stereo')
    return [{'reactant_smiles': sorted([canonical(pa), canonical(water)]),
             'product_smiles': sorted([canonical(dag), canonical(phosphate)])}]


def build(parent, catalog):
    source = next(r for r in catalog if r['rule_id'] == 'RHEA:27430')
    source_parts(source)
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    original = set(compounds)
    target_ids = {}
    for t in parent['targets']:
        target_ids.setdefault(t['compound_id'], []).append(t['cannabisdb_id'])
    reactions, rows, used = {}, [], set()

    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles)
        smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid

    # Scan the complete parent compound inventory, including required intermediates.
    for cid in sorted(original):
        mol = Chem.MolFromSmiles(compounds[cid]['smiles'])
        for candidate in instantiate(mol, source):
            sides = [[{'compound_id': c, 'coefficient': n}
                      for c, n in sorted(Counter(participant(s) for s in candidate[side]).items())]
                     for side in ('reactant_smiles', 'product_smiles')]
            if not balanced(sides, compounds):
                raise ValueError('Phosphatidate hydrolysis is not element/isotope/charge balanced')
            rid = stable_id('phosphatidate-hydrolysis-hypothesis', sides)
            reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                'hypothesis_type': 'phosphatidate-hydrolysis',
                'source_reaction_id': source['rule_id'], 'source_url': source['source_url'],
                'direction_status': 'forward-only-generic-source-hypothesis',
                'balance_status': 'independently-element-isotope-charge-balanced',
                'enzyme_evidence_ids': [], 'claim_boundary': BOUNDARY})
            required = {p['compound_id'] for p in sides[0]}
            rows.append({'compound_id': cid, 'cannabisdb_ids': target_ids.get(cid, []),
                'reaction_id': rid, 'required_precursor_ids': sorted(required),
                'precursors_absent_from_parent': sorted(required - original)})
    return {'schema': 'cannabis-carbon.phase1-phosphatidate-hydrolysis.v1',
        'claim_boundary': BOUNDARY, 'source_record': source, 'dag_candidates': rows,
        'reactions': list(reactions.values()), 'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'parent_compound_structures_scanned': len(original),
            'dag_structures_matched': len({r['compound_id'] for r in rows}),
            'direct_target_records': len({t for r in rows for t in r['cannabisdb_ids']}),
            'balanced_equations': len(reactions),
            'candidates_with_absent_precursors': sum(bool(r['precursors_absent_from_parent']) for r in rows),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-triglyceride-symmetry-net.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-phosphatidate-hydrolysis.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
