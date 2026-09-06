"""Explicit precursor stereochemistry for symmetric, achiral triglycerides."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches, precursors
from .phase1_triglyceride_acylation import source_parts
from .phase1_marts_completions import balanced
from .phase1_target_coverage import encoded_structure

BOUNDARY = ('Source-based stereochemical instantiation hypothesis, not demonstrated '
    'Cannabis activity. The encoded target glycerol center is achiral, not unknown. Unresolved acyl-chain stereochemistry remains unresolved. Only '
    'the newly asymmetric diacylglycerol precursor receives the configuration '
    'required by the source scaffold. Target identity, acyl-chain geometry, '
    'isotopes and all other stereocenters are unchanged. Unknown target glycerol '
    'stereochemistry is excluded. Forward acylation requires both explicit '
    'precursors; balance alone does not establish their availability or a CO2 route.')


def analyze(mol, source):
    template, coa, inputs, (ci, oi) = source_parts(source)
    centers = [a.GetIdx() for a in template.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
    if len(centers) != 1:
        raise ValueError('One source glycerol center required')
    relaxed = Chem.Mol(template); Chem.RemoveStereochemistry(relaxed)
    matches = exact_scaffold_matches(mol, relaxed)
    if not matches:
        return {'status': 'outside-source-scaffold', 'candidates': []}
    possible = {int(s.centeredOn) for s in Chem.FindPotentialStereo(mol) if str(s.type) == 'Atom_Tetrahedral'}
    candidate_rows = {}
    for match in matches:
        center = match[centers[0]]
        if center in possible:
            return {'status': 'glycerol-center-is-stereogenic; separate-identity-review', 'candidates': []}
        product = Chem.Mol(mol)
        for a in product.GetAtoms():
            if a.HasProp('_dgat_center'):
                a.ClearProp('_dgat_center')
        product.GetAtomWithIdx(center).SetBoolProp('_dgat_center', True)
        parts = precursors(product, coa, match[ci], match[oi])
        dag = next(m for m in parts if any(a.HasProp('_dgat_center') for a in m.GetAtoms()))
        donor = next(m for m in parts if m is not dag)
        dag_center = next(a.GetIdx() for a in dag.GetAtoms() if a.HasProp('_dgat_center'))
        for tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            variant = Chem.Mol(dag)
            variant.GetAtomWithIdx(dag_center).SetChiralTag(tag)
            Chem.AssignStereochemistry(variant, cleanIt=True, force=True)
            pair = [variant, donor]
            if any(sum(bool(exact_scaffold_matches(part, query)) for part in pair) != 1 for query in inputs):
                continue
            smiles = tuple(sorted(canonical(p) for p in pair))
            candidate_rows.setdefault(smiles, {'reactant_smiles': list(smiles),
                'product_smiles': sorted([canonical(mol), canonical(coa)]),
                'stereo_assumption': 'Source-required sn configuration introduced only in the newly asymmetric DAG precursor; no target stereo assignment.'})
    if not candidate_rows:
        raise ValueError('Achiral target has no source-compatible stereospecific precursor')
    return {'status': 'symmetric-achiral-target; source-stereospecific-precursor-hypothesis', 'candidates': list(candidate_rows.values())}


def build(review, parent):
    compounds = {c['id']: c for c in parent['compounds']}
    original_ids = set(compounds)
    targets, equations, used = [], {}, set()
    source = review['source_record']
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles, 'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol), 'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    for t in review['review_targets']:
        mol = Chem.MolFromSmiles(t['source_smiles'])
        if stable_id('structure', canonical(mol)) != t['compound_id']:
            raise ValueError('Source target identity mismatch')
        result = analyze(mol, source)
        ids, required = [], set()
        for candidate in result['candidates']:
            sides = []
            for side in ('reactant_smiles', 'product_smiles'):
                counts = Counter(participant(s) for s in candidate[side])
                sides.append([{'compound_id': c, 'coefficient': n} for c, n in sorted(counts.items())])
            if not balanced(sides, compounds):
                raise ValueError('Symmetric triglyceride equation is not balanced')
            rid = stable_id('triglyceride-symmetry-hypothesis', sides)
            equations.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                'hypothesis_type': 'sn3-acylation', 'source_reaction_id': source['rule_id'],
                'source_url': source['source_url'], 'enzyme_evidence_ids': [],
                'direction_status': 'forward-only-generic-source-hypothesis',
                'stereo_assumption': candidate['stereo_assumption'], 'claim_boundary': BOUNDARY})
            ids.append(rid); required.update(p['compound_id'] for p in sides[0])
        targets.append({**t, 'status': result['status'], 'encoded_structure_status': encoded_structure(t['source_smiles'])[1], 'reaction_ids': sorted(set(ids)),
            'required_precursor_ids': sorted(required), 'precursors_absent_from_parent': sorted(required - original_ids)})
    return {'schema': 'cannabis-carbon.phase1-triglyceride-symmetry.v1', 'claim_boundary': BOUNDARY,
        'source_record': source, 'targets': targets, 'reactions': list(equations.values()),
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'reviewed_target_records': len(targets), 'status_counts': dict(Counter(t['status'] for t in targets)),
            'balanced_equations': len(equations), 'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-triglyceride-acylation.json'), Path('data/reports/phase1-triglyceride-net.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-triglyceride-symmetry.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
