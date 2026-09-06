"""Source-forward DGAT reaction hypotheses with explicit sn-3 acyl transfer."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches, precursors
from .phase1_marts_completions import balanced
from .phase1_target_coverage import encoded_structure

BOUNDARY = ('Generic RHEA:10869 instantiation only, not a validated Cannabis '
    'reaction or enzyme assignment. Exact acyl chains, encoded stereochemistry '
    'and isotope labels are retained. Source-forward sn-3 acylation requires '
    'both the explicit 1,2-diacyl-sn-glycerol and acyl-CoA; neither is an external '
    'carbon seed. Precursor supply and complete CO2 conversion are not established '
    'by this report. No target stereo is inferred from its name. Symmetric or '
    'unassigned glycerol configurations that do not match the source scaffold '
    'remain a separate review set; they are not claimed chemically impossible.')


def source_parts(source):
    left, right = source['reaction_smarts'].split('>>')
    inputs = [Chem.MolFromSmiles(s) for s in left.split('.')]
    product, coa = [Chem.MolFromSmiles(s) for s in right.split('.')]
    third = [a for a in product.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == 3]
    if len(third) != 1 or third[0].GetDegree() != 1:
        raise ValueError('Expected one source sn-3 acyl substituent')
    carbon = third[0].GetNeighbors()[0]
    oxygen = [a for a in carbon.GetNeighbors() if a.GetAtomicNum() == 8 and a.GetDegree() == 2]
    if len(oxygen) != 1:
        raise ValueError('Expected source sn-3 ester oxygen')
    # Verify this is a terminal glycerol O-acyl group, not the central sn-2 group.
    glycerol = [a for a in oxygen[0].GetNeighbors() if a.GetIdx() != carbon.GetIdx()][0]
    if glycerol.GetAtomicNum() != 6 or glycerol.GetTotalNumHs() != 2:
        raise ValueError('Source acyl transfer must be on terminal glycerol carbon')
    return product, coa, inputs, (carbon.GetIdx(), oxygen[0].GetIdx())


def instantiate(mol, source):
    template, coa, input_templates, (ci, oi) = source_parts(source)
    candidates = {}
    for match in exact_scaffold_matches(mol, template):
        parts = precursors(mol, coa, match[ci], match[oi])
        if any(sum(bool(exact_scaffold_matches(part, query)) for part in parts) != 1
               for query in input_templates):
            raise ValueError('Generated DGAT inputs fail source stereo/scaffold validation')
        smiles = tuple(sorted(canonical(part) for part in parts))
        candidates.setdefault(smiles, {'reactant_smiles': list(smiles),
            'product_smiles': sorted([canonical(mol), canonical(coa)]),
            'transferred_target_bond': [match[ci], match[oi]]})
    return list(candidates.values())


def build(network, parent, catalog):
    source = next(r for r in catalog if r['rule_id'] == 'RHEA:10869')
    source_parts(source)
    compounds = {c['id']: c for c in parent['compounds']}
    parent_ids = set(compounds)
    statuses = {t['cannabisdb_id']: t for t in parent['targets']}
    used, equations, targets, reviews = set(), {}, [], []

    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid

    for t in network['targets']:
        before = statuses[t['cannabisdb_id']]
        if before['compound_id'] != t['compound_id']:
            raise ValueError('Target inventory identity mismatch')
        mol = Chem.MolFromSmiles(compounds[t['compound_id']]['smiles'])
        candidates = instantiate(mol, source)
        if not candidates:
            if t['label'].startswith('TG('):
                reviews.append({k: t[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'source_smiles')} | {
                    'status': 'source-scaffold-or-stereo-review-required',
                    'next_action': 'Distinguish symmetric achiral triacylglycerols from unspecified sn configuration and out-of-scope chains; do not assign stereochemistry by name.'})
            continue
        reaction_ids, precursor_ids = [], set()
        for candidate in candidates:
            sides = []
            for side in ('reactant_smiles', 'product_smiles'):
                counts = Counter(participant(s) for s in candidate[side])
                sides.append([{'compound_id': cid, 'coefficient': n} for cid, n in sorted(counts.items())])
            if not balanced(sides, compounds):
                raise ValueError('DGAT equation failed element/isotope/charge balance')
            rid = stable_id('triglyceride-acylation-hypothesis', sides)
            equations.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                'hypothesis_type': 'sn3-acylation', 'source_reaction_id': source['rule_id'],
                'source_url': source['source_url'], 'enzyme_evidence_ids': [],
                'direction_status': 'forward-only-generic-source-hypothesis',
                'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
            reaction_ids.append(rid)
            precursor_ids.update(p['compound_id'] for p in sides[0])
        targets.append({k: t[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'source_smiles')} | {
            'encoded_structure_status': encoded_structure(t['source_smiles'])[1],
            'reaction_ids': sorted(set(reaction_ids)), 'required_precursor_ids': sorted(precursor_ids),
            'precursors_absent_from_parent_network': sorted(precursor_ids - parent_ids),
            'parent_net_status': before['net_status'], 'parent_balanced_participant': before['balanced_participant'],
            'status': 'balanced-producing-reaction-hypotheses; upstream-supply-unresolved'})
    return {'schema': 'cannabis-carbon.phase1-triglyceride-acylation.v1',
        'claim_boundary': BOUNDARY, 'source_record': source, 'rdkit_version': rdBase.rdkitVersion,
        'targets': targets, 'review_targets': reviews, 'reactions': list(equations.values()),
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'inventory_records_screened': len(network['targets']), 'target_records_with_hypotheses': len(targets),
            'unique_target_structures': len({t['compound_id'] for t in targets}),
            'balanced_equations': len(equations),
            'previous_no_producer_records_with_hypotheses': sum(t['parent_net_status'] == 'no-net-producing-equation' for t in targets),
            'previous_nonparticipants_with_hypotheses': sum(not t['parent_balanced_participant'] for t in targets),
            'targets_with_precursors_absent_from_parent': sum(bool(t['precursors_absent_from_parent_network']) for t in targets),
            'TG_label_records_requiring_review': len(reviews), 'new_CO2_pathway_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-full-balanced-network.json'),
             Path('data/reports/phase1-lipid-acylation-net.json'), Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-triglyceride-acylation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
