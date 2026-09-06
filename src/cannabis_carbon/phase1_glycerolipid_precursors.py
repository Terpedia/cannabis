"""Exact PA/LPA precursor synthesis, including newly required intermediates."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches, precursors, sn2_bond
from .phase1_marts_completions import balanced

RULES = {'sn2-acylation': 'RHEA:19710', 'sn1-acylation': 'RHEA:15326'}
BOUNDARY = ('Source-forward generic acyltransferase instantiation hypotheses, not '
    'demonstrated Cannabis activity or exact-substrate specificity. Explicit '
    'acyl-CoA and lysophosphatidate/sn-glycerol-3-phosphate inputs are required. '
    'No organic precursor is seeded externally. Exact charge, sn configuration, '
    'isotopes and encoded chain geometry are retained, without identity merging '
    'or completion of unassigned stereo. These balanced reactions alone do not '
    'establish CO2 conversion, pool startup, energetics or compartment feasibility.')


def instantiate(mol, source):
    if source['rule_id'] not in RULES.values():
        raise ValueError('Require a supported source-forward acyltransferase')
    left, right = [[Chem.MolFromSmiles(s) for s in side.split('.')]
                   for side in source['reaction_smarts'].split('>>')]
    template = next(m for m in right if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
    coa = next(m for m in right if all(a.GetAtomicNum() != 0 for a in m.GetAtoms()))
    if source['rule_id'] == RULES['sn2-acylation']:
        carbon, oxygen = sn2_bond(template)
    else:
        dummy = [a for a in template.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummy) != 1 or dummy[0].GetDegree() != 1:
            raise ValueError('Require one source sn-1 acyl substituent')
        carbon = dummy[0].GetNeighbors()[0].GetIdx()
        ester = [a for a in template.GetAtomWithIdx(carbon).GetNeighbors()
                 if a.GetAtomicNum() == 8 and a.GetDegree() == 2]
        if len(ester) != 1:
            raise ValueError('Require one source ester oxygen')
        oxygen = ester[0].GetIdx()
        glycerol = next(a for a in ester[0].GetNeighbors() if a.GetIdx() != carbon)
        if glycerol.GetAtomicNum() != 6 or glycerol.GetTotalNumHs() != 2:
            raise ValueError('Source sn-1 acylation must act on terminal glycerol')
    candidates = {}
    for match in exact_scaffold_matches(mol, template):
        pair = precursors(mol, coa, match[carbon], match[oxygen])
        if any(sum(bool(exact_scaffold_matches(p, q)) for p in pair) != 1 for q in left):
            raise ValueError('Generated inputs fail source stereochemistry/scaffold')
        key = tuple(sorted(canonical(p) for p in pair))
        candidates[key] = {'reactant_smiles': list(key),
                           'product_smiles': sorted([canonical(mol), canonical(coa)])}
    return list(candidates.values())


def build(parent, hydrolysis, catalog):
    sources = {k: next(r for r in catalog if r['rule_id'] == rid) for k, rid in RULES.items()}
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    for c in hydrolysis['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Hydrolysis compound identity conflict')
        compounds.setdefault(c['id'], dict(c))
    original = set(compounds)
    used, equations, rows = set(), {}, []
    targets = {}
    for t in parent['targets']:
        targets.setdefault(t['compound_id'], []).append(t['cannabisdb_id'])

    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid

    # First construct LPA inputs to every exact PA, then G3P inputs to every LPA.
    # Include intermediates, not just CannabisDB-labelled target records.
    for kind, source in sources.items():
        for cid in sorted(compounds):
            for candidate in instantiate(Chem.MolFromSmiles(compounds[cid]['smiles']), source):
                sides = [[{'compound_id': c, 'coefficient': n}
                          for c, n in sorted(Counter(participant(s) for s in candidate[side]).items())]
                         for side in ('reactant_smiles', 'product_smiles')]
                if not balanced(sides, compounds):
                    raise ValueError('Precursor acylation failed element/isotope/charge balance')
                rid = stable_id('glycerolipid-precursor-hypothesis', [kind, sides])
                equations.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                    'hypothesis_type': kind, 'source_reaction_id': source['rule_id'],
                    'source_url': source['source_url'], 'enzyme_evidence_ids': [],
                    'direction_status': 'forward-only-generic-source-hypothesis',
                    'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
                required = {p['compound_id'] for p in sides[0]}
                rows.append({'compound_id': cid, 'cannabisdb_ids': targets.get(cid, []),
                    'reaction_id': rid, 'hypothesis_type': kind,
                    'required_precursor_ids': sorted(required),
                    'precursors_absent_from_input_inventory': sorted(required - original)})
    missing_pa = {c for r in hydrolysis['dag_candidates'] for c in r['precursors_absent_from_parent']}
    proposed_pa = {r['compound_id'] for r in rows if r['hypothesis_type'] == 'sn2-acylation'}
    return {'schema': 'cannabis-carbon.phase1-glycerolipid-precursors.v1',
        'claim_boundary': BOUNDARY, 'source_records': list(sources.values()),
        'precursor_candidates': rows, 'reactions': list(equations.values()),
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'input_compound_structures': len(original),
            'balanced_equations': len(equations),
            'equations_by_type': dict(Counter(r['hypothesis_type'] for r in equations.values())),
            'hydrolysis_missing_precursors_with_synthesis_hypotheses': len(missing_pa & proposed_pa),
            'hydrolysis_missing_precursors_without_synthesis_hypotheses': len(missing_pa - proposed_pa),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-triglyceride-symmetry-net.json'),
             Path('data/reports/phase1-phosphatidate-hydrolysis.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-glycerolipid-precursors.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
