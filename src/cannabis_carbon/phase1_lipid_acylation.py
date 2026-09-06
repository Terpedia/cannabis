"""Exact PA/PC sn-2 acylation hypotheses instantiated from Terpedia Rhea rows."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdqueries, rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced
from .phase1_target_coverage import encoded_structure

RULES = {'PA': 'RHEA:19710', 'PC': 'RHEA:12938'}
BOUNDARY = ('Generic-reaction instantiation hypothesis, not a verified Cannabis '
    'reaction or enzyme assignment. Preserve target acyl chains, isotope labels, '
    'sn stereochemistry and unspecified double-bond geometry. Acyl-CoA and '
    'lysophospholipid supply remain unresolved. PA protonation is an explicit '
    'separate balanced hypothesis, never identity merging. Forward acylation '
    'follows the generic source record but is unverified for the exact substrates '
    'in Cannabis. No new CO2 pathway is claimed by participation alone.')


def canonical(mol):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def query_for(template):
    query = Chem.RWMol(template)
    for atom in template.GetAtoms():
        if atom.GetAtomicNum() == 0:
            query.ReplaceAtom(atom.GetIdx(), rdqueries.AtomNumEqualsQueryAtom(6))
    return query.GetMol()


def exact_scaffold_matches(mol, template, *, uniquify=True):
    """Only R groups may extend outside the source scaffold; no headgroup drift."""
    result = []
    for match in mol.GetSubstructMatches(query_for(template), useChirality=True, uniquify=uniquify):
        if any(a.GetDegree() != mol.GetAtomWithIdx(match[a.GetIdx()]).GetDegree()
               or a.GetFormalCharge() != mol.GetAtomWithIdx(match[a.GetIdx()]).GetFormalCharge()
               or a.GetIsotope() != mol.GetAtomWithIdx(match[a.GetIdx()]).GetIsotope()
               for a in template.GetAtoms() if a.GetAtomicNum()):
            continue
        core = {match[a.GetIdx()] for a in template.GetAtoms() if a.GetAtomicNum()}
        tails = set(range(mol.GetNumAtoms())) - core
        # This initial class covers concrete acyclic hydrocarbon acyl chains only.
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).IsInRing()
               for i in tails):
            continue
        seen = set(core)
        for a in template.GetAtoms():
            if a.GetAtomicNum():
                continue
            stack = [match[a.GetIdx()]]
            while stack:
                i = stack.pop()
                if i in seen:
                    continue
                seen.add(i)
                stack.extend(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in core)
        if len(seen) == mol.GetNumAtoms():
            result.append(match)
    return result


def sn2_bond(template):
    choices = []
    for atom in template.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() != 1:
            continue
        for oxygen in atom.GetNeighbors():
            if oxygen.GetAtomicNum() != 8:
                continue
            for carbon in oxygen.GetNeighbors():
                if carbon.GetIdx() != atom.GetIdx() and carbon.GetAtomicNum() == 6:
                    choices.append((carbon.GetIdx(), oxygen.GetIdx()))
    if len(choices) != 1:
        raise ValueError('Source must specify exactly one sn-2 ester bond')
    return choices[0]


def precursors(product, coa, carbon, oxygen):
    """Transfer the unchanged sn-2 acyl group to the exact source CoA thiol."""
    sulfur = [a.GetIdx() for a in coa.GetAtoms() if a.GetAtomicNum() == 16 and a.GetDegree() == 1 and a.GetTotalNumHs() == 1]
    if len(sulfur) != 1 or product.GetBondBetweenAtoms(carbon, oxygen).IsInRing():
        raise ValueError('Unambiguous CoA thiol and acyclic ester required')
    edit = Chem.RWMol(Chem.CombineMols(product, coa))
    edit.RemoveBond(carbon, oxygen)
    edit.AddBond(carbon, product.GetNumAtoms() + sulfur[0], Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    parts = Chem.GetMolFrags(edit, asMols=True)
    if len(parts) != 2:
        raise ValueError('Acyl transfer must give exactly two precursors')
    return parts


def build(network, catalog):
    sources = {kind: next(r for r in catalog if r['rule_id'] == rid) for kind, rid in RULES.items()}
    compounds = {c['id']: c for c in network['compounds']}
    baseline_ids = set(compounds)
    original_participants = {m['compound_id'] for r in network['reactions'] for side in ('left', 'right') for m in r[side]}
    equations, rows, used = {}, [], set()

    def compound(mol):
        smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid

    def equation(left, right, kind, source_id):
        sides = [[{'compound_id': c, 'coefficient': n} for c, n in sorted(Counter(s).items())] for s in (left, right)]
        if not balanced(sides, compounds):
            raise ValueError('Instantiated equation failed isotope/element/charge balance')
        rid = stable_id('lipid-acylation-hypothesis', [kind, sides])
        equations.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
            'hypothesis_type': kind, 'source_reaction_id': source_id,
            'direction_status': 'forward-generic-source-instantiation; exact-Cannabis-direction-unverified' if kind == 'sn2-acylation' else 'acid-base-equilibrium-hypothesis',
            'enzyme_evidence_ids': [], 'balance_status': 'independently-element-isotope-charge-balanced',
            'claim_boundary': BOUNDARY})
        return rid

    for kind, source in sources.items():
        left, right = source['reaction_smarts'].split('>>')
        precursor_templates = [Chem.MolFromSmiles(s) for s in left.split('.')]
        products = [Chem.MolFromSmiles(s) for s in right.split('.')]
        template = next(m for m in products if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
        coa = next(m for m in products if all(a.GetAtomicNum() != 0 for a in m.GetAtoms()))
        coa_id = compound(coa)
        cidx, oidx = sn2_bond(template)
        neutral = Chem.RWMol(template)
        if kind == 'PA':
            for a in neutral.GetAtoms():
                if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
                    a.SetFormalCharge(0); a.SetNoImplicit(False)
            Chem.SanitizeMol(neutral)
        for target in network['targets']:
            mol = Chem.MolFromSmiles(compounds[target['compound_id']]['smiles'])
            matches = exact_scaffold_matches(mol, template)
            bridge = None
            if not matches and kind == 'PA':
                neutral_matches = exact_scaffold_matches(mol, neutral)
                if neutral_matches:
                    edit = Chem.RWMol(mol)
                    for a in template.GetAtoms():
                        if a.GetFormalCharge() == -1:
                            oxygen = edit.GetAtomWithIdx(neutral_matches[0][a.GetIdx()])
                            oxygen.SetFormalCharge(-1); oxygen.SetNumExplicitHs(0)
                    Chem.SanitizeMol(edit)
                    mol = edit.GetMol()
                    matches = exact_scaffold_matches(mol, template)
                    charged_id = compound(mol)
                    proton = compound(Chem.MolFromSmiles('[H+]'))
                    used.add(target['compound_id'])
                    bridge = equation([charged_id, proton, proton], [target['compound_id']], 'explicit-PA-protonation', source['rule_id'])
            if not matches:
                continue
            match = matches[0]
            product_id = compound(mol)
            parts = precursors(mol, coa, match[cidx], match[oidx])
            if len(parts) != len(precursor_templates) or any(
                sum(bool(exact_scaffold_matches(part, source_part)) for part in parts) != 1
                for source_part in precursor_templates):
                raise ValueError('Generated precursors fail exact source scaffold/stereo checks')
            precursor_ids = [compound(m) for m in parts]
            rid = equation(precursor_ids, [product_id, coa_id], 'sn2-acylation', source['rule_id'])
            rows.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id')} | {
                'lipid_class': kind, 'source_smiles': target['source_smiles'],
                'encoded_structure_status': encoded_structure(target['source_smiles'])[1],
                'baseline_balanced_participant': target['compound_id'] in original_participants,
                'acylation_reaction_id': rid, 'protonation_reaction_id': bridge,
                'reaction_product_id': product_id, 'required_precursor_ids': precursor_ids,
                'precursors_absent_from_baseline_network': sorted(set(precursor_ids) - baseline_ids),
                'status': 'balanced-producing-reaction-hypothesis; upstream-supply-unresolved'})
    return {'schema': 'cannabis-carbon.phase1-lipid-acylation.v1', 'claim_boundary': BOUNDARY,
        'source_records': list(sources.values()), 'targets': rows, 'reactions': list(equations.values()),
        'compounds': [compounds[c] for c in sorted(used)], 'rdkit_version': rdBase.rdkitVersion,
        'summary': {'target_records': len(rows), 'lipid_class_counts': dict(Counter(t['lipid_class'] for t in rows)),
            'balanced_equations': len(equations), 'reaction_type_counts': dict(Counter(r['hypothesis_type'] for r in equations.values())),
            'targets_without_baseline_balanced_participation': sum(not t['baseline_balanced_participant'] for t in rows),
            'targets_with_precursors_absent_from_baseline': sum(bool(t['precursors_absent_from_baseline_network']) for t in rows),
            'new_CO2_pathway_claims': 0},
        'scope': 'Full CannabisDB inventory screened by exact source scaffold, not label prefixes. Only PA and PC with acyclic hydrocarbon acyl tails and source-compatible sn configuration; other headgroups, ether lipids and unspecified sn configuration are outside this batch.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-full-balanced-network.json'), Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-lipid-acylation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
