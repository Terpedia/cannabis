"""Exact source-forward cardiolipin synthesis with separate protonation bridges."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches
from .phase1_marts_completions import balanced

BOUNDARY = ('Generic RHEA:32932 instantiation, not demonstrated Cannabis activity '
    'or exact-substrate specificity. Explicit CDP-diacylglycerol and phosphatidylglycerol '
    'inputs retain their source-required stereochemistry. CMP and H+ are explicit '
    'coproducts. Neutral cardiolipin is a separate identity reached only through an '
    'explicit balanced protonation hypothesis. Target central stereochemistry is '
    'never completed from names. Balance alone does not establish precursor supply, '
    'CO2 conversion, pool startup, energetics or compartments.')


def source_parts(source):
    if source['rule_id'] != 'RHEA:32932':
        raise ValueError('Require source-forward cardiolipin synthase record')
    left, right = [[Chem.MolFromSmiles(s) for s in side.split('.')]
                   for side in source['reaction_smarts'].split('>>')]
    template = next(m for m in right if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
    cmp = next(m for m in right if any(a.GetAtomicNum() == 7 for a in m.GetAtoms()))
    proton = next(m for m in right if m.GetNumAtoms() == 1)
    centers = [a for a in template.GetAtoms() if a.GetAtomicNum() == 6
               and any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1 for n in a.GetNeighbors())]
    if len(centers) != 1 or canonical(proton) != '[H+]':
        raise ValueError('Unexpected source central glycerol or proton')
    choices = []
    for carbon in centers[0].GetNeighbors():
        if carbon.GetAtomicNum() != 6:
            continue
        for oxygen in carbon.GetNeighbors():
            if oxygen.GetAtomicNum() != 8:
                continue
            for phosphorus in oxygen.GetNeighbors():
                if phosphorus.GetAtomicNum() != 15:
                    continue
                cut = Chem.RWMol(template)
                cut.RemoveBond(phosphorus.GetIdx(), oxygen.GetIdx())
                fragment = next(f for f in Chem.GetMolFrags(cut) if phosphorus.GetIdx() in f)
                if any(template.GetAtomWithIdx(i).GetAtomicNum() == 0
                       and template.GetAtomWithIdx(i).GetIsotope() == 1 for i in fragment):
                    choices.append((phosphorus.GetIdx(), oxygen.GetIdx()))
    if len(choices) != 1:
        raise ValueError('Require one source donor-phosphatidyl linkage')
    return template, cmp, proton, left, centers[0].GetIdx(), choices[0]


def instantiate(mol, source):
    template, cmp, proton, inputs, center, (pidx, oidx) = source_parts(source)
    neutral = Chem.RWMol(template)
    for a in neutral.GetAtoms():
        if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
            a.SetFormalCharge(0); a.SetNoImplicit(False)
    Chem.SanitizeMol(neutral)
    matches = exact_scaffold_matches(mol, template, uniquify=False)
    needs_bridge = not matches
    if needs_bridge:
        matches = exact_scaffold_matches(mol, neutral, uniquify=False)
    if not matches:
        return {'status': 'outside-source-scaffold', 'candidates': []}
    potential = {int(s.centeredOn): str(s.specified) for s in Chem.FindPotentialStereo(mol)}
    if any(potential.get(m[center], 'achiral') != 'Specified' for m in matches):
        return {'status': 'central-glycerol-stereo-review-required', 'candidates': []}
    candidates = {}
    for match in matches:
        charged = Chem.Mol(mol)
        if needs_bridge:
            for a in template.GetAtoms():
                if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
                    target = charged.GetAtomWithIdx(match[a.GetIdx()])
                    target.SetFormalCharge(-1); target.SetNumExplicitHs(0); target.SetNoImplicit(True)
            Chem.SanitizeMol(charged)
        oxygens = [a.GetIdx() for a in cmp.GetAtoms() if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
        if len(oxygens) != 2:
            raise ValueError('Require source CMP dianion')
        offset = charged.GetNumAtoms()
        edit = Chem.RWMol(Chem.CombineMols(charged, cmp))
        edit.RemoveBond(match[pidx], match[oidx])
        edit.GetAtomWithIdx(offset + oxygens[0]).SetFormalCharge(0)
        edit.AddBond(match[pidx], offset + oxygens[0], Chem.BondType.SINGLE)
        Chem.SanitizeMol(edit)
        pair = Chem.GetMolFrags(edit, asMols=True)
        if len(pair) != 2:
            raise ValueError('Expected CDP-DAG and PG fragments')
        if any(sum(bool(exact_scaffold_matches(p, q)) for p in pair) != 1 for q in inputs):
            continue
        key = tuple(sorted(canonical(p) for p in pair))
        candidates[key] = {'reactant_smiles': list(key),
            'product_smiles': sorted([canonical(charged), canonical(cmp), canonical(proton)]),
            'charged_product_smiles': canonical(charged),
            'neutral_product_smiles': canonical(mol) if needs_bridge else None}
    return {'status': 'source-compatible-exact-precursor-hypothesis' if candidates else
            'source-precursor-stereo-mismatch', 'candidates': list(candidates.values())}


def build(audit, parent):
    source = audit['source_record']
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    original = set(compounds)
    used, reactions, rows = set(), {}, []
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    def equation(left, right, kind):
        sides = [[{'compound_id': c, 'coefficient': n} for c, n in
                  sorted(Counter(participant(s) for s in side).items())] for side in (left, right)]
        if not balanced(sides, compounds):
            raise ValueError('Cardiolipin equation failed element/isotope/charge balance')
        rid = stable_id('cardiolipin-synthesis-hypothesis', [kind, sides])
        reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
            'hypothesis_type': kind, 'source_reaction_id': source['rule_id'],
            'source_url': source['source_url'], 'enzyme_evidence_ids': [],
            'direction_status': 'forward-only-generic-source-hypothesis' if kind == 'cardiolipin-synthesis'
                                else 'acid-base-equilibrium-hypothesis',
            'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
        return rid
    for t in audit['targets']:
        result = instantiate(Chem.MolFromSmiles(t['canonical_smiles']), source)
        ids, required = [], set()
        for c in result['candidates']:
            ids.append(equation(c['reactant_smiles'], c['product_smiles'], 'cardiolipin-synthesis'))
            required.update(participant(s) for s in c['reactant_smiles'])
            if c['neutral_product_smiles']:
                ids.append(equation([c['charged_product_smiles'], '[H+]', '[H+]'],
                                    [c['neutral_product_smiles']], 'cardiolipin-protonation'))
        rows.append({**t, 'status': result['status'], 'reaction_ids': sorted(set(ids)),
                     'required_precursor_ids': sorted(required),
                     'precursors_absent_from_parent': sorted(required - original)})
    return {'schema': 'cannabis-carbon.phase1-cardiolipin-synthesis.v1', 'source_record': source,
        'claim_boundary': BOUNDARY, 'targets': rows, 'reactions': list(reactions.values()),
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'reviewed_records': len(rows), 'status_counts': dict(Counter(t['status'] for t in rows)),
            'balanced_equations': len(reactions),
            'equations_by_type': dict(Counter(r['hypothesis_type'] for r in reactions.values())),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-cardiolipin-identity-audit.json'),
             Path('data/reports/phase1-phosphatidate-hydrolysis-net.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-cardiolipin-synthesis.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
