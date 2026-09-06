"""Source-exact PG, PGP and CDP-DAG precursor reaction proposals."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches
from .phase1_marts_completions import balanced

RULES = {'pgp-hydrolysis': 'RHEA:33752', 'pgp-synthesis': 'RHEA:12594',
         'cdp-dag-synthesis': 'RHEA:16230'}
BOUNDARY = ('Generic source-forward reaction instantiation only, not demonstrated '
    'Cannabis activity. Exact PG/PGP/CDP-DAG/PA identities, acyl chains, isotopes and '
    'encoded stereochemistry are preserved. CTP, CMP, glycerol phosphate, water, '
    'phosphate, diphosphate and H+ are explicit where present in the source equation. '
    'No carbon-containing currency is an external seed. Balanced participation '
    'does not establish precursor supply, CO2 conversion, pool startup or energetics.')


def instantiate(mol, source):
    kind = next(k for k, rid in RULES.items() if rid == source['rule_id'])
    left, right = [[Chem.MolFromSmiles(s) for s in side.split('.')]
                   for side in source['reaction_smarts'].split('>>')]
    template = next(m for m in right if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
    fixed = [m for m in right if all(a.GetAtomicNum() != 0 for a in m.GetAtoms())]
    matches = exact_scaffold_matches(mol, template, uniquify=False)
    if not matches:
        return []
    other = next(m for m in fixed if m.GetNumAtoms() > 1)
    candidates = {}
    for match in matches:
        edit = Chem.RWMol(Chem.CombineMols(mol, other))
        offset = mol.GetNumAtoms()
        if kind == 'pgp-hydrolysis':
            hydroxyl = [a for a in template.GetAtoms() if a.GetAtomicNum() == 8
                        and a.GetTotalNumHs() == 1 and a.GetDegree() == 1
                        and a.GetNeighbors()[0].GetAtomicNum() == 6
                        and a.GetNeighbors()[0].GetTotalNumHs() == 2]
            oh = [a for a in other.GetAtoms() if a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 1]
            if len(hydroxyl) != 1 or len(oh) != 1:
                raise ValueError('Require unique terminal PG OH and phosphate OH')
            p = oh[0].GetNeighbors()[0].GetIdx()
            edit.RemoveBond(offset + p, offset + oh[0].GetIdx())
            edit.AddBond(match[hydroxyl[0].GetIdx()], offset + p, Chem.BondType.SINGLE)
        elif kind == 'pgp-synthesis':
            centers = [a for a in template.GetAtoms() if a.GetAtomicNum() == 6
                       and any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1 for n in a.GetNeighbors())]
            choices = []
            for center in centers:
                for carbon in center.GetNeighbors():
                    if carbon.GetAtomicNum() != 6:
                        continue
                    for oxygen in carbon.GetNeighbors():
                        if oxygen.GetAtomicNum() != 8:
                            continue
                        for p in oxygen.GetNeighbors():
                            if p.GetAtomicNum() == 15 and sum(any(n.GetAtomicNum() == 6 for n in o.GetNeighbors())
                                for o in p.GetNeighbors() if o.GetAtomicNum() == 8) == 2:
                                choices.append((p.GetIdx(), oxygen.GetIdx()))
            oxygens = [a.GetIdx() for a in other.GetAtoms() if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
            if len(choices) != 1 or len(oxygens) != 2:
                raise ValueError('Require unique phosphatidyl linkage and source CMP')
            p, o = choices[0]
            edit.RemoveBond(match[p], match[o])
            edit.GetAtomWithIdx(offset + oxygens[0]).SetFormalCharge(0)
            edit.AddBond(match[p], offset + oxygens[0], Chem.BondType.SINGLE)
        else:
            bridges = [a for a in template.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 2
                       and all(n.GetAtomicNum() == 15 for n in a.GetNeighbors())]
            if len(bridges) != 1:
                raise ValueError('Require CDP-DAG pyrophosphate bridge')
            bridge = bridges[0]; choices = []
            for p in bridge.GetNeighbors():
                cut = Chem.RWMol(template); cut.RemoveBond(p.GetIdx(), bridge.GetIdx())
                frag = next(f for f in Chem.GetMolFrags(cut) if p.GetIdx() in f)
                if any(template.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in frag):
                    choices.append(p.GetIdx())
            oh = [a.GetIdx() for a in other.GetAtoms() if a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 1]
            if len(choices) != 1 or len(oh) != 1:
                raise ValueError('Require cytidine-side phosphate and source diphosphate OH')
            edit.RemoveBond(match[choices[0]], match[bridge.GetIdx()])
            edit.GetAtomWithIdx(match[bridge.GetIdx()]).SetFormalCharge(-1)
            edit.AddBond(match[choices[0]], offset + oh[0], Chem.BondType.SINGLE)
        Chem.SanitizeMol(edit)
        parts = list(Chem.GetMolFrags(edit, asMols=True))
        if kind == 'cdp-dag-synthesis':
            parts.append(Chem.MolFromSmiles('[H+]'))
        if len(parts) != len(left) or any(sum(bool(exact_scaffold_matches(p, q)) for p in parts) != 1 for q in left):
            raise ValueError('Generated inputs do not match source structures/stereochemistry')
        key = tuple(sorted(canonical(p) for p in parts))
        candidates[key] = {'reactant_smiles': list(key),
            'product_smiles': sorted([canonical(mol)] + [canonical(m) for m in fixed])}
    return list(candidates.values())


def build(parent, cardiolipin, catalog):
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    for c in cardiolipin['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Compound identity conflict')
        compounds.setdefault(c['id'], dict(c))
    original = set(compounds)
    sources = {k: next(r for r in catalog if r['rule_id'] == rid) for k, rid in RULES.items()}
    used, reactions, rows = set(), {}, []
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol)
        cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles,
            'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
            'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    # Ordered backward expansion: every PG, then every required PGP, then CDP-DAG.
    for kind, source in sources.items():
        for cid in sorted(compounds):
            for c in instantiate(Chem.MolFromSmiles(compounds[cid]['smiles']), source):
                sides = [[{'compound_id': x, 'coefficient': n} for x, n in
                          sorted(Counter(participant(s) for s in c[side]).items())]
                         for side in ('reactant_smiles', 'product_smiles')]
                if not balanced(sides, compounds):
                    raise ValueError('Precursor equation is not element/isotope/charge balanced')
                rid = stable_id('cardiolipin-precursor-hypothesis', [kind, sides])
                reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                    'hypothesis_type': kind, 'source_reaction_id': source['rule_id'], 'source_url': source['source_url'],
                    'direction_status': 'forward-only-generic-source-hypothesis', 'enzyme_evidence_ids': [],
                    'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
                required = {x['compound_id'] for x in sides[0]}
                rows.append({'compound_id': cid, 'reaction_id': rid, 'hypothesis_type': kind,
                    'required_precursor_ids': sorted(required), 'precursors_absent_from_input_inventory': sorted(required - original)})
    missing = {c for t in cardiolipin['targets'] for c in t['precursors_absent_from_parent']}
    proposed = {r['compound_id'] for r in rows}
    return {'schema': 'cannabis-carbon.phase1-cardiolipin-precursors.v1',
        'claim_boundary': BOUNDARY, 'source_records': list(sources.values()),
        'precursor_candidates': rows, 'reactions': list(reactions.values()),
        'compounds': [compounds[c] for c in sorted(used)],
        'summary': {'input_compound_structures': len(original), 'balanced_equations': len(reactions),
            'equations_by_type': dict(Counter(r['hypothesis_type'] for r in reactions.values())),
            'missing_cardiolipin_precursors_with_synthesis_hypotheses': len(missing & proposed),
            'missing_cardiolipin_precursors_without_synthesis_hypotheses': len(missing - proposed),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-phosphatidate-hydrolysis-net.json'),
             Path('data/reports/phase1-cardiolipin-synthesis.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-cardiolipin-precursors.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
