"""Exact PE/PS synthesis and PE methylation proposals with separate speciation."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches
from .phase1_marts_completions import balanced

KINDS = {'PE': 'pe-synthesis', 'PS': 'ps-synthesis',
         'PE-NMe': 'pe-first-methylation', 'PE-NMe2': 'pe-second-methylation'}
BOUNDARY = ('Generic-source forward synthesis hypothesis, not established exact-substrate '
    'activity in Cannabis. Complete source reactants and coproducts are explicit. '
    'Acyl-chain identity, isotopes and encoded stereochemistry are preserved. Charge-state '
    'conversion is a separate reversible sensitivity assumption, including explicitly '
    'identified net-zero internal proton relocation; identities are never merged. '
    'No organic carbon currency is supplied externally. Precursor synthesis, physiological '
    'direction, pH/pKa, compartments, energetics and enzyme activity remain unresolved '
    'unless separately evidenced. Balanced participation alone adds no CO2 route claim.')


def instantiate(mol, source, kind):
    left, right = [[Chem.MolFromSmiles(s) for s in side.split('.')]
                   for side in source['reaction_smarts'].split('>>')]
    template = next(m for m in right if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
    fixed = [m for m in right if all(a.GetAtomicNum() != 0 for a in m.GetAtoms())]
    other = next(m for m in fixed if m.GetNumAtoms() > 1)
    result = {}
    for match in exact_scaffold_matches(mol, template, uniquify=False):
        edit = Chem.RWMol(Chem.CombineMols(mol, other)); offset = mol.GetNumAtoms()
        if kind in ('pe-synthesis', 'ps-synthesis'):
            p, = [a for a in template.GetAtoms() if a.GetAtomicNum() == 15]
            choices = []
            for o in p.GetNeighbors():
                if o.GetAtomicNum() != 8 or o.GetDegree() != 2:
                    continue
                # Remove the glycerol-side ester for PE; the serine-side ester for PS.
                cut = Chem.RWMol(template); cut.RemoveBond(p.GetIdx(), o.GetIdx())
                fragment = next(f for f in Chem.GetMolFrags(cut) if o.GetIdx() in f)
                has_tail = any(template.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in fragment)
                has_n = any(template.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in fragment)
                if (kind == 'pe-synthesis' and has_tail and not has_n) or (kind == 'ps-synthesis' and has_n and not has_tail):
                    choices.append(o.GetIdx())
            oxygen, = choices
            cmp_oxygens = [a.GetIdx() for a in other.GetAtoms() if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
            if len(cmp_oxygens) != 2:
                raise ValueError('Expected exact CMP dianion')
            edit.RemoveBond(match[p.GetIdx()], match[oxygen])
            edit.GetAtomWithIdx(offset + cmp_oxygens[0]).SetFormalCharge(0)
            edit.AddBond(match[p.GetIdx()], offset + cmp_oxygens[0], Chem.BondType.SINGLE)
        elif kind in ('pe-first-methylation', 'pe-second-methylation'):
            n, = [a for a in template.GetAtoms() if a.GetAtomicNum() == 7]
            methyls = [a.GetIdx() for a in n.GetNeighbors() if a.GetAtomicNum() == 6 and a.GetDegree() == 1]
            if len(methyls) != (1 if kind == 'pe-first-methylation' else 2):
                raise ValueError('Expected exact N-methyl substitution')
            sulfur, = [a.GetIdx() for a in other.GetAtoms() if a.GetAtomicNum() == 16 and a.GetFormalCharge() == 0 and a.GetDegree() == 2]
            edit.RemoveBond(match[n.GetIdx()], match[methyls[0]])
            atom = edit.GetAtomWithIdx(match[n.GetIdx()])
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            edit.GetAtomWithIdx(offset + sulfur).SetFormalCharge(1)
            edit.AddBond(match[methyls[0]], offset + sulfur, Chem.BondType.SINGLE)
        else:
            raise ValueError('Unknown synthesis class')
        Chem.SanitizeMol(edit)
        parts = list(Chem.GetMolFrags(edit, asMols=True))
        if len(parts) != len(left) or any(sum(bool(exact_scaffold_matches(p, q)) for p in parts) != 1 for q in left):
            raise ValueError('Generated precursors fail exact source scaffold/stereo checks')
        reactants = sorted(canonical(m) for m in parts)
        products = sorted([canonical(mol)] + [canonical(m) for m in fixed])
        result[tuple(reactants)] = {'reactant_smiles': reactants, 'product_smiles': products}
    return [result[k] for k in sorted(result)]


def build(parent, audit):
    compounds = {c['id']: dict(c) for c in parent['compounds']}
    known = set(compounds)
    for c in audit['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Speciation identity conflict')
        compounds.setdefault(c['id'], dict(c))
    sources = {s['rule_id']: s for s in audit['source_records']}
    reactions, used, targets = {}, set(), []
    def participant(smiles):
        mol = Chem.MolFromSmiles(smiles); smiles = canonical(mol); cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles, 'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol), 'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    for t in audit['targets']:
        kind = KINDS[t['lipid_class']]
        source = sources[t['source_reaction_id']]
        mol = Chem.MolFromSmiles(compounds[t['source_form_compound_id']]['smiles'])
        synthesis_ids, required = [], set()
        for candidate in instantiate(mol, source, kind):
            sides = [[{'compound_id': cid, 'coefficient': n} for cid, n in
                      sorted(Counter(participant(s) for s in candidate[side]).items())]
                     for side in ('reactant_smiles', 'product_smiles')]
            if not balanced(sides, compounds):
                raise ValueError('Synthesis is not element/isotope/charge balanced')
            rid = stable_id('amino-phospholipid-synthesis', [kind, sides])
            reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                'hypothesis_type': kind, 'source_reaction_id': source['rule_id'], 'source_url': source['source_url'],
                'direction_status': 'forward-only-generic-source-hypothesis', 'enzyme_evidence_ids': [],
                'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY})
            synthesis_ids.append(rid); required.update(m['compound_id'] for m in sides[0])
        if not synthesis_ids:
            raise ValueError('Audited source form lacks synthesis instantiation')
        if not balanced([t['left'], t['right']], compounds):
            raise ValueError('Speciation failed independent balance')
        rid = stable_id('amino-phospholipid-speciation', [t['left'], t['right']])
        reactions.setdefault(rid, {'id': rid, 'left': t['left'], 'right': t['right'],
            'hypothesis_type': 'amino-phospholipid-speciation', 'speciation_type': t['speciation_type'],
            'source_reaction_id': source['rule_id'], 'source_url': source['source_url'],
            'direction_status': 'explicit-reversible-speciation-assumption-not-source-reaction',
            'enzyme_evidence_ids': [], 'balance_status': 'independently-element-isotope-charge-balanced',
            'claim_boundary': t['claim_boundary'] + ' Reversibility is an explicit sensitivity assumption.'})
        used.update(m['compound_id'] for side in (t['left'], t['right']) for m in side)
        targets.append({k:t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'lipid_class')} | {
            'synthesis_reaction_ids': synthesis_ids, 'speciation_reaction_id': rid,
            'required_precursor_ids': sorted(required), 'precursors_absent_from_parent': sorted(required - known)})
    return {'schema': 'cannabis-carbon.phase1-amino-phospholipid-synthesis.v1', 'source_records': list(sources.values()),
        'reactions': list(reactions.values()), 'compounds': [compounds[c] for c in sorted(used)],
        'targets': targets, 'claim_boundary': BOUNDARY,
        'summary': {'target_records': len(targets), 'balanced_equations': len(reactions),
            'equations_by_type': dict(Counter(r['hypothesis_type'] for r in reactions.values())),
            'targets_with_precursors_absent_from_parent': sum(bool(t['precursors_absent_from_parent']) for t in targets),
            'distinct_missing_precursors': len({p for t in targets for p in t['precursors_absent_from_parent']}),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-source-mapped-protonation-net.json'), Path('data/reports/phase1-amino-phospholipid-speciation.json')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-amino-phospholipid-synthesis.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
