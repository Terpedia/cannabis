"""Exact scaffold audit of amino-phospholipid charge-state gaps, not new routes."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical, exact_scaffold_matches
from .phase1_marts_completions import balanced

RULES = {'PE': 'RHEA:32944', 'PE-NMe': 'RHEA:11165',
         'PE-NMe2': 'RHEA:32736', 'PS': 'RHEA:16914'}
BOUNDARY = ('Review-only exact scaffold speciation proposal, not an exact ChEBI mapping '
    'or curated reaction. Neutral and zwitterionic identities remain distinct. Only '
    'local proton/charge changes on source headgroup N/O atoms are considered; heavy '
    'atoms, bonds, isotopes and encoded stereochemistry remain unchanged. Net-zero '
    'proton relocation is explicitly distinguished from exchange with H+. Source '
    'scaffolds do not establish exact lipid specificity, tissue pH, pKa, transport, '
    'enzyme activity or a CO2 route. No model coverage is changed by this audit.')


def neutral_template(template):
    edit = Chem.RWMol(template)
    for a in edit.GetAtoms():
        q = a.GetFormalCharge()
        if not q:
            continue
        if not ((a.GetAtomicNum() == 7 and q == 1 and a.GetTotalNumHs() >= 1)
                or (a.GetAtomicNum() == 8 and q == -1)):
            raise ValueError('Unexpected charged scaffold atom')
        a.SetNumExplicitHs(a.GetTotalNumHs() - q)
        a.SetNoImplicit(True)
        a.SetFormalCharge(0)
    Chem.SanitizeMol(edit)
    return edit.GetMol()


def source_forms(mol, template):
    """Instantiate only source headgroup charges after an exact neutral match."""
    result = {}
    for match in exact_scaffold_matches(mol, neutral_template(template), uniquify=False):
        edit = Chem.RWMol(mol)
        for atom in template.GetAtoms():
            if atom.GetFormalCharge():
                a = edit.GetAtomWithIdx(match[atom.GetIdx()])
                a.SetFormalCharge(atom.GetFormalCharge())
                a.SetNumExplicitHs(atom.GetTotalNumHs())
                a.SetNoImplicit(True)
        Chem.SanitizeMol(edit)
        generated = edit.GetMol()
        if not exact_scaffold_matches(generated, template):
            raise ValueError('Generated source form fails exact scaffold')
        # Replay the local change, including the net-zero proton-relocation case.
        for a, b in zip(mol.GetAtoms(), generated.GetAtoms()):
            if (a.GetAtomicNum(), a.GetIsotope(), a.GetChiralTag()) != (b.GetAtomicNum(), b.GetIsotope(), b.GetChiralTag()):
                raise ValueError('Heavy atom identity changed')
            if b.GetTotalNumHs() - a.GetTotalNumHs() != b.GetFormalCharge() - a.GetFormalCharge():
                raise ValueError('Not a proton-only change')
        bonds = lambda m: [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType(), b.GetStereo(), b.GetBondDir()) for b in m.GetBonds()]
        if bonds(mol) != bonds(generated):
            raise ValueError('Bond or stereochemistry changed')
        result[canonical(generated)] = generated
    return list(result.values())


def build(parent, catalog, *, rules=None, boundary=BOUNDARY,
          schema='cannabis-carbon.phase1-amino-phospholipid-speciation.v1'):
    compounds = {c['id']: c for c in parent['compounds']}
    sources = {kind: next(r for r in catalog if r['rule_id'] == rid)
               for kind, rid in (RULES if rules is None else rules).items()}
    templates = {kind: next(m for m in (Chem.MolFromSmiles(s) for s in r['reaction_smarts'].split('>>')[1].split('.'))
                           if any(a.GetAtomicNum() == 0 for a in m.GetAtoms())) for kind, r in sources.items()}
    known = set(compounds)
    rows, used = [], set()
    def participant(mol):
        smiles = canonical(mol); cid = stable_id('structure', smiles)
        compounds.setdefault(cid, {'id': cid, 'smiles': smiles, 'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol), 'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
        used.add(cid)
        return cid
    for t in parent['targets']:
        mol = Chem.MolFromSmiles(compounds[t['compound_id']]['smiles'])
        for kind, template in templates.items():
            for source_form in source_forms(mol, template):
                cid = participant(source_form)
                used.add(t['compound_id'])
                delta = Chem.GetFormalCharge(mol) - Chem.GetFormalCharge(source_form)
                left = [{'compound_id': cid, 'coefficient': 1}]
                right = [{'compound_id': t['compound_id'], 'coefficient': 1}]
                if delta:
                    (left if delta > 0 else right).append({'compound_id': participant(Chem.MolFromSmiles('[H+]')), 'coefficient': abs(delta)})
                if not balanced([left, right], compounds):
                    raise ValueError('Speciation proposal is not balanced')
                rows.append({**{k:t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'net_status')},
                    'lipid_class': kind, 'source_form_compound_id': cid,
                    'source_form_in_parent_inventory': cid in known,
                    'speciation_type': 'net-zero-intramolecular-proton-relocation' if delta == 0 else 'explicit-proton-exchange',
                    'protons_consumed': delta, 'left': left, 'right': right,
                    'source_reaction_id': sources[kind]['rule_id'], 'source_url': sources[kind]['source_url'],
                    'claim_boundary': boundary})
    return {'schema': schema,
        'source_records': list(sources.values()), 'targets': rows,
        'compounds': [compounds[c] for c in sorted(used)], 'claim_boundary': boundary,
        'summary': {'inventory_records_screened': len(parent['targets']), 'matched_records': len({t['cannabisdb_id'] for t in rows}),
            'proposal_count': len(rows), 'lipid_class_counts': dict(Counter(t['lipid_class'] for t in rows)),
            'speciation_type_counts': dict(Counter(t['speciation_type'] for t in rows)),
            'parent_status_counts': dict(Counter(t['net_status'] for t in rows)),
            'source_forms_in_parent_inventory': sum(t['source_form_in_parent_inventory'] for t in rows),
            'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-source-mapped-protonation-net.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_bytes()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-amino-phospholipid-speciation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
