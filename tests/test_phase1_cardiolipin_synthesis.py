import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_cardiolipin_synthesis import instantiate
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced


def forward_products(reactants):
    molecules = [Chem.MolFromSmiles(s) for s in reactants]
    donor = next(m for m in molecules if any(a.GetAtomicNum() == 7 for a in m.GetAtoms()))
    pg = next(m for m in molecules if m is not donor)
    bridges = [a for a in donor.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 2
               and all(n.GetAtomicNum() == 15 for n in a.GetNeighbors())]
    assert len(bridges) == 1
    bridge = bridges[0]
    choices = []
    for p in bridge.GetNeighbors():
        cut = Chem.RWMol(donor); cut.RemoveBond(p.GetIdx(), bridge.GetIdx())
        f = next(f for f in Chem.GetMolFrags(cut) if p.GetIdx() in f)
        if not any(donor.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in f):
            choices.append(p.GetIdx())
    hydroxyls = [a.GetIdx() for a in pg.GetAtoms() if a.GetAtomicNum() == 8
                 and a.GetDegree() == 1 and a.GetTotalNumHs() == 1
                 and a.GetNeighbors()[0].GetAtomicNum() == 6
                 and a.GetNeighbors()[0].GetTotalNumHs() == 2]
    assert len(choices) == len(hydroxyls) == 1
    edit = Chem.RWMol(Chem.CombineMols(donor, pg))
    edit.RemoveBond(choices[0], bridge.GetIdx())
    edit.GetAtomWithIdx(bridge.GetIdx()).SetFormalCharge(-1)
    edit.AddBond(choices[0], donor.GetNumAtoms() + hydroxyls[0], Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted([canonical(m) for m in Chem.GetMolFrags(edit, asMols=True)] + ['[H+]'])


class CardiolipinSynthesisTests(unittest.TestCase):
    def test_all_synthesis_steps_replay_exactly_and_protonation_is_separate(self):
        report = json.loads(Path('data/reports/phase1-cardiolipin-synthesis.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in report['compounds']}
        self.assertEqual(report['summary']['equations_by_type'],
                         {'cardiolipin-synthesis': 92, 'cardiolipin-protonation': 92})
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            sides = [[compounds[p['compound_id']]['smiles'] for p in r[side]
                      for _ in range(p['coefficient'])] for side in ('left', 'right')]
            if r['hypothesis_type'] == 'cardiolipin-synthesis':
                self.assertEqual(forward_products(sides[0]), sorted(sides[1]))
            else:
                self.assertEqual(sides[0].count('[H+]'), 2)
                mol = Chem.MolFromSmiles(next(s for s in sides[0] if s != '[H+]'))
                for a in mol.GetAtoms():
                    if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
                        a.SetFormalCharge(0); a.SetNoImplicit(False)
                Chem.SanitizeMol(mol)
                self.assertEqual([canonical(mol)], sides[1])
        unknown = [t for t in report['targets'] if t['central_glycerol_stereo'] == 'Unspecified']
        self.assertEqual(len(unknown), 1273)
        self.assertTrue(all(not t['reaction_ids'] for t in unknown))
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)

    def test_atom_order_does_not_select_a_different_precursor_orientation(self):
        report = json.loads(Path('data/reports/phase1-cardiolipin-synthesis.json').read_text())
        target = next(t for t in report['targets'] if t['reaction_ids'])
        mol = Chem.MolFromSmiles(target['canonical_smiles'])
        reversed_mol = Chem.RenumberAtoms(mol, list(range(mol.GetNumAtoms()))[::-1])
        self.assertEqual(instantiate(mol, report['source_record']),
                         instantiate(reversed_mol, report['source_record']))
