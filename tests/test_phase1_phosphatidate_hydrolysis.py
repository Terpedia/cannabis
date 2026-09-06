import hashlib
import json
import unittest
from pathlib import Path

from rdkit import Chem
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_phosphatidate_hydrolysis import instantiate


def hydrolyze(reactants):
    molecules = [Chem.MolFromSmiles(s) for s in reactants]
    pa = next(m for m in molecules if m.GetNumAtoms() > 1)
    water = next(m for m in molecules if m.GetNumAtoms() == 1)
    p = [a for a in pa.GetAtoms() if a.GetAtomicNum() == 15]
    assert len(p) == 1 and canonical(water) == 'O'
    ester = [a for a in p[0].GetNeighbors() if a.GetAtomicNum() == 8
             and any(n.GetAtomicNum() == 6 for n in a.GetNeighbors())]
    assert len(ester) == 1
    edit = Chem.RWMol(Chem.CombineMols(pa, water))
    edit.RemoveBond(p[0].GetIdx(), ester[0].GetIdx())
    edit.AddBond(p[0].GetIdx(), pa.GetNumAtoms(), Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted(canonical(m) for m in Chem.GetMolFrags(edit, asMols=True))


class PhosphatidateHydrolysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = next(r for r in json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text())
                          if r['rule_id'] == 'RHEA:27430')

    def test_exact_forward_replay_preserves_tail_geometry_and_isotope(self):
        mol = Chem.MolFromSmiles('[13CH3]CCC/C=C/CC(=O)OC[C@H](CO)OC(=O)CCCC')
        candidates = instantiate(mol, self.source)
        self.assertEqual(len(candidates), 1)
        row = candidates[0]
        self.assertEqual(hydrolyze(row['reactant_smiles']), row['product_smiles'])
        self.assertIn(canonical(mol), row['product_smiles'])
        pa = next(s for s in row['reactant_smiles'] if s != 'O')
        self.assertEqual(Chem.GetFormalCharge(Chem.MolFromSmiles(pa)), -2)
        self.assertIn('[13CH3]', pa)
        self.assertIn('/', pa)

    def test_unknown_wrong_sn_and_other_headgroups_not_promoted(self):
        for smiles in ('CCCC(=O)OCC(CO)OC(=O)CCCC',
                       'CCCC(=O)OC[C@@H](CO)OC(=O)CCCC',
                       'CCCC(=O)OC[C@H](COP(=O)(O)O)OC(=O)CCCC',
                       'CCCC(=O)OC[C@H](COC(=O)CCCC)O'):
            with self.subTest(smiles=smiles):
                self.assertEqual(instantiate(Chem.MolFromSmiles(smiles), self.source), [])

    def test_all_reported_equations_balance_and_replay_exactly(self):
        report = json.loads(Path('data/reports/phase1-phosphatidate-hydrolysis.json').read_text())
        compounds = {c['id']: c for c in report['compounds']}
        self.assertGreater(len(report['reactions']), 0)
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            smiles = [[compounds[p['compound_id']]['smiles'] for p in r[side]
                       for _ in range(p['coefficient'])] for side in ('left', 'right')]
            self.assertEqual(hydrolyze(smiles[0]), sorted(smiles[1]))
            self.assertEqual(r['source_reaction_id'], 'RHEA:27430')
            self.assertEqual(r['enzyme_evidence_ids'], [])
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)


if __name__ == '__main__':
    unittest.main()
