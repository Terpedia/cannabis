import hashlib
import json
import unittest
from pathlib import Path

from rdkit import Chem
from cannabis_carbon.phase1_glycerolipid_precursors import instantiate, RULES
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced


def forward_products(reactants, kind):
    molecules = [Chem.MolFromSmiles(s) for s in reactants]
    donor = next(m for m in molecules if any(a.GetAtomicNum() == 16 for a in m.GetAtoms()))
    acceptor = next(m for m in molecules if m is not donor)
    hydrogen_count = 1 if kind == 'sn2-acylation' else 2
    oxygen = [a.GetIdx() for a in acceptor.GetAtoms() if a.GetAtomicNum() == 8
              and a.GetTotalNumHs() == 1 and a.GetDegree() == 1
              and a.GetNeighbors()[0].GetAtomicNum() == 6
              and a.GetNeighbors()[0].GetTotalNumHs() == hydrogen_count]
    sulfur = [a for a in donor.GetAtoms() if a.GetAtomicNum() == 16]
    assert len(oxygen) == len(sulfur) == 1
    carbon = [a.GetIdx() for a in sulfur[0].GetNeighbors()
              if any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetAtomicNum() == 8
                     for b in a.GetBonds())]
    assert len(carbon) == 1
    offset = acceptor.GetNumAtoms()
    edit = Chem.RWMol(Chem.CombineMols(acceptor, donor))
    edit.RemoveBond(offset + sulfur[0].GetIdx(), offset + carbon[0])
    edit.AddBond(oxygen[0], offset + carbon[0], Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted(canonical(m) for m in Chem.GetMolFrags(edit, asMols=True))


class GlycerolipidPrecursorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text())
        cls.sources = {k: next(r for r in catalog if r['rule_id'] == rid) for k, rid in RULES.items()}

    def test_both_regiospecific_steps_preserve_exact_product(self):
        examples = {'sn2-acylation': '[13CH3]CCC/C=C/CC(=O)OC[C@H](COP(=O)([O-])[O-])OC(=O)CCCC',
                    'sn1-acylation': '[13CH3]CCC/C=C/CC(=O)OC[C@@H](O)COP(=O)([O-])[O-]'}
        for kind, smiles in examples.items():
            mol = Chem.MolFromSmiles(smiles)
            rows = instantiate(mol, self.sources[kind])
            self.assertEqual(len(rows), 1)
            self.assertEqual(forward_products(rows[0]['reactant_smiles'], kind), rows[0]['product_smiles'])
            self.assertIn(canonical(mol), rows[0]['product_smiles'])
            unknown = Chem.Mol(mol); Chem.RemoveStereochemistry(unknown)
            self.assertEqual(instantiate(unknown, self.sources[kind]), [])
            wrong = Chem.Mol(mol)
            for a in wrong.GetAtoms():
                if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    a.InvertChirality()
            self.assertEqual(instantiate(wrong, self.sources[kind]), [])

    def test_complete_report_balance_roundtrip_and_missing_precursor_coverage(self):
        report = json.loads(Path('data/reports/phase1-glycerolipid-precursors.json').read_text())
        compounds = {c['id']: c for c in report['compounds']}
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            smiles = [[compounds[p['compound_id']]['smiles'] for p in r[side]
                       for _ in range(p['coefficient'])] for side in ('left', 'right')]
            self.assertEqual(forward_products(smiles[0], r['hypothesis_type']), sorted(smiles[1]))
            self.assertEqual(r['source_reaction_id'], RULES[r['hypothesis_type']])
        hydrolysis = json.loads(Path('data/reports/phase1-phosphatidate-hydrolysis.json').read_text())
        required = {c for r in hydrolysis['dag_candidates'] for c in r['precursors_absent_from_parent']}
        proposed = {r['compound_id'] for r in report['precursor_candidates'] if r['hypothesis_type'] == 'sn2-acylation'}
        self.assertEqual(len(required), 166)
        self.assertTrue(required <= proposed)
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)


if __name__ == '__main__':
    unittest.main()
