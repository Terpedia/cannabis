import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_triglyceride_symmetry import analyze
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced


def forward_products(reactants):
    molecules = [Chem.MolFromSmiles(s) for s in reactants]
    donor = next(m for m in molecules if any(a.GetAtomicNum() == 16 for a in m.GetAtoms()))
    dag = next(m for m in molecules if m is not donor)
    oxygen = [a.GetIdx() for a in dag.GetAtoms() if a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 1]
    sulfur = [a.GetIdx() for a in donor.GetAtoms() if a.GetAtomicNum() == 16]
    assert len(oxygen) == len(sulfur) == 1
    carbons = [a.GetIdx() for a in donor.GetAtomWithIdx(sulfur[0]).GetNeighbors()
               if any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetAtomicNum() == 8 for b in a.GetBonds())]
    assert len(carbons) == 1
    offset = dag.GetNumAtoms()
    edit = Chem.RWMol(Chem.CombineMols(dag, donor))
    edit.RemoveBond(offset + sulfur[0], offset + carbons[0])
    edit.AddBond(oxygen[0], offset + carbons[0], Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    parts = Chem.GetMolFrags(edit, asMols=True)
    for m in parts:
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return sorted(canonical(m) for m in parts)


class SymmetricTriglycerideTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = json.loads(Path('data/reports/phase1-triglyceride-acylation.json').read_text())['source_record']

    def test_symmetric_product_has_stereospecific_precursor_and_exact_round_trip(self):
        target = Chem.MolFromSmiles('CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCCCC')
        result = analyze(target, self.source)
        self.assertEqual(len(result['candidates']), 1)
        candidate = result['candidates'][0]
        dag = next(s for s in candidate['reactant_smiles'] if 'S' not in s)
        self.assertIn('@', dag)
        self.assertNotIn('@', canonical(target))
        self.assertEqual(forward_products(candidate['reactant_smiles']), candidate['product_smiles'])

    def test_unknown_asymmetric_target_is_not_promoted(self):
        result = analyze(Chem.MolFromSmiles('CCCC(=O)OCC(COC(=O)CCCCCC)OC(=O)CCCCC'), self.source)
        self.assertEqual(result['candidates'], [])
        self.assertIn('separate-identity-review', result['status'])

    def test_all_source_targets_reconstruct_without_stereo_or_identity_loss(self):
        report = json.loads(Path('data/reports/phase1-triglyceride-symmetry.json').read_text())
        compounds = {c['id']: c for c in report['compounds']}
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            reactants = [compounds[p['compound_id']]['smiles'] for p in r['left'] for _ in range(p['coefficient'])]
            products = sorted(compounds[p['compound_id']]['smiles'] for p in r['right'] for _ in range(p['coefficient']))
            self.assertEqual(forward_products(reactants), products)
        self.assertEqual(report['summary']['balanced_equations'], 129)
        self.assertEqual(sum(not t['reaction_ids'] for t in report['targets']), 14)
