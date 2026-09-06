import copy
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_lipid_acylation import (
    build, canonical, exact_scaffold_matches, precursors, sn2_bond)
from cannabis_carbon.phase1_marts_completions import balanced


class LipidAcylationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text())
        cls.network = json.loads(Path('data/reports/phase1-full-balanced-network.json').read_text())

    def small_network(self):
        return {**self.network, 'targets': [t for t in self.network['targets']
            if t['cannabisdb_id'] in ('CDB000803', 'CDB000806', 'CDB000988')]}

    def test_exact_balance_source_scaffold_and_identity_preservation(self):
        network = self.small_network()
        before = copy.deepcopy(network)
        report = build(network, self.raw)
        self.assertEqual(before, network)
        self.assertEqual(len(report['targets']), 3)
        self.assertEqual(report['summary']['reaction_type_counts'], {'explicit-PA-protonation': 2, 'sn2-acylation': 3})
        compounds = {c['id']: c for c in report['compounds']}
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
        unspecified = next(t for t in report['targets'] if t['cannabisdb_id'] == 'CDB000988')
        self.assertEqual(unspecified['encoded_structure_status'], 'stereo-unspecified-or-unknown')
        self.assertNotIn('/', compounds[unspecified['reaction_product_id']]['smiles'])
        self.assertNotIn('\\', compounds[unspecified['reaction_product_id']]['smiles'])

    def test_wrong_sn_configuration_and_headgroup_rejected(self):
        source = next(r for r in self.raw if r['rule_id'] == 'RHEA:12938')
        template = Chem.MolFromSmiles(source['reaction_smarts'].split('>>')[1].split('.')[0])
        product = Chem.MolFromSmiles('CCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCC')
        self.assertTrue(exact_scaffold_matches(product, template))
        wrong = Chem.MolFromSmiles('CCCC(=O)OC[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCC')
        self.assertFalse(exact_scaffold_matches(wrong, template))
        extra = Chem.MolFromSmiles('CCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)CC)OC(=O)CCC')
        self.assertFalse(exact_scaffold_matches(extra, template))

    def test_isotopic_and_double_bond_information_survives_transfer(self):
        source = next(r for r in self.raw if r['rule_id'] == 'RHEA:12938')
        template, coa = [Chem.MolFromSmiles(s) for s in source['reaction_smarts'].split('>>')[1].split('.')]
        product = Chem.MolFromSmiles('CCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)[13CH2]/C=C/CC')
        match = exact_scaffold_matches(product, template)[0]
        ci, oi = sn2_bond(template)
        parts = precursors(product, coa, match[ci], match[oi])
        donor = next(m for m in parts if any(a.GetAtomicNum() == 16 for a in m.GetAtoms()))
        self.assertIn('[13CH2]', canonical(donor))
        self.assertIn('/', canonical(donor))
        lyso = next(m for m in parts if m is not donor)
        self.assertIn('CCCC(=O)O', canonical(lyso))

    def test_cyclic_acyl_and_unassigned_sn_excluded(self):
        source = next(r for r in self.raw if r['rule_id'] == 'RHEA:12938')
        template = Chem.MolFromSmiles(source['reaction_smarts'].split('>>')[1].split('.')[0])
        for smiles in ('CCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCC',
                       'CCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)C1CCCCC1'):
            self.assertFalse(exact_scaffold_matches(Chem.MolFromSmiles(smiles), template))
