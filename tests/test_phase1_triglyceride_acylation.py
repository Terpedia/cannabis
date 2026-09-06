import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_triglyceride_acylation import instantiate, source_parts
from cannabis_carbon.phase1_marts_completions import balanced


class TriglycerideAcylationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = next(r for r in json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text()) if r['rule_id'] == 'RHEA:10869')

    def test_sn3_not_sn2_and_enantiomer_changes_transferred_chain(self):
        a = instantiate(Chem.MolFromSmiles('CCCC(=O)OC[C@H](COC(=O)CCCCCC)OC(=O)CCCCC'), self.source)
        b = instantiate(Chem.MolFromSmiles('CCCC(=O)OC[C@@H](COC(=O)CCCCCC)OC(=O)CCCCC'), self.source)
        self.assertEqual(len(a), 1); self.assertEqual(len(b), 1)
        donor_a = next(s for s in a[0]['reactant_smiles'] if 'S' in s)
        donor_b = next(s for s in b[0]['reactant_smiles'] if 'S' in s)
        self.assertTrue(donor_a.startswith('CCCCCCC(=O)S'))
        self.assertTrue(donor_b.startswith('CCCC(=O)S'))
        # The central six-carbon acyl group remains on glycerol in both cases.
        self.assertFalse(donor_a.startswith('CCCCCC(=O)S'))
        self.assertFalse(donor_b.startswith('CCCCCC(=O)S'))

    def test_unassigned_and_symmetric_targets_not_silently_stereospecified(self):
        for s in ('CCCC(=O)OCC(COC(=O)CCCCCC)OC(=O)CCCCC',
                  'CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCCCC'):
            self.assertEqual(instantiate(Chem.MolFromSmiles(s), self.source), [])

    def test_isotope_and_double_bond_geometry_retained(self):
        rows = instantiate(Chem.MolFromSmiles('CCCC(=O)OC[C@H](COC(=O)[13CH2]/C=C/CCC)OC(=O)CCCCC'), self.source)
        self.assertEqual(len(rows), 1)
        donor = next(s for s in rows[0]['reactant_smiles'] if 'S' in s)
        self.assertIn('[13CH2]', donor); self.assertIn('/', donor)

    def test_all_generated_equations_replay_and_sources_match(self):
        report = json.loads(Path('data/reports/phase1-triglyceride-acylation.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in report['compounds']}
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            self.assertEqual(r['direction_status'], 'forward-only-generic-source-hypothesis')
            self.assertEqual(r['enzyme_evidence_ids'], [])
        self.assertEqual(len(report['reactions']), 1574)
        self.assertEqual(report['summary']['previous_no_producer_records_with_hypotheses'], 1571)
        self.assertEqual(report['summary']['new_CO2_pathway_claims'], 0)
