import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_amino_phospholipid_speciation import source_forms
from cannabis_carbon.phase1_lipid_acylation import canonical, exact_scaffold_matches
from cannabis_carbon.phase1_marts_completions import balanced


class AminoPhospholipidSpeciationTests(unittest.TestCase):
    def setUp(self):
        self.report = json.loads(Path('data/reports/phase1-amino-phospholipid-speciation.json').read_bytes())

    def test_all_proposals_replay_without_merging_identities_or_claiming_coverage(self):
        r = self.report
        parent = json.loads(Path('data/reports/phase1-source-mapped-protonation-net.json').read_bytes())
        original = {t['cannabisdb_id']: t for t in parent['targets']}
        cs = {c['id']: c for c in r['compounds']}
        templates = {s['rule_id']: next(m for m in map(Chem.MolFromSmiles, s['reaction_smarts'].split('>>')[1].split('.'))
            if any(a.GetAtomicNum() == 0 for a in m.GetAtoms())) for s in r['source_records']}
        for path, sha in r['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), sha)
        self.assertEqual(r['summary']['inventory_records_screened'], 6220)
        self.assertEqual(r['summary']['matched_records'], 162)
        self.assertEqual(r['summary']['new_CO2_route_claims'], 0)
        self.assertEqual(Counter(t['speciation_type'] for t in r['targets']),
                         {'explicit-proton-exchange': 78, 'net-zero-intramolecular-proton-relocation': 84})
        for t in r['targets']:
            self.assertEqual(t['compound_id'], original[t['cannabisdb_id']]['compound_id'])
            self.assertEqual(t['net_status'], original[t['cannabisdb_id']]['net_status'])
            self.assertNotEqual(t['compound_id'], t['source_form_compound_id'])
            self.assertTrue(balanced([t['left'], t['right']], cs))
            mol = Chem.MolFromSmiles(cs[t['compound_id']]['smiles'])
            template = templates[t['source_reaction_id']]
            forms = source_forms(mol, template)
            self.assertIn(cs[t['source_form_compound_id']]['smiles'], {canonical(m) for m in forms})
            self.assertTrue(all(exact_scaffold_matches(m, template) for m in forms))
            if t['speciation_type'] == 'net-zero-intramolecular-proton-relocation':
                self.assertEqual(t['protons_consumed'], 0)
                self.assertEqual(len(t['left']), 1)
                self.assertEqual(len(t['right']), 1)

    def test_unassigned_or_inverted_glycerol_stereo_is_not_promoted(self):
        r = self.report
        t = next(t for t in r['targets'] if t['lipid_class'] == 'PE')
        cs = {c['id']: c for c in r['compounds']}
        src = next(s for s in r['source_records'] if s['rule_id'] == t['source_reaction_id'])
        template = next(m for m in map(Chem.MolFromSmiles, src['reaction_smarts'].split('>>')[1].split('.'))
                        if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
        mol = Chem.MolFromSmiles(cs[t['compound_id']]['smiles'])
        self.assertTrue(source_forms(mol, template))
        unassigned = Chem.Mol(mol)
        Chem.RemoveStereochemistry(unassigned)
        self.assertEqual(source_forms(unassigned, template), [])
        inverted = Chem.Mol(mol)
        for a in inverted.GetAtoms():
            if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                a.InvertChirality()
        self.assertEqual(source_forms(inverted, template), [])
