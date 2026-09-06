import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_cardiolipin_identity_audit import build
from cannabis_carbon.phase1_lipid_acylation import canonical


class CardiolipinIdentityAuditTests(unittest.TestCase):
    def test_full_inventory_identity_review_does_not_assign_stereo_or_reactions(self):
        report = json.loads(Path('data/reports/phase1-cardiolipin-identity-audit.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        network = json.loads(Path('data/reports/phase1-full-balanced-network.json').read_text())
        catalog = json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text())
        for t in network['targets']:
            t['label'] = 'label not used in matching'
        rebuilt = build(network, catalog)
        self.assertEqual(rebuilt['summary'], report['summary'])
        self.assertEqual(rebuilt['summary']['central_glycerol_stereo_counts'],
                         {'Unspecified': 1273, 'Specified': 92})
        self.assertEqual(rebuilt['summary']['inventory_target_records_scanned'], 6220)
        self.assertEqual({t['cannabisdb_id'] for t in rebuilt['targets']},
                         {t['cannabisdb_id'] for t in report['targets']})
        for t in report['targets']:
            self.assertEqual(canonical(Chem.MolFromSmiles(t['source_smiles'])), t['canonical_smiles'])
            self.assertEqual(Chem.GetFormalCharge(Chem.MolFromSmiles(t['canonical_smiles'])), 0)
        self.assertNotIn('reactions', report)
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)
