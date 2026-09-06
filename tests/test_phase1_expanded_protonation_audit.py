import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_protonation_audit import fingerprint, bridge
from cannabis_carbon.phase1_row_export import encode, decode


class ExpandedProtonationAuditTests(unittest.TestCase):
    def test_exact_review_bridges_and_current_inventory(self):
        path = Path('data/reports/phase1-expanded-protonation-audit.json')
        report = json.loads(path.read_text())
        parent = json.loads(Path('data/reports/phase1-cardiolipin-net.json').read_text())
        for source, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(source).read_bytes()).hexdigest(), sha)
        self.assertEqual([(t['cannabisdb_id'], t['compound_id']) for t in report['targets']],
                         [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']])
        self.assertEqual(report['summary']['expanded_balanced_equations'], 17120)
        self.assertEqual(report['summary']['no_producer_records_with_review_bridges'], 381)
        compounds = {c['id']: c for c in report['compounds']}
        for row in report['bridges']:
            left, _ = fingerprint(compounds[row['target_compound_id']]['smiles'])
            right, _ = fingerprint(compounds[row['reaction_participant_compound_id']]['smiles'])
            self.assertEqual(left, row['target_identity_check'])
            self.assertEqual(right, row['participant_identity_check'])
            self.assertNotEqual(left['canonical_smiles'], right['canonical_smiles'])
            replay = bridge(left, right)
            for key, value in replay.items():
                self.assertEqual(row[key], value)
            self.assertEqual(row['status'], 'structurally-inferred-review-required')
            self.assertFalse(row['enzyme_evidence_ids'])
            self.assertTrue(row['source_reaction_participation'])
        self.assertEqual(decode(encode(report, hashlib.sha256(path.read_bytes()).hexdigest())), report)
