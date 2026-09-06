import csv
import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_protonation_source_join import exact_key


class VerifiedJoinTests(unittest.TestCase):
    def test_all_matched_endpoints_and_mapping_origins_replay(self):
        report = json.loads(Path('data/reports/phase1-protonation-verified-join.json').read_text())
        audit = json.loads(Path('data/reports/phase1-expanded-protonation-audit.json').read_text())
        for source, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(source).read_bytes()).hexdigest(), sha)
        with Path('data/raw/rhea-chebi-smiles-20260906.tsv').open() as handle:
            structures = dict(csv.reader(handle, delimiter='\t'))
        for check in report['retrieved_identity_checks']:
            if check['status'] == 'exact-id-structure-available':
                data = json.loads(Path(check['source_path']).read_text())
                self.assertEqual(data['chebi_accession'], check['chebi_id'])
                structures[check['chebi_id']] = data['default_structure']['smiles']
        with Path('data/raw/rhea-chebi-ph73-mapping-20260906.tsv').open() as handle:
            mappings = {('CHEBI:' + r['CHEBI'], 'CHEBI:' + r['CHEBI_PH7_3'], r['ORIGIN'])
                        for r in csv.DictReader(handle, delimiter='\t')}
        compounds = {c['id']: c['smiles'] for c in audit['compounds']}
        self.assertEqual([r['id'] for r in report['rows']], [r['id'] for r in audit['bridges']])
        matched = 0
        for row in report['rows']:
            for field, cid in [('target_exact_chebi_ids', row['target_compound_id']),
                               ('participant_exact_chebi_ids', row['reaction_participant_compound_id'])]:
                for chebi in row[field]:
                    self.assertEqual(exact_key(structures[chebi]), exact_key(compounds[cid]))
            for evidence in row['mapping_evidence']:
                self.assertIn((evidence['source_chebi_id'], evidence['major_chebi_id'], evidence['origin']), mappings)
                a, b = evidence['source_chebi_id'], evidence['major_chebi_id']
                self.assertTrue((a in row['target_exact_chebi_ids'] and b in row['participant_exact_chebi_ids']) or
                                (b in row['target_exact_chebi_ids'] and a in row['participant_exact_chebi_ids']))
            if row['status'] == 'exact-endpoints-with-source-mapping':
                matched += 1
                self.assertTrue(row['mapping_evidence'])
        self.assertEqual(matched, 302)
        self.assertEqual(len(report['retrieved_identity_checks']), 393)
