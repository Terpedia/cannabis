import unittest
from cannabis_carbon.phase1_protonation_source_join import build, exact_key


class SourceJoinTests(unittest.TestCase):
    def test_incoming_mapping_is_not_an_exact_target_assignment(self):
        audit = {'compounds': [{'id': 'acid', 'smiles': 'CC(=O)O'},
                              {'id': 'base', 'smiles': 'CC(=O)[O-]'}],
                 'bridges': [{'id': 'b', 'target_compound_id': 'acid',
                              'reaction_participant_compound_id': 'base', 'cannabisdb_ids': ['t']}]}
        mappings = [{'CHEBI': '1', 'CHEBI_PH7_3': '2', 'ORIGIN': 'curation'}]
        result = build(audit, [('CHEBI:2', 'CC(=O)[O-]')], mappings)['rows'][0]
        self.assertEqual(result['status'], 'missing-exact-target-endpoint')
        self.assertFalse(result['mapping_evidence'])
        self.assertFalse(result['target_exact_chebi_ids'])
        self.assertEqual(result['unresolved_source_structure_leads'][0]['source_chebi_id'], 'CHEBI:1')
        result = build(audit, [('CHEBI:1', 'CC(=O)O'), ('CHEBI:2', 'CC(=O)[O-]')], mappings)['rows'][0]
        self.assertEqual(result['status'], 'exact-endpoints-with-source-mapping')
        self.assertEqual(result['mapping_evidence'][0]['origin'], 'curation')

    def test_stereo_charge_isotope_and_fragments_remain_distinct(self):
        for left, right in [('CC(=O)O', 'CC(=O)[O-]'), ('C', '[13CH4]'),
                            ('C[C@H](O)C(=O)O', 'C[C@@H](O)C(=O)O'),
                            ('C[C@H](O)C(=O)O', 'CC(O)C(=O)O'), ('CCO', 'CCO.[Na+]')]:
            self.assertNotEqual(exact_key(left), exact_key(right))
        self.assertIsNone(exact_key('*CC'))
        self.assertEqual(exact_key('OC(C)=O'), exact_key('CC(=O)O'))
