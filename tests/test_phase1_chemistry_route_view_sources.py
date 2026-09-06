import json
import unittest
from pathlib import Path

from cannabis_carbon.phase1_chemistry_route_view import reaction_sources


class ChemistryRouteSourceTests(unittest.TestCase):
    def test_amino_speciation_is_not_presented_as_exact_source_mapping(self):
        report = json.loads(Path('data/reports/phase1-amino-phospholipid-synthesis.json').read_text())
        reactions = [r for r in report['reactions'] if r['hypothesis_type'] == 'amino-phospholipid-speciation']
        self.assertEqual(len(reactions), 162)
        for reaction in reactions:
            source, = reaction_sources(reaction, {})
            self.assertEqual(source['speciation_type'], reaction['speciation_type'])
            self.assertEqual(source['claim_boundary'], reaction['claim_boundary'])
            self.assertIn('not-curated-reaction-or-exact-ChEBI-mapping', source['evidence_type'])
            self.assertNotIn('mapping_evidence', source)
            with self.assertRaisesRegex(ValueError, 'speciation classification'):
                reaction_sources({**reaction, 'speciation_type': None}, {})

    def test_all_protonation_sources_keep_mapping_provenance(self):
        report = json.loads(Path('data/reports/phase1-source-mapped-protonation.json').read_text())
        self.assertEqual(len(report['reactions']), 302)
        for reaction in report['reactions']:
            source, = reaction_sources(reaction, {})
            self.assertEqual(source['source_urls'], [reaction['source_url']])
            self.assertEqual(source['mapping_evidence'], reaction['source_mapping_evidence'])
            self.assertEqual(source['claim_boundary'], reaction['claim_boundary'])
            self.assertEqual(source['evidence_type'],
                             'exact-endpoint-protonation-mapping-not-curated-reaction')
            with self.assertRaisesRegex(ValueError, 'Missing exact protonation'):
                reaction_sources({**reaction, 'source_mapping_evidence': []}, {})

    def test_reaction_template_source_and_existing_evidence_are_preserved(self):
        existing = {'source_urls': ['https://example.org/existing']}
        reaction = {'source_reaction_id': 'RHEA:32932', 'sources': [existing]}
        index = {'RHEA:32932': {'source_url': 'https://www.rhea-db.org/rhea/32932'}}
        self.assertEqual(reaction_sources(reaction, index),
                         [existing, {'source_urls': ['https://www.rhea-db.org/rhea/32932']}])
        self.assertEqual(reaction['sources'], [existing])
