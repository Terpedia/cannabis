import json
import unittest
from pathlib import Path
from unittest.mock import patch
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_lipid_acylation_net import build


class SourceMappedProtonationTests(unittest.TestCase):
    def test_balance_source_gate_and_explicit_reversible_boundary(self):
        report = json.loads(Path('data/reports/phase1-source-mapped-protonation.json').read_text())
        joined = json.loads(Path('data/reports/phase1-protonation-verified-join.json').read_text())
        expected = {r['id'] for r in joined['rows'] if r['status'] == 'exact-endpoints-with-source-mapping'}
        self.assertEqual({r['source_bridge_id'] for r in report['reactions']}, expected)
        self.assertEqual(len(report['reactions']), 302)
        compounds = {c['id']: c for c in report['compounds']}
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            self.assertEqual(r['hypothesis_type'], 'source-mapped-protonation')
            self.assertEqual(r['direction_status'], 'explicit-reversible-acid-base-sensitivity-assumption')
            self.assertTrue(r['source_mapping_evidence'])
            self.assertFalse(r['enzyme_evidence_ids'])
            self.assertEqual({e['origin'] for e in r['source_mapping_evidence']}, {'computation'})
        compounds['test-co2'] = {'id': 'test-co2', 'smiles': 'O=C=O', 'carbon_count': 1}
        parent = {'targets': [], 'external_exchange_compound_ids': ['test-co2'], 'co2_compound_id': 'test-co2'}
        layer = {'compounds': [], 'added_reactions': [], 'forbidden_step_ids': ['inherited-exclusion']}
        module = 'cannabis_carbon.phase1_lipid_acylation_net.'
        with patch(module + 'assemble', return_value=({}, compounds, [])), \
             patch(module + 'NetModel', side_effect=RuntimeError('boundary')) as model:
            with self.assertRaisesRegex(RuntimeError, 'boundary'):
                build({'targets': []}, {}, parent, parent, report, prior_layers=(layer,))
        self.assertEqual(model.call_args.kwargs['forbidden_step_ids'], ['inherited-exclusion'])
