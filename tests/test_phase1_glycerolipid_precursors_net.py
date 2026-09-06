import copy
import json
import unittest
from unittest.mock import patch
from pathlib import Path
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build


class PrecursorMergeTests(unittest.TestCase):
    def test_exact_merge_and_identity_conflict_rejection(self):
        reports = [json.loads(Path('data/reports/phase1-' + name + '.json').read_text())
                   for name in ('glycerolipid-precursors', 'triglyceride-inventory-supplement')]
        merged = merge_hypotheses(*reports)
        self.assertEqual(len(merged['reactions']), 318)
        self.assertEqual({r['hypothesis_type'] for r in merged['reactions']},
                         {'sn1-acylation', 'sn2-acylation', 'sn3-acylation'})
        self.assertEqual(merge_hypotheses(*reports, *reports), merged)
        bad = copy.deepcopy(reports[0])
        bad['compounds'][0]['smiles'] = 'C'
        with self.assertRaisesRegex(ValueError, 'compound identity'):
            merge_hypotheses(reports[0], bad)
        bad = copy.deepcopy(reports[0])
        bad['reactions'][0]['left'][0]['coefficient'] += 1
        with self.assertRaises(ValueError):
            merge_hypotheses(bad)

    def test_all_three_acylation_types_are_forward_only_at_solver_boundary(self):
        reports = [json.loads(Path('data/reports/phase1-' + name + '.json').read_text())
                   for name in ('glycerolipid-precursors', 'triglyceride-inventory-supplement')]
        merged = merge_hypotheses(*reports)
        compounds = {c['id']: c for c in merged['compounds']}
        compounds['test-co2'] = {'id': 'test-co2', 'smiles': 'O=C=O', 'carbon_count': 1}
        parent = {'targets': [], 'external_exchange_compound_ids': ['test-co2'],
                  'co2_compound_id': 'test-co2'}
        module = 'cannabis_carbon.phase1_lipid_acylation_net.'
        with patch(module + 'assemble', return_value=({}, compounds, [])), \
             patch(module + 'NetModel', side_effect=RuntimeError('captured solver boundary')) as model:
            with self.assertRaisesRegex(RuntimeError, 'captured solver boundary'):
                build({'targets': []}, {}, parent, parent, merged)
        self.assertEqual(set(model.call_args.kwargs['forbidden_step_ids']),
                         {r['id'] + ':hypothetical-right-to-left' for r in merged['reactions']})
