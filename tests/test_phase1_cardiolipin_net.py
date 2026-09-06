import json
import unittest
from pathlib import Path
from unittest.mock import patch
from cannabis_carbon.phase1_cardiolipin_net import TYPES
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build


class CardiolipinNetBoundaryTests(unittest.TestCase):
    def test_full_package_balances_and_only_explicit_equilibria_can_reverse(self):
        reports = [json.loads(Path('data/reports/phase1-' + n + '.json').read_text())
                   for n in ('cardiolipin-synthesis', 'cardiolipin-precursors')]
        merged = merge_hypotheses(*reports, allowed_types=TYPES)
        self.assertEqual(len(merged['reactions']), 305)
        self.assertEqual(merge_hypotheses(*reports, *reports, allowed_types=TYPES), merged)
        compounds = {c['id']: c for c in merged['compounds']}
        compounds['test-co2'] = {'id': 'test-co2', 'smiles': 'O=C=O', 'carbon_count': 1}
        parent = {'targets': [], 'external_exchange_compound_ids': ['test-co2'], 'co2_compound_id': 'test-co2'}
        prefix = 'cannabis_carbon.phase1_lipid_acylation_net.'
        with patch(prefix + 'assemble', return_value=({}, compounds, [])), \
             patch(prefix + 'NetModel', side_effect=RuntimeError('captured boundary')) as model:
            with self.assertRaisesRegex(RuntimeError, 'captured boundary'):
                build({'targets': []}, {}, parent, parent, merged)
        expected = {r['id'] + ':hypothetical-right-to-left' for r in merged['reactions']
                    if r['hypothesis_type'] != 'cardiolipin-protonation'}
        self.assertEqual(len(expected), 213)
        self.assertEqual(set(model.call_args.kwargs['forbidden_step_ids']), expected)
