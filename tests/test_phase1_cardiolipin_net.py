import json
import hashlib
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch
from cannabis_carbon.phase1_cardiolipin_net import TYPES
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode, decode


class CardiolipinNetBoundaryTests(unittest.TestCase):
    def test_completed_report_replays_all_new_certificates_and_preserves_inventory(self):
        path = Path('data/reports/phase1-cardiolipin-net.json')
        if not path.exists():
            self.skipTest('Full-inventory cardiolipin report still pending; no new coverage verified')
        report = json.loads(path.read_text())
        parent = json.loads(Path('data/reports/phase1-glycerolipid-precursors-net.json').read_text())
        for source, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(source).read_bytes()).hexdigest(), sha)
        self.assertEqual(report['external_exchange_compound_ids'], parent['external_exchange_compound_ids'])
        self.assertEqual(report['co2_compound_id'], parent['co2_compound_id'])
        self.assertEqual([(t['cannabisdb_id'], t['compound_id']) for t in report['targets']],
                         [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']])
        compounds = {c['id']: c for c in report['compounds']}
        reactions = {r['id']: r for r in report['certificate_reactions'] + report['added_reactions']}
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']) |
                         {r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']
                          if r['hypothesis_type'] != 'cardiolipin-protonation'})
        for reaction in report['added_reactions']:
            self.assertIn(reaction['hypothesis_type'], TYPES)
        for reaction in reactions.values():
            self.assertTrue(balanced([reaction['left'], reaction['right']], compounds))
        added = {r['id'] for r in report['added_reactions']}
        for cert in report['new_certificates']:
            cert_steps = {s['step_id'] for s in cert['steps']}
            self.assertFalse(forbidden & cert_steps)
            self.assertTrue(any(step.startswith(rid + ':') for step in cert_steps for rid in added))
            validate_certificate(cert, steps, compounds,
                                 set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        old = {t['cannabisdb_id'] for t in parent['targets']
               if t['net_status'] == 'exact-net-conversion-hypothesis'}
        covered = {t['cannabisdb_id'] for t in report['targets']
                   if t['net_status'] == 'exact-net-conversion-hypothesis'}
        self.assertTrue(old <= covered)
        self.assertEqual(report['summary']['new_certificate_records'], len(covered - old))
        self.assertEqual(report['summary']['net_status_counts'], dict(Counter(t['net_status'] for t in report['targets'])))
        self.assertEqual(sum(report['summary']['net_status_counts'].values()), 6220)
        self.assertEqual(decode(encode(report, hashlib.sha256(path.read_bytes()).hexdigest())[::-1]), report)

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
