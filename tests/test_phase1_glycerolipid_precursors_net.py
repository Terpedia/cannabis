import copy
import hashlib
import json
import unittest
from collections import Counter
from unittest.mock import patch
from pathlib import Path
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode, decode


class PrecursorMergeTests(unittest.TestCase):
    def test_completed_report_exact_replay_and_whole_inventory(self):
        path = Path('data/reports/phase1-glycerolipid-precursors-net.json')
        if not path.exists():
            self.skipTest('Full-inventory solver report not yet available; coverage remains unverified')
        report = json.loads(path.read_text())
        parent = json.loads(Path('data/reports/phase1-phosphatidate-hydrolysis-net.json').read_text())
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
                         {r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']})
        added = {r['id'] for r in report['added_reactions']}
        for reaction in report['added_reactions']:
            self.assertIn(reaction['hypothesis_type'], ('sn1-acylation', 'sn2-acylation', 'sn3-acylation'))
        for reaction in reactions.values():
            self.assertTrue(balanced([reaction['left'], reaction['right']], compounds))
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
