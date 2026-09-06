import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


class TriglycerideNetTests(unittest.TestCase):
    def test_complete_batch_exact_replay_and_inherited_direction_restrictions(self):
        report = json.loads(Path('data/reports/phase1-triglyceride-net.json').read_text())
        parent = json.loads(Path('data/reports/phase1-lipid-acylation-net.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in report['compounds']}
        reactions = {r['id']: r for r in report['certificate_reactions'] + report['added_reactions']}
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']) |
            {r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']})
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        for r in reactions.values():
            self.assertTrue(balanced([r['left'], r['right']], compounds))
        for cert in report['new_certificates']:
            self.assertFalse(forbidden & {s['step_id'] for s in cert['steps']})
            validate_certificate(cert, steps, compounds, set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        self.assertEqual(len(report['new_certificates']), 111)
        self.assertEqual(report['summary']['net_status_counts']['exact-net-conversion-hypothesis'], 470)
        self.assertEqual(sum(report['summary']['net_status_counts'].values()), 6220)
        for t in report['targets']:
            if t['parent_net_status'] == 'exact-net-conversion-hypothesis':
                self.assertEqual(t['net_status'], t['parent_net_status'])
