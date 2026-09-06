import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path

from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode, decode


class PhosphatidateHydrolysisNetTests(unittest.TestCase):
    def test_exact_certificates_replay_with_inherited_and_hydrolysis_bounds(self):
        path = Path('data/reports/phase1-phosphatidate-hydrolysis-net.json')
        report = json.loads(path.read_text())
        parent = json.loads(Path('data/reports/phase1-triglyceride-symmetry-net.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        self.assertEqual(report['external_exchange_compound_ids'], parent['external_exchange_compound_ids'])
        self.assertEqual([(t['cannabisdb_id'], t['compound_id']) for t in report['targets']],
                         [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']])
        compounds = {c['id']: c for c in report['compounds']}
        reactions = {r['id']: r for r in report['certificate_reactions'] + report['added_reactions']}
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']) |
                         {r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']})
        for r in report['added_reactions']:
            self.assertEqual(r['hypothesis_type'], 'phosphatidate-hydrolysis')
        for r in reactions.values():
            self.assertTrue(balanced([r['left'], r['right']], compounds))
        for cert in report['new_certificates']:
            self.assertFalse(forbidden & {s['step_id'] for s in cert['steps']})
            validate_certificate(cert, steps, compounds,
                                 set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        old_covered = {t['cannabisdb_id'] for t in parent['targets']
                       if t['net_status'] == 'exact-net-conversion-hypothesis'}
        covered = {t['cannabisdb_id'] for t in report['targets']
                   if t['net_status'] == 'exact-net-conversion-hypothesis'}
        self.assertTrue(old_covered <= covered)
        self.assertEqual(report['summary']['new_certificate_records'], len(covered - old_covered))
        self.assertEqual(report['summary']['net_status_counts'], dict(Counter(t['net_status'] for t in report['targets'])))
        self.assertEqual(sum(report['summary']['net_status_counts'].values()), 6220)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        self.assertEqual(decode(encode(report, sha)[::-1]), report)
