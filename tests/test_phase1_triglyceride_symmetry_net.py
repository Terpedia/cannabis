import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import decode


class SymmetryNetTests(unittest.TestCase):
    def test_exact_certificates_replay_with_all_inherited_exclusions(self):
        path = Path('data/reports/phase1-triglyceride-symmetry-net.json')
        report = json.loads(path.read_text())
        parent = json.loads(Path('data/reports/phase1-triglyceride-net.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in report['compounds']}
        reactions = {r['id']: r for r in report['certificate_reactions'] + report['added_reactions']}
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']) | {r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']})
        for r in reactions.values():
            self.assertTrue(balanced([r['left'], r['right']], compounds))
        for cert in report['new_certificates']:
            self.assertFalse(forbidden & {s['step_id'] for s in cert['steps']})
            validate_certificate(cert, steps, compounds, set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        self.assertEqual(len(report['new_certificates']), 7)
        self.assertEqual(report['summary']['net_status_counts']['exact-net-conversion-hypothesis'], 477)
        self.assertEqual(sum(report['summary']['net_status_counts'].values()), 6220)
        rows = [json.loads(line) for line in Path('data/derived/phase1-triglyceride-symmetry-net.ndjson').read_text().splitlines()]
        self.assertEqual(decode(rows[::-1]), report)
        self.assertEqual({r['report_sha256'] for r in rows}, {hashlib.sha256(path.read_bytes()).hexdigest()})
