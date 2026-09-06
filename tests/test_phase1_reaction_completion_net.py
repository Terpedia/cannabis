import copy
import hashlib
import json
import unittest
from pathlib import Path

from cannabis_carbon.phase1_balance_reference import concrete_participants
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_reaction_completion_net import assemble, build, validate_certificate
from cannabis_carbon.phase1_scope import orientations


class ReactionCompletionNetTests(unittest.TestCase):
    def fixture(self):
        compounds, sides = {}, []
        for side in concrete_participants('O=C=O.O>>O=C(O)O'):
            members = []
            for p in side:
                cid = stable_id('structure', p['smiles'])
                compounds[cid] = {'id': cid, **{k: v for k, v in p.items() if k != 'coefficient'}}
                members.append({'compound_id': cid, 'coefficient': p['coefficient']})
            sides.append(members)
        target = sides[1][0]['compound_id']
        co2 = next(c for c in compounds if compounds[c]['smiles'] == 'O=C=O')
        exchange = [p['compound_id'] for p in sides[0]]
        row = {'cannabisdb_id': 'fixture-only', 'compound_id': target, 'label': 'carbonic acid'}
        network = {'compounds': list(compounds.values()), 'reactions': [], 'targets': [row]}
        h = {'id': 'fixture-completion', 'balanced_equation_id': stable_id('balanced-equation', sides),
             'left': sides[0], 'right': sides[1], 'enzyme_evidence_ids': []}
        completions = {'compounds': list(compounds.values()), 'completions': [h]}
        baseline = {'targets': [{**row, 'net_status': 'no-net-producing-catalog-equation'}],
                    'certificates': [], 'external_exchange_compound_ids': exchange, 'co2_compound_id': co2}
        return network, completions, baseline

    def test_balanced_hypothesis_does_not_require_enzyme(self):
        inputs = self.fixture()
        original = copy.deepcopy(inputs)
        report = build(*inputs)
        self.assertEqual(inputs, original)
        self.assertEqual(report['summary']['new_certificate_records'], 1)
        self.assertEqual(report['summary']['added_completion_equations'], 1)
        self.assertEqual(report['targets'][0]['startup_status'], 'structural-scope-reachable')
        self.assertEqual(len(report['new_certificates'][0]['added_reaction_ids']), 1)
        self.assertEqual(report['added_reactions'][0]['enzyme_evidence_ids'], [])

    def test_balance_failure_even_with_recomputed_id(self):
        network, completions, _ = self.fixture()
        h = completions['completions'][0]
        h['right'][0]['coefficient'] = 2
        h['balanced_equation_id'] = stable_id('balanced-equation', [h['left'], h['right']])
        with self.assertRaisesRegex(ValueError, 'balance'):
            assemble(network, completions)

    def test_duplicate_hypothesis_does_not_duplicate_equation(self):
        network, completions, _ = self.fixture()
        h = copy.deepcopy(completions['completions'][0]); h['id'] = 'second-reference'
        completions['completions'].append(h)
        reactions, _, added = assemble(network, completions)
        self.assertEqual(len(reactions), 1)
        self.assertEqual(len(next(iter(added.values()))['completion_ids']), 2)

    def test_organic_exchange_rejected(self):
        network, completions, baseline = self.fixture()
        baseline['external_exchange_compound_ids'].append(network['targets'][0]['compound_id'])
        with self.assertRaisesRegex(ValueError, 'sole carbon'):
            build(network, completions, baseline)

    def test_target_identity_mismatch_rejected(self):
        network, completions, baseline = self.fixture()
        baseline['targets'][0]['compound_id'] = 'different-stereoisomer'
        with self.assertRaisesRegex(ValueError, 'inventory'):
            build(network, completions, baseline)

    def test_published_certificates_replay_and_source_hashes(self):
        path = Path('data/reports/phase1-reaction-completion-net.json')
        report = json.loads(path.read_text())
        self.assertEqual(path.read_bytes(), Path('docs/data/reaction-completion-net.json').read_bytes())
        for name, digest in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(name).read_bytes()).hexdigest(), digest)
        steps = {s['id']: s for s in orientations(report['certificate_reactions'])}
        compounds = {c['id']: c for c in report['compounds']}
        for cert in report['new_certificates']:
            validate_certificate(cert, steps, compounds,
                set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        self.assertEqual(len(report['targets']), 6220)
        self.assertEqual(len(report['new_certificates']), 10)
        self.assertEqual(sum(t['new_certificate'] for t in report['targets']), 10)
