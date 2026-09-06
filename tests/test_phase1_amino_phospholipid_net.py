import json
import hashlib
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch
from cannabis_carbon.phase1_amino_phospholipid_net import TYPES
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode, decode


class AminoPhospholipidNetBoundaryTests(unittest.TestCase):
    def test_completed_inventory_replays_every_new_certificate(self):
        path = Path('data/reports/phase1-amino-phospholipid-net.json')
        if not path.exists():
            self.skipTest('Full calculation pending; no new net coverage verified')
        report = json.loads(path.read_bytes())
        parent = json.loads(Path('data/reports/phase1-source-mapped-protonation-net.json').read_bytes())
        self.assertEqual(report['schema'], 'cannabis-carbon.phase1-amino-phospholipid-net.v1')
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        self.assertEqual([(t['cannabisdb_id'], t['compound_id']) for t in report['targets']],
                         [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']])
        self.assertEqual(len(report['targets']), 6220)
        self.assertEqual(len({t['cannabisdb_id'] for t in report['targets']}), 6220)
        for k in ('co2_compound_id', 'external_exchange_compound_ids'):
            self.assertEqual(report[k], parent[k])
        cs = {c['id']: c for c in report['compounds']}
        self.assertEqual({cid for cid in report['external_exchange_compound_ids'] if cs[cid]['carbon_count']},
                         {report['co2_compound_id']})
        self.assertEqual(cs[report['co2_compound_id']]['smiles'], 'O=C=O')
        reactions = {r['id']: r for r in report['added_reactions'] + report['certificate_reactions']}
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']) | {
            r['id'] + ':hypothetical-right-to-left' for r in report['added_reactions']
            if r['hypothesis_type'] != 'amino-phospholipid-speciation'})
        evidence = [json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
                    for n in ('amino-phospholipid-synthesis', 'amino-phospholipid-precursors')]
        sources = {r['id']: r for e in evidence for r in e['reactions']}
        added = {r['id'] for r in report['added_reactions']}
        for r in report['added_reactions']:
            self.assertIn(r['hypothesis_type'], TYPES)
            self.assertEqual(r, sources[r['id']])
        for r in reactions.values():
            self.assertTrue(balanced([r['left'], r['right']], cs), r['id'])
        for cert in report['new_certificates']:
            self.assertFalse(forbidden & {s['step_id'] for s in cert['steps']})
            self.assertTrue(added & {s['reaction_id'] for s in cert['steps']})
            validate_certificate(cert, steps, cs, set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        old = {t['cannabisdb_id'] for t in parent['targets'] if t['net_status'] == 'exact-net-conversion-hypothesis'}
        covered = {t['cannabisdb_id'] for t in report['targets'] if t['net_status'] == 'exact-net-conversion-hypothesis'}
        self.assertTrue(old <= covered)
        new_structures = {t['compound_id'] for t in report['targets'] if t['cannabisdb_id'] in covered - old}
        self.assertEqual({c['compound_id'] for c in report['new_certificates']}, new_structures)
        self.assertEqual(len(report['new_certificates']), len(new_structures))
        self.assertEqual(report['summary']['new_certificate_records'], len(covered - old))
        self.assertEqual(report['summary']['new_certificate_structures'], len(new_structures))
        self.assertEqual(report['summary']['net_status_counts'], dict(Counter(t['net_status'] for t in report['targets'])))
        self.assertEqual(decode(encode(report, hashlib.sha256(path.read_bytes()).hexdigest())[::-1]), report)

    def test_new_synthesis_cannot_reverse_and_only_speciation_is_reversible(self):
        reports = [json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
                   for n in ('amino-phospholipid-synthesis', 'amino-phospholipid-precursors')]
        merged = merge_hypotheses(*reports, allowed_types=TYPES)
        compounds = {c['id']: c for c in merged['compounds']}
        compounds['co2'] = {'id': 'co2', 'smiles': 'O=C=O', 'carbon_count': 1}
        parent = {'targets': [], 'external_exchange_compound_ids': ['co2'], 'co2_compound_id': 'co2'}
        prior = {'compounds': [], 'added_reactions': [], 'forbidden_step_ids': ['inherited:reverse']}
        prefix = 'cannabis_carbon.phase1_lipid_acylation_net.'
        with patch(prefix + 'assemble', return_value=({}, compounds, [])), \
             patch(prefix + 'NetModel', side_effect=RuntimeError('captured')) as model:
            with self.assertRaisesRegex(RuntimeError, 'captured'):
                build({'targets': []}, {}, parent, parent, merged, prior_layers=(prior,))
        # Equation joins share the first equation's direction constraints.
        from cannabis_carbon.phase1_lipid_acylation_net import equation_key
        unique = {}
        for r in merged['reactions']:
            unique.setdefault(equation_key(r), r)
        expected = {'inherited:reverse'} | {r['id'] + ':hypothetical-right-to-left'
            for r in unique.values() if r['hypothesis_type'] != 'amino-phospholipid-speciation'}
        self.assertEqual(set(model.call_args.kwargs['forbidden_step_ids']), expected)
