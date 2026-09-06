"""Audit the completed sensitivity run; a pending run is not verified coverage."""
import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path

from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import decode, encode
from cannabis_carbon.phase1_scope import orientations


class SourceMappedProtonationNetTests(unittest.TestCase):
    def test_full_inventory_and_every_new_certificate(self):
        path = Path('data/reports/phase1-source-mapped-protonation-net.json')
        if not path.exists():
            self.skipTest('Sensitivity calculation pending; no new coverage verified')
        report = json.loads(path.read_bytes())
        parent = json.loads(Path('data/reports/phase1-cardiolipin-net.json').read_bytes())
        evidence = json.loads(Path('data/reports/phase1-source-mapped-protonation.json').read_bytes())
        for source, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(source).read_bytes()).hexdigest(), sha)
        self.assertEqual(report['schema'], 'cannabis-carbon.phase1-source-mapped-protonation-net.v1')
        self.assertIn('separate full-inventory sensitivity scenario', report['claim_boundary'])
        for key in ('external_exchange_compound_ids', 'co2_compound_id'):
            self.assertEqual(report[key], parent[key])
        forbidden = set(report['forbidden_step_ids'])
        self.assertEqual(forbidden, set(parent['forbidden_step_ids']))
        self.assertEqual([(t['cannabisdb_id'], t['compound_id']) for t in report['targets']],
                         [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']])
        self.assertEqual(len(report['targets']), 6220)
        self.assertEqual(len({t['cannabisdb_id'] for t in report['targets']}), 6220)
        source_reactions = {r['id']: r for r in evidence['reactions']}
        compounds = {c['id']: c for c in report['compounds']}
        carbon_exchanges = {cid for cid in report['external_exchange_compound_ids']
                            if compounds[cid]['carbon_count']}
        self.assertEqual(carbon_exchanges, {report['co2_compound_id']})
        self.assertEqual(compounds[report['co2_compound_id']]['smiles'], 'O=C=O')
        reactions = {r['id']: r for r in report['certificate_reactions'] + report['added_reactions']}
        steps = {s['id']: s for s in orientations(list(reactions.values()))}
        added = {r['id'] for r in report['added_reactions']}
        for reaction in report['added_reactions']:
            source = source_reactions[reaction['id']]
            for key in ('left', 'right', 'hypothesis_type', 'source_mapping_evidence',
                        'direction_status', 'source_bridge_id', 'claim_boundary'):
                self.assertEqual(reaction[key], source[key])
            self.assertEqual(reaction['hypothesis_type'], 'source-mapped-protonation')
            self.assertTrue(reaction['source_mapping_evidence'])
            self.assertFalse(reaction['enzyme_evidence_ids'])
        for reaction in reactions.values():
            self.assertTrue(balanced([reaction['left'], reaction['right']], compounds), reaction['id'])
        for certificate in report['new_certificates']:
            used = {s['step_id'] for s in certificate['steps']}
            self.assertFalse(used & forbidden)
            self.assertTrue(any(s.startswith(r + ':') for s in used for r in added))
            validate_certificate(certificate, steps, compounds,
                                 set(report['external_exchange_compound_ids']), report['co2_compound_id'])
        old = {t['cannabisdb_id'] for t in parent['targets']
               if t['net_status'] == 'exact-net-conversion-hypothesis'}
        covered = {t['cannabisdb_id'] for t in report['targets']
                   if t['net_status'] == 'exact-net-conversion-hypothesis'}
        self.assertTrue(old <= covered)
        new_structures = {t['compound_id'] for t in report['targets']
                          if t['cannabisdb_id'] in covered - old}
        self.assertEqual({c['compound_id'] for c in report['new_certificates']}, new_structures)
        self.assertEqual(len(report['new_certificates']), len(new_structures))
        self.assertEqual(report['summary']['new_certificate_structures'], len(new_structures))
        self.assertEqual(report['summary']['new_certificate_records'], len(covered - old))
        self.assertEqual(report['summary']['net_status_counts'],
                         dict(Counter(t['net_status'] for t in report['targets'])))
        self.assertEqual(decode(encode(report, hashlib.sha256(path.read_bytes()).hexdigest())[::-1]), report)
