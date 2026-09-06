import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_triglyceride_inventory_supplement import build
from cannabis_carbon.phase1_marts_completions import balanced
from test_phase1_triglyceride_symmetry import forward_products


class TriglycerideInventorySupplementTests(unittest.TestCase):
    def test_full_inventory_not_label_prefix_selects_exact_additional_targets(self):
        names = ('full-balanced-network', 'triglyceride-acylation',
                 'triglyceride-symmetry', 'triglyceride-symmetry-net')
        inputs = [json.loads(Path('data/reports/phase1-' + n + '.json').read_text()) for n in names]
        report = json.loads(Path('data/reports/phase1-triglyceride-inventory-supplement.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        for t in inputs[0]['targets']:
            t['label'] = 'label deliberately removed'
        rebuilt = build(*inputs)
        self.assertEqual(rebuilt['reactions'], report['reactions'])
        self.assertEqual({t['cannabisdb_id'] for t in report['targets']},
                         {'CDB000794', 'CDB000795', 'CDB000796', 'CDB000798', 'CDB000799'})
        self.assertEqual(report['summary']['inventory_target_records'],
                         report['summary']['previously_proposed_target_records'] + len(report['inventory_scan']))
        compounds = {c['id']: c for c in report['compounds']}
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            self.assertEqual(r['enzyme_evidence_ids'], [])
            sides = [[compounds[p['compound_id']]['smiles'] for p in r[side]
                      for _ in range(p['coefficient'])] for side in ('left', 'right')]
            self.assertEqual(forward_products(sides[0]), sorted(sides[1]))
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)
