import hashlib
import json
import unittest
from pathlib import Path
from cannabis_carbon.phase1_marts_completions import balanced
from test_phase1_amino_phospholipid_synthesis import forward_transfer
from test_phase1_cardiolipin_precursors import cdp_synthesis
from test_phase1_glycerolipid_precursors import forward_products as acylation


class AminoPhospholipidPrecursorTests(unittest.TestCase):
    def test_all_precursor_products_replay_and_frontier_is_explicit(self):
        r = json.loads(Path('data/reports/phase1-amino-phospholipid-precursors.json').read_bytes())
        s = json.loads(Path('data/reports/phase1-amino-phospholipid-synthesis.json').read_bytes())
        for path, sha in r['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), sha)
        cs = {c['id']: c for c in r['compounds']}
        self.assertEqual(len(r['reactions']), 177)
        for reaction in r['reactions']:
            self.assertTrue(balanced([reaction['left'], reaction['right']], cs))
            sides = [[cs[p['compound_id']]['smiles'] for p in reaction[side]
                      for _ in range(p['coefficient'])] for side in ('left', 'right')]
            kind = reaction['hypothesis_type']
            products = (cdp_synthesis(sides[0]) if kind == 'cdp-dag-synthesis' else
                        acylation(sides[0], kind) if kind.startswith('sn') else
                        forward_transfer(sides[0], 'methylation' in kind))
            self.assertEqual(products, sorted(sides[1]))
            self.assertEqual(reaction['direction_status'], 'forward-only-generic-source-hypothesis')
            self.assertFalse(reaction['enzyme_evidence_ids'])
        missing = {p for t in s['targets'] for p in t['precursors_absent_from_parent']}
        proposed = {p['compound_id'] for p in r['precursor_candidates']}
        self.assertEqual(len(missing), 56)
        self.assertTrue(missing <= proposed)
        frontier = {p['compound_id'] for p in r['frontier']}
        roots = {p for t in s['targets'] for p in t['required_precursor_ids']}
        required = {p for row in r['precursor_candidates'] for p in row['required_precursor_ids']}
        self.assertEqual(frontier, (roots | required) - proposed)
        self.assertEqual(len(frontier), 101)
        self.assertTrue(all(p['in_parent_inventory'] for p in r['frontier']))
        self.assertEqual(r['summary']['new_CO2_route_claims'], 0)
