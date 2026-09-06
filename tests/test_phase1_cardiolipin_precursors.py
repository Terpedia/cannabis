import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced
from test_phase1_cardiolipin_synthesis import forward_products as phosphatidyl_transfer


def hydrolyze(reactants):
    molecules = [Chem.MolFromSmiles(s) for s in reactants]
    pgp = next(m for m in molecules if m.GetNumAtoms() > 1)
    water = next(m for m in molecules if m.GetNumAtoms() == 1)
    assert canonical(water) == 'O'
    phosphorus = [a for a in pgp.GetAtoms() if a.GetAtomicNum() == 15
                  and sum(n.GetFormalCharge() == -1 for n in a.GetNeighbors()) == 2]
    assert len(phosphorus) == 1
    p = phosphorus[0]
    ester = [a for a in p.GetNeighbors() if a.GetAtomicNum() == 8
             and any(n.GetAtomicNum() == 6 for n in a.GetNeighbors())]
    assert len(ester) == 1
    edit = Chem.RWMol(Chem.CombineMols(pgp, water))
    edit.RemoveBond(p.GetIdx(), ester[0].GetIdx())
    edit.AddBond(p.GetIdx(), pgp.GetNumAtoms(), Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted(canonical(m) for m in Chem.GetMolFrags(edit, asMols=True))


def cdp_synthesis(reactants):
    assert reactants.count('[H+]') == 1
    molecules = [Chem.MolFromSmiles(s) for s in reactants if s != '[H+]']
    ctp = next(m for m in molecules if any(a.GetAtomicNum() == 7 for a in m.GetAtoms()))
    pa = next(m for m in molecules if m is not ctp)
    alpha = [a for a in ctp.GetAtoms() if a.GetAtomicNum() == 15
             and any(n.GetAtomicNum() == 8 and any(c.GetAtomicNum() == 6 for c in n.GetNeighbors())
                     for n in a.GetNeighbors())]
    assert len(alpha) == 1
    bridge = [a for a in alpha[0].GetNeighbors() if a.GetAtomicNum() == 8
              and a.GetDegree() == 2 and all(n.GetAtomicNum() == 15 for n in a.GetNeighbors())]
    nucleophiles = [a for a in pa.GetAtoms() if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
    assert len(bridge) == 1 and len(nucleophiles) == 2
    edit = Chem.RWMol(Chem.CombineMols(pa, ctp)); offset = pa.GetNumAtoms()
    edit.RemoveBond(offset + alpha[0].GetIdx(), offset + bridge[0].GetIdx())
    edit.GetAtomWithIdx(offset + bridge[0].GetIdx()).SetNoImplicit(False)
    edit.GetAtomWithIdx(nucleophiles[0].GetIdx()).SetFormalCharge(0)
    edit.AddBond(nucleophiles[0].GetIdx(), offset + alpha[0].GetIdx(), Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted(canonical(m) for m in Chem.GetMolFrags(edit, asMols=True))


class CardiolipinPrecursorTests(unittest.TestCase):
    def test_every_equation_balances_and_reconstructs_exact_products(self):
        report = json.loads(Path('data/reports/phase1-cardiolipin-precursors.json').read_text())
        for p, sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(p).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in report['compounds']}
        replay = {'pgp-hydrolysis': hydrolyze, 'pgp-synthesis': phosphatidyl_transfer,
                  'cdp-dag-synthesis': cdp_synthesis}
        for r in report['reactions']:
            self.assertTrue(balanced([r['left'], r['right']], compounds))
            sides = [[compounds[p['compound_id']]['smiles'] for p in r[side]
                      for _ in range(p['coefficient'])] for side in ('left', 'right')]
            self.assertEqual(replay[r['hypothesis_type']](sides[0]), sorted(sides[1]))
            self.assertEqual(r['direction_status'], 'forward-only-generic-source-hypothesis')
            self.assertEqual(r['enzyme_evidence_ids'], [])
        self.assertEqual(report['summary']['balanced_equations'], 121)
        cardiolipin = json.loads(Path('data/reports/phase1-cardiolipin-synthesis.json').read_text())
        missing = {c for t in cardiolipin['targets'] for c in t['precursors_absent_from_parent']}
        self.assertEqual(len(missing), 36)
        self.assertTrue(missing <= {r['compound_id'] for r in report['precursor_candidates']})
        self.assertEqual(report['summary']['new_CO2_route_claims'], 0)
