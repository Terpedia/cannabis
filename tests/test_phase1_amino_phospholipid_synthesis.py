import hashlib
import json
import unittest
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_amino_phospholipid_synthesis import instantiate
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced


def forward_transfer(smiles, methylation):
    molecules = [Chem.MolFromSmiles(s) for s in smiles]
    if methylation:
        donor = next(m for m in molecules if any(a.GetAtomicNum() == 16 for a in m.GetAtoms()))
        acceptor = next(m for m in molecules if m is not donor)
        sulfur, = [a for a in donor.GetAtoms() if a.GetAtomicNum() == 16 and a.GetFormalCharge() == 1]
        carbon, = [a for a in sulfur.GetNeighbors() if a.GetAtomicNum() == 6 and a.GetDegree() == 1]
        nitrogen, = [a for a in acceptor.GetAtoms() if a.GetAtomicNum() == 7 and a.GetFormalCharge() == 1]
        edit = Chem.RWMol(Chem.CombineMols(donor, acceptor)); offset = donor.GetNumAtoms()
        edit.RemoveBond(sulfur.GetIdx(), carbon.GetIdx())
        edit.GetAtomWithIdx(sulfur.GetIdx()).SetFormalCharge(0)
        n = edit.GetAtomWithIdx(offset + nitrogen.GetIdx())
        n.SetNumExplicitHs(n.GetTotalNumHs() - 1)
        edit.AddBond(carbon.GetIdx(), n.GetIdx(), Chem.BondType.SINGLE)
    else:
        donor = next(m for m in molecules if sum(a.GetAtomicNum() == 15 for a in m.GetAtoms()) == 2)
        acceptor = next(m for m in molecules if m is not donor)
        bridge, = [a for a in donor.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 2
                   and all(n.GetAtomicNum() == 15 for n in a.GetNeighbors())]
        candidates = []
        for p in bridge.GetNeighbors():
            cut = Chem.RWMol(donor); cut.RemoveBond(p.GetIdx(), bridge.GetIdx())
            fragment = next(f for f in Chem.GetMolFrags(cut) if p.GetIdx() in f)
            if not any(donor.GetAtomWithIdx(i).GetAtomicNum() == 7 and donor.GetAtomWithIdx(i).IsInRing() for i in fragment):
                candidates.append(p.GetIdx())
        phosphorus, = candidates
        oxygen, = [a.GetIdx() for a in acceptor.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 1
                   and a.GetTotalNumHs() == 1 and a.GetNeighbors()[0].GetAtomicNum() == 6
                   and a.GetNeighbors()[0].GetTotalNumHs() == 2]
        edit = Chem.RWMol(Chem.CombineMols(donor, acceptor))
        edit.RemoveBond(phosphorus, bridge.GetIdx())
        edit.GetAtomWithIdx(bridge.GetIdx()).SetFormalCharge(-1)
        edit.AddBond(phosphorus, donor.GetNumAtoms() + oxygen, Chem.BondType.SINGLE)
    Chem.SanitizeMol(edit)
    return sorted([canonical(m) for m in Chem.GetMolFrags(edit, asMols=True)] + ['[H+]'])


class AminoPhospholipidSynthesisTests(unittest.TestCase):
    def test_every_synthesis_reconstructs_exact_products(self):
        r = json.loads(Path('data/reports/phase1-amino-phospholipid-synthesis.json').read_bytes())
        for path, sha in r['source_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), sha)
        compounds = {c['id']: c for c in r['compounds']}
        sources = {s['rule_id']: s for s in r['source_records']}
        self.assertEqual(r['summary']['balanced_equations'], 324)
        self.assertEqual(r['summary']['new_CO2_route_claims'], 0)
        count = 0
        for reaction in r['reactions']:
            self.assertTrue(balanced([reaction['left'], reaction['right']], compounds))
            self.assertFalse(reaction['enzyme_evidence_ids'])
            if reaction['hypothesis_type'] == 'amino-phospholipid-speciation':
                self.assertEqual(reaction['direction_status'], 'explicit-reversible-speciation-assumption-not-source-reaction')
                continue
            count += 1
            sides = [[compounds[p['compound_id']]['smiles'] for p in reaction[side]
                      for _ in range(p['coefficient'])] for side in ('left', 'right')]
            self.assertEqual(forward_transfer(sides[0], 'methylation' in reaction['hypothesis_type']), sorted(sides[1]))
            self.assertEqual(reaction['direction_status'], 'forward-only-generic-source-hypothesis')
            target = next(Chem.MolFromSmiles(s) for s in sides[1]
                          if s != '[H+]' and not any(a.IsInRing() for a in Chem.MolFromSmiles(s).GetAtoms()))
            reversed_mol = Chem.RenumberAtoms(target, list(range(target.GetNumAtoms()))[::-1])
            self.assertEqual(instantiate(target, sources[reaction['source_reaction_id']], reaction['hypothesis_type']),
                             instantiate(reversed_mol, sources[reaction['source_reaction_id']], reaction['hypothesis_type']))
        self.assertEqual(count, 162)
