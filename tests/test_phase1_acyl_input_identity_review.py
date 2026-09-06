import hashlib
import json
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_acyl_input_identity_review import double_stereo_key, build


def test_comparison_removes_only_double_stereo_and_does_not_mutate_input():
    mol = Chem.MolFromSmiles('C/C=C/[C@H](O)[13CH3]')
    original = Chem.MolToSmiles(mol)
    key = double_stereo_key(mol)
    assert Chem.MolToSmiles(mol) == original
    assert key == double_stereo_key(Chem.MolFromSmiles('C/C=C\\[C@H](O)[13CH3]'))
    assert key != double_stereo_key(Chem.MolFromSmiles('C/C=C/[C@@H](O)[13CH3]'))
    assert key != double_stereo_key(Chem.MolFromSmiles('C/C=C/[C@H](O)C'))
    assert key != double_stereo_key(Chem.MolFromSmiles('C/C=C/[C@H]([O-])[13CH3]'))


def test_full_inventory_review_is_reproducible_without_identity_promotion():
    r = json.loads(Path('data/reports/phase1-acyl-input-identity-review.json').read_bytes())
    audit = json.loads(Path('data/reports/phase1-glycerophospholipid-input-audit.json').read_bytes())
    current = json.loads(Path('data/reports/phase1-glycerophospholipid-net.json').read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    assert build(audit,current) == {k:v for k,v in r.items() if k!='source_sha256'}
    assert r['summary']['unresolved_input_structures'] == 5
    assert r['summary']['affected_target_records'] == 23
    assert r['summary']['inputs_with_unassigned_alkene_geometry'] == 5
    assert r['summary']['new_CO2_route_claims'] == r['summary']['new_reactions'] == 0
    for row in r['inputs']:
        assert row['unassigned_alkene_bonds'] > 0
        assert row['candidate_structures']
        for candidate in row['candidate_structures']:
            assert candidate['compound_id'] != row['compound_id']
            assert candidate['smiles'] != row['smiles']
            assert double_stereo_key(Chem.MolFromSmiles(candidate['smiles'])) == row['comparison_key']
            assert candidate['identity_relationship'].endswith('not-exact-identity')
