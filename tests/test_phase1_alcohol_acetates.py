import hashlib
import json
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_alcohol_acetates import build, precursor
from cannabis_carbon.phase1_marts_completions import balanced


def test_scope_exclusions_and_exact_stereo_retention():
    assert precursor('CCO')[1] == 'not-an-acetate-ester'
    assert precursor('CC(=O)Oc1ccccc1')[1] == 'phenolic-or-unsaturated-oxygen-substrate-separate-review'
    assert precursor('CC(=O)OC(C)(C)C')[1] == 'tertiary-alcohol-substrate-separate-review'
    assert precursor('CC(=O)OCCOC(C)=O')[1] == 'multiple-esters-separate-site-specific-review'
    assert precursor('[13CH3]C(=O)OCC')[1] == 'charged-or-isotope-specific-target-separate-review'
    a = precursor('CC(=O)O[C@H](C)CC')[0]
    b = precursor('CC(=O)O[C@@H](C)CC')[0]
    assert a != b


def test_every_proposal_roundtrips_to_exact_target_and_preserves_cofactors():
    read = lambda n: json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
    current, network, report = map(read, ('c17-elongation-net', 'full-balanced-network', 'alcohol-acetates'))
    assert build(current, network) == {k: v for k, v in report.items() if k != 'source_sha256'}
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    assert report['summary'] == {'acetate_target_records_reviewed': 14, 'balanced_proposed_equations': 6,
                                 'new_compound_structures': 0, 'coverage_gain_claimed': 0}
    compounds = {c['id']: c for c in report['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    for t in report['targets']:
        if 'alcohol_smiles' not in t:
            assert 'hypothesis_id' not in t
            continue
        alcohol = Chem.MolFromSmiles(t['alcohol_smiles'])
        products = set()
        for atom in alcohol.GetAtoms():
            if atom.GetAtomicNum() != 8 or atom.GetTotalNumHs() != 1:
                continue
            edit = Chem.RWMol(Chem.CombineMols(alcohol, Chem.MolFromSmiles('CC=O')))
            edit.AddBond(atom.GetIdx(), alcohol.GetNumAtoms()+1, Chem.BondType.SINGLE)
            mol = edit.GetMol(); Chem.SanitizeMol(mol)
            products.add(Chem.MolToSmiles(mol,isomericSmiles=True))
        assert t['target_smiles'] in products
        if 'hypothesis_id' not in t:
            assert not t['alcohol_present_in_model']
            continue
        r = reactions[t['hypothesis_id']]
        assert balanced([r['left'], r['right']], compounds)
        assert len(r['left']) == len(r['right']) == 2
        assert all(p['coefficient'] == 1 for s in ('left','right') for p in r[s])
        assert r['enzyme_evidence_ids'] == [] and r['direction_status'] == 'proposed-forward-only'
        unchanged = {p['compound_id'] for s in ('left','right') for p in r[s]} - {t['compound_id'],t['alcohol_compound_id']}
        assert len(unchanged) == 2
        assert any(compounds[c]['smiles'].startswith('CC(=O)SCC') for c in unchanged)
        assert any(compounds[c]['smiles'].endswith('NCCS') for c in unchanged)
