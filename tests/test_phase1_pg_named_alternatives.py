import hashlib
import json
from pathlib import Path
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_pg_named_alternatives import derive,parse_name
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_acyl_input_identity_review import double_stereo_key


def test_all_named_alternatives_keep_originals_and_replay_exact_geometry():
    r=json.loads(Path('data/reports/phase1-pg-named-alternatives.json').read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    assert len(r['targets'])==23
    assert r['summary']['changed_historical_identities']==r['summary']['new_CO2_route_claims']==0
    for row in r['targets']:
        original=Chem.MolFromSmiles(row['original_smiles'])
        alternative=Chem.MolFromSmiles(row['alternative_smiles'])
        assert row['original_compound_id']!=row['alternative_compound_id']
        assert double_stereo_key(original)==double_stereo_key(alternative)
        assert canonical(derive(original,row['source_name'])[0])==row['alternative_smiles']
        reversed_mol=Chem.RenumberAtoms(original,list(range(original.GetNumAtoms()))[::-1])
        assert canonical(derive(reversed_mol,row['source_name'])[0])==row['alternative_smiles']
        # Independent distance-based positional audit on reparsed output.
        named=parse_name(row['source_name'])
        for carbonyl,oxygen,ester,glycerol in alternative.GetSubstructMatches(Chem.MolFromSmarts('[C](=[O])-[O]-[C]')):
            sn={2:1,1:2}[alternative.GetAtomWithIdx(glycerol).GetTotalNumHs()]
            cut=Chem.RWMol(alternative);cut.RemoveBond(carbonyl,ester)
            fragment=next(set(f) for f in Chem.GetMolFrags(cut) if carbonyl in f)
            carbons={i for i in fragment if alternative.GetAtomWithIdx(i).GetAtomicNum()==6}
            assert len(carbons)==named[sn-1][0]
            found={}
            for bond in alternative.GetBonds():
                a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                if bond.GetBondType()==Chem.BondType.DOUBLE and a in carbons and b in carbons:
                    position=min(len(Chem.GetShortestPath(alternative,carbonyl,x)) if x!=carbonyl else 1 for x in (a,b))
                    assert bond.GetStereo() in (Chem.BondStereo.STEREOZ,Chem.BondStereo.STEREOE)
                    found[position]='Z' if bond.GetStereo()==Chem.BondStereo.STEREOZ else 'E'
            assert found==named[sn-1][1]


def test_incomplete_name_or_conflicting_structure_is_rejected():
    row=json.loads(Path('data/reports/phase1-pg-named-alternatives.json').read_bytes())['targets'][0]
    original=Chem.MolFromSmiles(row['original_smiles'])
    with pytest.raises(ValueError):parse_name('PG(16:0/18:1)')
    with pytest.raises(ValueError):derive(original,row['source_name'].replace('16:0','17:0'))
    with pytest.raises(ValueError):derive(Chem.MolFromSmiles(row['alternative_smiles']),row['source_name'].replace('Z','E'))
