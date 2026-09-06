import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from cannabis_carbon.phase1_catalog import stable_id


def test_mixture_evidence_remains_separate_from_exact_isomer_channels():
    report = json.loads(Path('data/reports/phase1-citral-reductase-review.json').read_bytes())
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    current = json.loads(Path('data/reports/phase1-alcohol-acetates-net.json').read_bytes())
    compounds = {c['id']: c for c in current['compounds']}
    evidence = report['evidence']
    assert evidence == json.loads(Path('data/curation/citral-reductase-evidence.json').read_bytes())
    assert not evidence['reaction_added'] and not evidence['single_isomer_reaction_supported_by_this_review']
    assert not evidence['enantiopure_S_production_supported'] and not evidence['Cannabis_activity_supported']
    assert [r['S_percent'] for r in evidence['NADH_product_distributions']] == [92, 97, 56]
    assert all(r['S_percent'] + r['R_percent_by_complement'] == 100 for r in evidence['NADH_product_distributions'])
    assert len(report['identities']) == 3 and len({r['compound_id'] for r in report['identities']}) == 3
    formulas = []
    for row in report['identities']:
        mol = Chem.MolFromSmiles(row['smiles'])
        assert stable_id('structure',Chem.MolToSmiles(mol,isomericSmiles=True)) == row['compound_id']
        assert row['present_in_model'] and compounds[row['compound_id']]['smiles'] == row['smiles']
        assert row['tetrahedral_CIP_labels'] == [label for _,label in Chem.FindMolChiralCenters(mol,includeUnassigned=True)]
        assert row['double_bond_stereo'] == [str(b.GetStereo()) for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        formulas.append(rdMolDescriptors.CalcMolFormula(mol))
    assert formulas == ['C10H18O', 'C10H16O', 'C10H16O']
    assert report['identities'][0]['tetrahedral_CIP_labels'] == ['S']
    assert 'STEREOE' in report['identities'][1]['double_bond_stereo']
    assert 'STEREOZ' in report['identities'][2]['double_bond_stereo']
    assert report['summary'] == {'resolved_model_identity_candidates': 3, 'new_reactions': 0,
        'coverage_change': 0, 'new_Cannabis_enzyme_assignments': 0}
