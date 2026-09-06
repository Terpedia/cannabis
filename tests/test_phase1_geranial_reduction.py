import hashlib
import json
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_geranial_reduction import build, NADPH, NADP, PROTON
from cannabis_carbon.phase1_marts_completions import balanced


def test_exact_geranial_S_product_cofactors_and_evidence():
    read = lambda name: json.loads(Path('data/reports/phase1-' + name + '.json').read_bytes())
    report = read('geranial-reduction')
    rebuilt = build(read('alcohol-acetates-net'), read('full-balanced-network'))
    assert {k: v for k, v in report.items() if k != 'source_sha256'} == rebuilt
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    r, = report['reactions']; c = {c['id']: c for c in report['compounds']}
    assert balanced([r['left'], r['right']], c)
    assert all(p['coefficient'] == 1 for side in ('left', 'right') for p in r[side])
    assert {p['compound_id'] for p in r['left'][1:]} == {NADPH, PROTON}
    assert r['right'][1]['compound_id'] == NADP
    substrate = Chem.MolFromSmiles(c[r['left'][0]['compound_id']]['smiles'])
    product = Chem.MolFromSmiles(c[r['right'][0]['compound_id']]['smiles'])
    assert [b.GetStereo() for b in substrate.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE] == [Chem.BondStereo.STEREOE]
    assert [a.GetProp('_CIPCode') for a in product.GetAtoms() if a.HasProp('_CIPCode')] == ['S']
    assert r['enzyme_evidence_ids'] == []
    assert r['direction_status'] == 'proposed-forward-only'
    e, = report['biochemical_evidence']
    assert e['reducing_cofactor'] == 'NADPH'
    assert not e['Cannabis_activity_established'] and not e['exclusive_S_product_claimed']
    assert report['summary']['coverage_gain_claimed'] == 0
