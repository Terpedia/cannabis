import json
import hashlib
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_alkane_precursors import build,REFERENCE
from cannabis_carbon.phase1_marts_completions import balanced


def test_exact_acyl_carrier_and_cofactor_preservation():
    read=lambda n:json.loads(Path('data/reports/phase1-'+n+'.json').read_bytes())
    current,alkanes,network,report=map(read,('ketone-stereo-net','alkane-hypotheses','full-balanced-network','alkane-precursors'))
    ref=next(r for r in network['reactions'] if r['id']==REFERENCE)
    original={c['id']:c for c in current['compounds']}
    refsm=next(original[p['compound_id']]['smiles'] for p in ref['left'] if '(=O)SCC' in original[p['compound_id']]['smiles'])
    ref_chirality=[cip for _,cip in Chem.FindMolChiralCenters(Chem.MolFromSmiles(refsm),includeUnassigned=True)]
    assert {k:v for k,v in report.items() if k!='source_sha256'}==build(current,alkanes,ref)
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    cs={c['id']:c for c in report['compounds']}; rs={r['id']:r for r in report['reactions']}
    assert len(rs)==8 and len(report['pairs'])==4
    assert report['summary']['new_acyl_CoA_structures']==1
    for pair in report['pairs']:
        r=rs[pair['aldehyde_reduction_hypothesis_id']]
        assert balanced([r['left'],r['right']],cs)
        assert len(r['left'])==len(r['right'])==3
        assert all(p['coefficient']==1 for s in ('left','right') for p in r[s])
        mol=Chem.MolFromSmiles(cs[pair['acyl_CoA_compound_id']]['smiles'])
        assert [cip for _,cip in Chem.FindMolChiralCenters(mol,includeUnassigned=True)]==ref_chirality
        assert pair['acyl_CoA_present_in_parent_model']==(pair['aldehyde_carbon_count']!=21)
        assert r['enzyme_evidence_ids']==[]
