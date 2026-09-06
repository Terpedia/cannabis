import json
import hashlib
from pathlib import Path
from rdkit import Chem


def test_reviews_join_exact_hypotheses_without_promoting_candidate_support():
    review=json.loads(Path('data/curation/ketone-substrate-evidence-review.json').read_bytes())
    hypotheses=json.loads(Path('data/reports/phase1-ketone-stereo-hypotheses.json').read_bytes())
    compounds={c['id']:c for c in hypotheses['compounds']}
    reactions={r['id']:r for r in hypotheses['reactions']}
    assert len(review['reviews'])==2
    for row in review['reviews']:
        assert row['exact_target_smiles']==compounds[row['compound_id']]['smiles']
        mol=Chem.MolFromSmiles(row['exact_target_smiles'])
        assert [cip for _,cip in Chem.FindMolChiralCenters(mol,includeUnassigned=True)]==row['target_CIP']
        reaction=reactions[row['hypothesis_id']]
        assert reaction['hypothesis_type']=='stereoselective-ketone-reduction'
        assert row['compound_id'] in {p['compound_id'] for p in reaction['right']}
        assert not reaction['enzyme_evidence_ids']
        assert row['assignment_decision'].startswith('do-not-')
        assert row['scope_limit'] and row['next_test']
        assert all(s['url'].startswith('https://') and s['evidence_level'] for s in row['sources'])


def test_export_retains_complete_review_and_pins_sources():
    review=json.loads(Path('data/curation/ketone-substrate-evidence-review.json').read_bytes())
    report=json.loads(Path('data/reports/phase1-ketone-substrate-review.json').read_bytes())
    assert {k:v for k,v in report.items() if k!='source_sha256'}==review
    assert len(report['source_sha256'])==2
    for path,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
    assert Path('docs/data/ketone-substrate-review.json').read_bytes()==Path('data/reports/phase1-ketone-substrate-review.json').read_bytes()
