import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_fulltext_review_preserves_source_exact_identity_and_no_model_promotion():
    review = json.loads((ROOT / 'data/curation/glucuronolactone-fulltext-review.json').read_text())
    assert hashlib.sha256((ROOT / review['snapshot']).read_bytes()).hexdigest() == review['snapshot_sha256']
    catalog = json.loads((ROOT / 'data/reports/phase1-full-balanced-network.json').read_text())
    reaction = next(r for r in catalog['reactions'] if r['id'] == review['reaction_id'])
    identity = review['identity_review']
    compound = next(c for c in catalog['compounds'] if c['id'] == identity['catalog_product_id'])
    assert compound['smiles'] == identity['catalog_product_smiles']
    assert identity['catalog_product_id'] in {m['compound_id'] for m in reaction['right']}
    assert identity['exact_structure_identity_established'] is False
    assert review['candidate_model_changed'] is False
    baseline = json.loads((ROOT / 'data/reports/phase1-remaining-candidate-net.json').read_text())
    assert review['reaction_id'] not in baseline['candidate_reaction_evidence_ids']
    assert review['observations'][0]['starting_material'] == 'lactone'
    assert review['observations'][0]['elapsed_unit'] == 'days'
