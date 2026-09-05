import copy
import hashlib
import json
from pathlib import Path

import pytest
from cannabis_carbon.phase1_catalog_evidence import build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_supplement_exact_rebuild_preserves_every_frozen_certificate():
    parent, search, supplement = [read(n) for n in (
        'phase1-catalog-net-gaps', 'phase1-catalog-protein-search', 'phase1-catalog-evidence')]
    before = copy.deepcopy(parent)
    for path, digest in supplement['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    assert build(parent, search) == {k:v for k,v in supplement.items() if k != 'source_sha256'}
    assert parent == before
    added = {e['reaction_id']: e for e in supplement['enzyme_evidence']}
    assert len(added) == 97
    updates = {c['compound_id']: c for c in supplement['certificate_updates']}
    for cert in parent['certificates']:
        expected = [rid for rid in cert['missing_candidate_reaction_ids'] if rid not in added]
        if expected != cert['missing_candidate_reaction_ids']:
            assert updates[cert['compound_id']]['missing_candidate_reaction_ids'] == expected
            assert updates[cert['compound_id']]['baseline_missing_candidate_reaction_ids'] == cert['missing_candidate_reaction_ids']
        else:
            assert cert['compound_id'] not in updates
    assert supplement['summary']['selected_certificate_targets_with_candidates_for_all_steps'] == 102
    assert supplement['summary']['remaining_missing_candidate_equations'] == 368
    assert supplement['summary']['newly_candidate_linked_selected_certificate_target_ids'] == ['CDB004839']
    payload = (ROOT / 'data/reports/phase1-catalog-evidence.json').read_bytes()
    assert payload == (ROOT / 'docs/data/catalog-net-view/evidence.json').read_bytes()
    manifest = json.loads((ROOT / 'docs/data/catalog-net-view/index.json').read_text())
    assert manifest['evidence']['sha256'] == hashlib.sha256(payload).hexdigest()
    assert manifest['evidence']['bytes'] == len(payload)


def test_supplement_rejects_existing_evidence_and_mismatched_proteins():
    parent, search = read('phase1-catalog-net-gaps'), read('phase1-catalog-protein-search')
    row = next(r for r in search['rows'] if r['passing_alignment_ids'])
    next(r for r in parent['reactions'] if r['id'] == row['reaction_id'])['enzyme_evidence_ids'] = ['existing']
    with pytest.raises(ValueError, match='overlaps'):
        build(parent, search)
    parent = read('phase1-catalog-net-gaps')
    row['screened_cannabis_proteins'] = []
    with pytest.raises(ValueError, match='do not match'):
        build(parent, search)
