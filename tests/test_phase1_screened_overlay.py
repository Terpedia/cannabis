import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_screened_overlay import apply_overlay, build_overlay


def test_overlay_only_changes_candidate_evidence_and_missing_candidate_blocker():
    parent = {'compounds': [], 'reactions': [{'id': 'r', 'left': ['all-inputs'], 'right': ['all-products']}],
        'targets': [{'cannabisdb_id': 'A', 'carbon_count': 1}], 'summary': {}, 'enzyme_evidence': [],
        'hypotheses': [{'id': 'h', 'reaction_id': 'r', 'cannabisdb_id': 'A', 'has_candidate_enzyme_evidence': False,
            'evidence_ids': [], 'status': 'blocked', 'direction_mode': 'hypothetical-right-to-left',
            'required_inputs': ['water', 'substrate'], 'net_target_coefficient': 2,
            'blockers': ['no-candidate-enzyme-evidence-attached', 'Cannabis-enzyme-activity-unconfirmed', 'requires-target-bootstrap']}]}
    before = copy.deepcopy(parent)
    overlay = {'enzyme_evidence': [{'reaction_id': 'r', 'id': 'ev', 'claim_boundary': 'Homology only'}]}
    result = apply_overlay(parent, overlay)
    assert parent == before
    assert result['reactions'] == before['reactions']
    h = result['hypotheses'][0]
    assert h['status'] == 'blocked'
    assert h['required_inputs'] == ['water', 'substrate']
    assert h['net_target_coefficient'] == 2
    assert h['direction_mode'] == 'hypothetical-right-to-left'
    assert h['blockers'] == ['Cannabis-enzyme-activity-unconfirmed', 'requires-target-bootstrap']
    assert result['summary']['carbon_bearing_targets_with_candidate_enzyme_evidence'] == 1
    with pytest.raises(ValueError, match='collision'):
        apply_overlay(result, overlay)


def test_published_integration_preserves_every_chemical_and_pathway_claim():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-screened-enzyme-overlay.json'
    overlay = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-screened-enzyme-overlay.json').read_bytes()
    for source, digest in overlay['source_sha256'].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest
    parent = json.loads((root / 'data/reports/phase1-target-hypotheses.json').read_text())
    search = json.loads((root / 'data/reports/phase1-new-protein-search.json').read_text())
    rebuilt = build_overlay(parent, search)
    assert rebuilt['enzyme_evidence'] == overlay['enzyme_evidence']
    report = apply_overlay(parent, overlay)
    assert report['summary'] == overlay['integrated_summary']
    assert report['reactions'] == parent['reactions']
    assert report['compounds'] == parent['compounds']
    assert report['targets'] == parent['targets']
    supported_reactions = {r['reaction_id'] for r in search['rows'] if r['screened_cannabis_proteins']}
    assert {r['reaction_id'] for r in overlay['enzyme_evidence']} == supported_reactions
    for old, new in zip(parent['hypotheses'], report['hypotheses']):
        for key in ['id', 'reaction_id', 'cannabisdb_id', 'compound_id', 'required_inputs', 'outputs', 'net_target_coefficient', 'status', 'direction_mode']:
            assert old[key] == new[key]
        assert new['has_candidate_enzyme_evidence'] == (old['has_candidate_enzyme_evidence'] or old['reaction_id'] in supported_reactions)
        assert 'Cannabis-enzyme-activity-unconfirmed' in new['blockers']
        assert set(old['blockers']) - set(new['blockers']) <= {'no-candidate-enzyme-evidence-attached'}
