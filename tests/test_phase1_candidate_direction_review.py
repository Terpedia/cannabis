import copy
import hashlib
import json
from pathlib import Path

import pytest

from cannabis_carbon.phase1_candidate_direction_review import build
from cannabis_carbon.phase1_gap_annotations import assemble
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def inputs():
    report = json.loads((ROOT/'data/reports/phase1-expanded-candidate-net.json').read_text())
    review = json.loads((ROOT/'data/reports/phase1-candidate-direction-review.json').read_text())
    families = direction_families((ROOT/'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    return report, review, families


def test_source_orientation_replay_and_unchanged_certificates(monkeypatch):
    monkeypatch.chdir(ROOT)
    report, review, families = inputs()
    annotations = assemble([], review['lookups'])['source_annotations']
    before = copy.deepcopy(report)
    assert build(report, families, annotations) == {k:v for k,v in review.items() if k not in ('lookups','source_sha256')}
    assert report == before
    for path, sha in review['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha
    assert len(review['reviews']) == 5
    assert len(review['targets']) == 7
    assert all(t['direction_review_ids'] for t in review['targets'])
    for r in review['reviews']:
        lr = next(s for s in r['source_joins'] if s['source_reaction_id'] == r['published_family']['RHEA_ID_LR'])
        for use in r['uses']:
            source_input = 'left' if use['direction_mode']=='hypothetical-left-to-right' else 'right'
            assert source_input != lr['source_left_corresponds_to']
            assert use['required_inputs'] == r[source_input]
            assert use['source_orientation'] == 'opposite-to-source-written'
        assert len(r['discriminating_tests']) == 4
    bundle = json.loads((ROOT/'docs/data/expanded-net-view/bundle.json').read_text())
    assert bundle['certificates'] == report['certificates']
    assert {r['direction_review']['id'] for r in bundle['reactions'] if 'direction_review' in r} == {r['id'] for r in review['reviews']}


def test_no_automatic_reverse_classification_and_conflict_rejected(monkeypatch):
    monkeypatch.chdir(ROOT)
    report, review, families = inputs()
    annotations = assemble([], review['lookups'])['source_annotations']
    for cert in report['certificates']:
        for step in cert['steps']:
            step['direction_mode'] = 'hypothetical-right-to-left' if step['direction_mode']=='hypothetical-left-to-right' else 'hypothetical-left-to-right'
    changed = build(report, families, annotations)
    assert not any(t['direction_review_ids'] for t in changed['targets'])
    assert all('follow source-written' in r['warning'] for r in changed['reviews'])
    reaction = next(r for r in report['reactions'] if r['new_catalog_candidate'])
    reaction['sources'][1]['source_left_corresponds_to'] = reaction['sources'][0]['source_left_corresponds_to']
    with pytest.raises(ValueError, match='Conflicting'):
        build(report, families, annotations)
