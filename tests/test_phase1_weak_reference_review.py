import copy
import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_missing_reference_review import prepare
from cannabis_carbon.phase1_new_references import attach


def test_weak_search_followups_preserve_exact_scope_and_review_status(monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.chdir(root)
    discovery = json.loads(Path('data/reports/phase1-remaining-gap-references.json').read_text())
    prior = json.loads(Path('data/reports/phase1-remaining-gap-search.json').read_text())
    selected = prepare(discovery, prior, ('weak-hits-only', 'no-hits'))
    assert len(selected) == 4
    for scope in ('plant', 'nonplant'):
        report = json.loads(Path(f'data/reports/phase1-weak-{scope}-reference-review.json').read_text())
        for path, digest in report['source_sha256'].items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
        lookup = report['lookups'][0]
        assert Path(lookup['snapshot']).is_file()
        assert 'reviewed:false' in lookup['query']
        assert ('NOT taxonomy_id:33090' in lookup['query']) == (scope == 'nonplant')
        rows = copy.deepcopy(selected)
        proteins = attach(rows, report['lookups'])
        assert {r['reaction_id'] for r in rows} == {r['reaction_id'] for r in report['rows']}
        for expected, actual in zip(rows, report['rows']):
            assert expected['prior_reviewed_search'] == actual['prior_reviewed_search']
            assert expected['rhea_families'] == actual['rhea_families']
            assert expected['target_ids'] == actual['target_ids']
            assert [{**m, 'review_status': 'unreviewed'} for m in expected['reference_matches']] == actual['reference_matches']
        assert len(proteins) == report['summary'][f'unreviewed_{scope}_reference_records']
