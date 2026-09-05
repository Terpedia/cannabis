import copy
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_missing_reference_review import prepare
from cannabis_carbon.phase1_new_references import attach


def test_ureide_followups_preserve_cofactors_scope_and_negative_results(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    discovery = json.loads(Path('data/reports/phase1-ureide-gap-references.json').read_text())
    search = json.loads(Path('data/reports/phase1-ureide-gap-search.json').read_text())
    baseline = prepare(discovery, search, {'no-reference-sequence', 'weak-hits-only', 'no-hits'})
    assert len(baseline) == 3
    for scope in ('plant', 'nonplant'):
        report = json.loads(Path(f'data/reports/phase1-ureide-{scope}-reference-review.json').read_text())
        lookup = report['lookups'][0]
        assert lookup['requested_master_ids'] == ['RHEA:15329', 'RHEA:15333', 'RHEA:33867']
        assert 'reviewed:false' in lookup['query'] and 'fragment:false' in lookup['query']
        assert ('NOT taxonomy_id:33090' in lookup['query']) == (scope == 'nonplant')
        assert 'taxonomy_id:33090' in lookup['query']
        assert attach(copy.deepcopy(baseline), report['lookups']) == {}
        assert report['reference_proteins'] == []
        assert {r['reaction_id'] for r in report['rows']} == {r['reaction_id'] for r in baseline}
        for row in report['rows']:
            before = next(r for r in baseline if r['reaction_id'] == row['reaction_id'])
            for key in ('left', 'right', 'rhea_families', 'selected_uses', 'target_ids', 'prior_reviewed_search'):
                assert row[key] == before[key]
            assert row['reference_matches'] == []
        for p, sha in report['source_sha256'].items():
            assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
        assert hashlib.sha256(Path(lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']


def test_cofactor_review_keeps_exact_equations_separate(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    review = json.loads(Path('data/curation/ureidoglycolate-cofactor-review.json').read_text())
    report = json.loads(Path('data/reports/phase1-ureide-gap-references.json').read_text())
    rows = {r['reaction_id']: r for r in report['rows']}
    nad, nadp = [rows[review[k + '_reaction_id']] for k in ('nad', 'nadp')]
    for key, row in (('nad', nad), ('nadp', nadp)):
        assert {f['RHEA_ID_MASTER'] for f in row['rhea_families'].values()} == {review[key + '_rhea_master']}
    assert nad['left'] != nadp['left'] and nad['right'] != nadp['right']
    assert {m['accession'] for m in nad['reference_matches']} == {'P58408', 'P77555'}
    assert nadp['reference_matches'] == []
    assert review['candidate_model_changed'] is False
    assert len(review['claims']) == 4
    assert all(c['source_url'].startswith('https://') and c['access_scope'] for c in review['claims'])
