import copy
import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_producer_screen_audit import build

ROOT = Path(__file__).resolve().parents[1]


def test_all_exact_producer_search_rows_are_retained_without_promotion():
    report = json.loads((ROOT / 'data/reports/phase1-producer-screen-audit.json').read_text())
    audit = json.loads((ROOT / 'data/reports/phase1-no-producer-audit.json').read_text())
    searches = {p: json.loads((ROOT / p).read_text()) for p in report['search_report_paths']}
    before = copy.deepcopy(searches)
    core = json.loads((ROOT / 'docs/data/networkdb.json').read_text())
    assert build(audit, searches, core) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert searches == before
    assert len(report['rows']) == 208 and len(report['targets']) == 120
    assert all(r['prior_searches'] and not r['admitted_to_candidate_model'] for r in report['rows'])
    assert sum(len(r['prior_searches']) for r in report['rows']) == 318
    assert sum(r['disposition'].startswith('prior-candidate-lead') for r in report['rows']) == 2
    core_by_id = {r['id']: r for r in core['reactions']}
    for row in report['rows']:
        expected = sorted({s['source_reaction_id'] for s in row['reaction']['sources']} & core_by_id.keys())
        assert row['core_source_reactions'] == [core_by_id[s] for s in expected]
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_no_reference_is_not_a_negative_alignment_and_candidates_keep_context():
    base = {'targets': [{'cannabisdb_id': 'CDB1', 'exact_catalog_net_producer_reaction_ids': ['r'],
                         'source_identity_status': 'source-structure-disagreement'}], 'reactions': [{'id': 'r'}]}
    for status, expected in [('no-reference-sequence', 'reference-sequence-gap'), ('no-hits', 'prior-no-hits'), ('weak-hits-only', 'prior-weak-hits')]:
        result = build(base, {'search': {'rows': [{'reaction_id': 'r', 'search_status': status}]}})
        assert result['rows'][0]['disposition'].startswith(expected)
        assert result['rows'][0]['identity_conflict_target_ids'] == ['CDB1']
    assert build(base, {})['rows'][0]['disposition'] == 'no-prior-search-record'
