import hashlib
import json
from pathlib import Path
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_alternative_substrate_screen_does_not_close_delta6_gap():
    root = Path(__file__).resolve().parents[1]
    verify_search('phase1-sphingolipid-alternative-search')
    report = json.loads((root / 'data/reports/phase1-sphingolipid-alternative-search.json').read_text())
    discovery = json.loads((root / report['source_discovery']).read_text())
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['passing_alignments'] == 4
    assert {p['accession'] for p in report['cannabis_candidates']} == {'A0A7J6DP00', 'A0A7J6F905'}
    assert all(h['identity_percent'] >= 64 and h['reference_coverage_percent'] >= 99 for h in report['passing_alignments'])
    assert report['rows'][0]['reaction_id'] == 'alternative-hypothesis:sphingolipid-delta8'
    assert report['rows'][0]['target_ids'] == []
    assert all(m['model_eligible'] is False for m in discovery['rows'][0]['reference_matches'])
    for key in ('proteome', 'reference', 'hits'):
        assert (root / report[key + '_path']).is_file()
