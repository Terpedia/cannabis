import hashlib
import json
from pathlib import Path

from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_full_proteome_generic_reference_search_preserves_negative_result():
    root = Path(__file__).resolve().parents[1]
    verify_search('phase1-biopterin-lead-search')
    report = json.loads((root / 'data/reports/phase1-biopterin-lead-search.json').read_text())
    discovery = json.loads((root / report['source_discovery']).read_text())
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    for name in ('proteome', 'reference', 'hits'):
        assert (root / report[name + '_path']).is_file()
    raw = json.loads((root / report['retrievals'][0]['snapshot']).read_text())
    assert report['reference_sequences'][0]['sequence'] == raw['sequence']['value']
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['retrieved_references'] == 1
    assert report['summary']['raw_alignments'] == 0
    assert report['summary']['distinct_cannabis_candidates'] == 0
    row = report['rows'][0]
    assert row['search_status'] == 'no-hits'
    assert row['reference_matches'][0]['model_eligible'] is False
    assert row['reference_matches'][0]['exact_reaction_annotation_match'] is False
    assert 'not-eligible-for-exact-reaction-model' in row['validation_blockers']
