import hashlib
import json
from pathlib import Path

from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_unreviewed_nonplant_full_proteome_screen_and_short_hit():
    root = Path(__file__).resolve().parents[1]
    verify_search('phase1-weak-nonplant-search')
    report = json.loads((root / 'data/reports/phase1-weak-nonplant-search.json').read_text())
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['requested_references'] == report['summary']['retrieved_references'] == 6629
    assert report['summary']['raw_alignments'] == 2444
    assert report['summary']['passing_alignments'] == 1
    for retrieval in report['retrievals']:
        assert retrieval['status'] == 'retrieved'
        assert not retrieval['missing_accessions']
        assert hashlib.sha256((root / retrieval['snapshot']).read_bytes()).hexdigest() == retrieval['sha256']
    for name in ('proteome', 'reference', 'hits'):
        assert (root / report[name + '_path']).is_file()
    hit = report['passing_alignments'][0]
    assert hit['cannabis_accession'] == 'A0A7J6FB06'
    assert hit['alignment_length'] == 84
    assert hit['reference_accession'] == 'A0ACB9UWW7'
    row = next(r for r in report['rows'] if r['screened_cannabis_proteins'])
    assert row['evidence_class'] == 'unreviewed-nonplant-reference-homology-candidate'
    assert 'unreviewed-reference-activity-unverified' in row['validation_blockers']
    assert 'catalytic-residues-and-domains-not-reviewed' in row['validation_blockers']
