import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_cpr_full_proteome_search_retains_carrier_and_partnership_gaps(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-cpr-search')
    report = json.loads(Path('data/reports/phase1-cpr-search.json').read_text())
    discovery = json.loads(Path('data/reports/phase1-cpr-references.json').read_text())
    review = json.loads(Path('data/curation/cpr-fnsii-carrier-interface-review.json').read_text())
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    row = report['rows'][0]
    assert row['reaction_id'] == 'RHEA:24040'
    assert row['carrier_interface_review'] == review
    assert set(row['screened_cannabis_proteins']) == set(review['candidate_accessions'])
    assert row['compatible_fnsii_partners'] == []
    assert row['partner_status'] == 'unverified'
    assert report['model_eligible'] is row['model_eligible'] is discovery['model_eligible'] is False
    assert all(r['model_eligible'] is False and not r['exact_reaction_annotation_match'] for r in row['reference_matches'])
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['retrieved_references'] == 11
    assert report['summary']['raw_alignments'] == 113
    assert report['summary']['passing_alignments'] == 33
    assert report['summary']['distinct_cannabis_candidates'] == 3
    raw = Path(report['hits_path']).read_text().splitlines()
    assert len({line.split('\t')[0] for line in raw}) == 11
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay({'reactions': [{'id': row['reaction_id']}], 'hypotheses': []}, report)
