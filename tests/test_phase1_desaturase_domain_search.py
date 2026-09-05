import hashlib
import json
from pathlib import Path
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_source_defined_domain_search_preserves_boundaries_and_evidence():
    root = Path(__file__).resolve().parents[1]
    verify_search('phase1-desaturase-domain-search')
    report = json.loads((root / 'data/reports/phase1-desaturase-domain-search.json').read_text())
    discovery = json.loads((root / report['source_discovery']).read_text())
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    protein = json.loads((root / 'data/raw/desaturase-catalytic-domain/O95864.json').read_text())
    domain = json.loads((root / 'data/raw/desaturase-catalytic-domain/O95864-PF00487.json').read_text())
    fragment = domain['proteins'][0]['entry_protein_locations'][0]['fragments'][0]
    assert (fragment['start'], fragment['end']) == (157, 416)
    assert report['reference_sequences'][0]['sequence'] == protein['sequence']['value'][156:416]
    assert len(report['reference_sequences'][0]['sequence']) == 260
    match = discovery['rows'][0]['reference_matches'][0]
    assert match['model_eligible'] is False
    assert match['source_activity']['reaction']['evidences'] == [
        {'evidenceCode': 'ECO:0000250', 'source': 'UniProtKB', 'id': 'Q9Z0R9'}]
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['raw_alignments'] == 2
    assert report['summary']['passing_alignments'] == 0
    for key in ('proteome', 'reference', 'hits'):
        assert (root / report[key + '_path']).is_file()
