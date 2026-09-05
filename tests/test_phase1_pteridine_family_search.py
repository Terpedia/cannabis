import csv
import hashlib
import json
from pathlib import Path

from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_reviewed_generic_family_screen_retains_all_references_and_no_hits():
    root = Path(__file__).resolve().parents[1]
    verify_search('phase1-pteridine-family-search')
    report = json.loads((root / 'data/reports/phase1-pteridine-family-search.json').read_text())
    discovery = json.loads((root / report['source_discovery']).read_text())
    for path, digest in discovery['source_sha256'].items():
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == digest
    lookup = discovery['lookups'][0]
    assert 'rhea:17865' in lookup['query'] and 'rhea:17869' in lookup['query']
    assert 'reviewed:true' in lookup['query'] and 'fragment:false' in lookup['query']
    with (root / lookup['snapshot']).open() as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    matches = discovery['rows'][0]['reference_matches']
    assert {m['accession'] for m in matches} == {r['Entry'] for r in records}
    assert len(matches) == 7
    assert all(not m['model_eligible'] and not m['exact_reaction_annotation_match'] for m in matches)
    directions = (root / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text().splitlines()
    families = {line.split('\t')[0]: line.split('\t') for line in directions}
    allowed = {f'RHEA:{v}' for master in ('17865', '17869') for v in families[master]}
    assert all(set(m['matched_generic_rhea_ids']) <= allowed for m in matches)
    for name in ('proteome', 'reference', 'hits'):
        assert (root / report[name + '_path']).is_file()
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['retrieved_references'] == 7
    assert report['summary']['raw_alignments'] == 0
    assert report['rows'][0]['search_status'] == 'no-hits'
