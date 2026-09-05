import hashlib
import json
from pathlib import Path

from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_cannabinoid_revalidation import prepare
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search

ROOT = Path(__file__).resolve().parents[1]


def test_full_proteome_revalidation_and_exact_cbdas_sequence():
    verify_search('phase1-cannabinoid-revalidation-search')
    report = json.loads((ROOT / 'data/reports/phase1-cannabinoid-revalidation-search.json').read_text())
    assert report['summary']['proteome_sequences'] == 30304
    assert report['summary']['raw_alignments'] == 79
    assert report['summary']['passing_alignments'] == 67
    references = {r['accession']: r['sequence'] for r in report['reference_sequences']}
    proteins = _fasta(ROOT / report['proteome_path'])
    assert proteins['A0A7J6G9C8'] == references['A6P6V9']
    assert len(proteins['A0A7J6G9C8']) == 544
    assert len(proteins['A0A7J6DJS8']) == 368 < len(references['Q8GTB6']) == 545


def test_legacy_reference_discovery_replays_without_repairing_old_claim():
    path = ROOT / 'data/reports/phase1-cannabinoid-revalidation-references.json'
    report = json.loads(path.read_text())
    audit = json.loads((ROOT / 'data/reports/phase1-producer-screen-audit.json').read_text())
    references = {a: json.loads((ROOT / f'data/raw/phase1-cannabinoid-revalidation/{a}.json').read_text()) for a in ('A6P6V9', 'Q8GTB6')}
    sequences = _fasta(ROOT / 'data/raw/UP000583929.fasta')
    assert prepare(audit, sequences, references) == {k: v for k, v in report.items() if k not in ('lookups', 'source_sha256')}
    assert [c['legacy_query_lengths_match_pinned_sequence'] for c in report['legacy_checks']] == [False, True]
    for source, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / source).read_bytes()).hexdigest() == sha
    for lookup in report['lookups']:
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']
