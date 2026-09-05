import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_archived_references import archive_links, validate_entry, build
from cannabis_carbon.phase1_archived_evidence import build as build_evidence, apply
from cannabis_carbon.phase1_new_protein_search import screen

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_archive_sequence_identity_length_checksum_are_required():
    report = read('phase1-archived-references')
    for record in report['retrievals']:
        assert record['status'] == 'retrieved' and record['http_status'] == 200
        assert hashlib.sha256(record['response_text'].encode()).hexdigest() == record['response_sha256']
        entry = json.loads(record['response_text'])
        assert validate_entry(record['requested_accession'], entry) == record['reference_sequence']
        assert record['cross_references'] == entry['uniParcCrossReferences']
    record = report['retrievals'][0]; entry = json.loads(record['response_text'])
    with pytest.raises(ValueError, match='identity'):
        validate_entry('UPI0000000000', entry)
    bad = copy.deepcopy(entry); bad['sequence']['length'] += 1
    with pytest.raises(ValueError, match='length'):
        validate_entry(record['requested_accession'], bad)
    bad = copy.deepcopy(entry); bad['sequence']['md5'] = '0' * 32
    with pytest.raises(ValueError, match='MD5'):
        validate_entry(record['requested_accession'], bad)


def test_exact_original_source_links_and_inactive_reciprocity():
    discovery, review, report = [read(n) for n in ['phase1-completion-protein-discovery', 'phase1-marts-gap-references', 'phase1-archived-references']]
    before = [json.dumps(x, sort_keys=True) for x in (discovery, review, report)]
    assert build(discovery, review, report['retrievals']) == {k: v for k, v in report.items() if k != 'source_sha256'}
    links = archive_links(discovery, review)
    assert len(links) == 104 and len({r['archive_accession'] for r in links}) == 47
    assert sum(r['source_uniprot_id'] == r['archive_accession'] for r in links) == 102
    inactive = [r for r in links if r['source_uniprot_id'] != r['archive_accession']]
    assert len(inactive) == 2
    assert {r['source_uniprot_id'] for r in inactive} == {'A0A4P2VJ76'}
    assert {r['archive_accession'] for r in inactive} == {'UPI00113D540E'}
    broken = copy.deepcopy(report['retrievals'])
    next(r for r in broken if r['requested_accession'] == 'UPI00113D540E')['cross_references'] = []
    gated = build(discovery, review, broken)
    assert 'UPI00113D540E' not in {r['accession'] for r in gated['reference_sequences']}
    assert any(l['resolution_status'] == 'inactive-link-not-reciprocal; not screened' for r in gated['rows'] for l in r['archive_resolution_links'])
    assert before == [json.dumps(x, sort_keys=True) for x in (discovery, review, report)]


def test_archived_overlay_replays_scope_and_preserves_prior_evidence():
    discovery, search, parent, supplement = [read(n) for n in ['phase1-archived-references', 'phase1-archived-protein-search', 'phase1-completion-protein-evidence', 'phase1-archived-evidence']]
    before = [json.dumps(x, sort_keys=True) for x in (discovery, search, parent, supplement)]
    assert build_evidence(discovery, search, parent) == {k: v for k, v in supplement.items() if k != 'source_sha256'}
    combined = apply(parent, supplement)
    assert sum(r['has_candidate_lead'] for r in combined['rows']) == 387
    assert sum(r['has_candidate_lead'] for r in parent['rows']) == 385
    assert sum(r['has_archived_candidate_lead'] for r in supplement['rows']) == 3
    assert supplement['summary']['new_equations_with_candidate_lead'] == 2
    extra = {r['id']: r for r in supplement['rows']}
    for old, new in zip(parent['rows'], combined['rows']):
        assert old['id'] == new['id'] and old['reaction_id'] == new['reaction_id']
        assert set(old['screened_cannabis_proteins']) <= set(new['screened_cannabis_proteins'])
        if old['id'] not in extra:
            assert old == new
    assert before == [json.dumps(x, sort_keys=True) for x in (discovery, search, parent, supplement)]
    for name in ['phase1-archived-references', 'phase1-archived-protein-search', 'phase1-archived-evidence']:
        path = ROOT / 'data/reports' / (name + '.json')
        assert path.read_bytes() == (ROOT / 'docs/data' / path.name).read_bytes()
        report = json.loads(path.read_text())
        for source, digest in report.get('source_sha256', {}).items():
            assert hashlib.sha256((ROOT / source).read_bytes()).hexdigest() == digest
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['raw_alignments'] == 107 and search['summary']['passing_alignments'] == 65
    assert search['summary']['distinct_cannabis_candidates'] == 44
    assert all('archived-sequence-is-not-functional-annotation' in r['validation_blockers'] for r in search['rows'])


def test_prepared_screen_rejects_out_of_scope_and_corrupted_sequences(tmp_path):
    discovery = {'rows': [{'reference_matches': [{'accession': 'UPI0000000001'}]}]}
    ref = {'accession': 'UPI0000000001', 'header': 'UPI0000000001', 'sequence': 'M', 'sequence_sha256': 'bad'}
    with pytest.raises(ValueError, match='checksum'):
        screen(discovery, tmp_path / 'source.json', tmp_path, tmp_path / 'out.json', {ref['accession']: ref}, [])
    with pytest.raises(ValueError, match='scope'):
        screen(discovery, tmp_path / 'source.json', tmp_path, tmp_path / 'out.json', {'wrong': ref}, [])
