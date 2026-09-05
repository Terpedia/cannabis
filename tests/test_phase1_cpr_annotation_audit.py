import hashlib
import json
from pathlib import Path
from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_cpr_annotation_audit import NAME


def test_cpr_inventory_preserves_annotation_only_scope_and_carrier_mismatch(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-cpr-annotation-audit.json').read_text())
    parent = json.loads(Path('data/reports/phase1-fnsii-search.json').read_text())
    proteome = Path(parent['proteome_path'])
    sequences = _fasta(proteome)
    selected = [h[1:] for h in proteome.read_text().splitlines() if h.startswith('>') and NAME in h]
    assert [r['source_header'] for r in report['rows']] == selected
    assert len(sequences) == report['summary']['proteome_sequences'] == 30304
    assert len(report['rows']) == 3
    assert report['related_fnsii_candidate_accessions'] == sorted(r['accession'] for r in parent['cannabis_candidates'])
    assert len(report['related_fnsii_candidate_accessions']) == 95
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    review = json.loads(Path('data/curation/cpr-fnsii-carrier-interface-review.json').read_text())
    assert {r['accession'] for r in report['rows']} == set(review['candidate_accessions'])
    for row in report['rows']:
        assert row['annotation'] == json.loads(Path(row['snapshot']).read_text())
        assert row['sequence'] == sequences[row['accession']] == row['annotation']['sequence']['value']
        assert row['compatible_fnsii_partners'] == []
        assert row['partner_status'] == 'unverified'
        assert row['model_eligible'] is False
        reaction = next(c['reaction'] for c in row['annotation']['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY')
        assert reaction['ecNumber'] == '1.6.2.4'
        ids = {r['id'] for r in reaction['reactionCrossReferences']}
        assert review['cpr_annotation_master'] in ids
        assert set(review['cpr_annotation_carriers']) <= ids
        assert {e['evidenceCode'] for e in reaction['evidences']} == {'ECO:0000256'}
    fnsii = json.loads(Path(review['fnsii_review']).read_text())
    assert set(review['fnsii_annotation_carriers']) == {fnsii['reduced_carrier']['participant'], fnsii['oxidized_carrier']['participant']}
    assert not set(review['cpr_annotation_carriers']) & set(review['fnsii_annotation_carriers'])
    assert report['model_eligible'] is review['model_eligible'] is False
