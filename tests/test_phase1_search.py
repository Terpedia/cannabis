import pytest
import json
from pathlib import Path
from cannabis_carbon.phase1_search import fasta_accession, validate_single_reference, annotated_cds_reference


def test_archive_and_knowledgebase_accessions_keep_their_namespaces():
    assert fasta_accession('UPI003D48F166') == 'UPI003D48F166'
    assert fasta_accession('sp|Q09765|UBA3_SCHPO') == 'Q09765'
    assert fasta_accession('tr|A0A7J6DQR0|A0A7J6DQR0_CANSA') == 'A0A7J6DQR0'
    assert fasta_accession('WP_169336908.1') == 'WP_169336908.1'


def test_genbank_reference_requires_exact_version_and_single_sequence():
    assert validate_single_reference('WP_169336908.1', b'>WP_169336908.1 enzyme [organism]\nMALW\n') == 'WP_169336908.1 enzyme [organism]'
    for data in [b'>WP_169336908.2 enzyme\nMALW', b'>WP_169336908.1 enzyme\n',
                 b'>WP_169336908.1 enzyme\nMALW\n>other enzyme\nMALW', b'<html>error</html>']:
        with pytest.raises(ValueError):
            validate_single_reference('WP_169336908.1', data)


def test_nucleotide_translation_keeps_protein_and_nucleotide_namespaces_separate():
    feature = '''<GBFeature><GBFeature_key>CDS</GBFeature_key><GBFeature_location>1..12</GBFeature_location>
      <GBFeature_quals>
      <GBQualifier><GBQualifier_name>protein_id</GBQualifier_name><GBQualifier_value>QTEST123.1</GBQualifier_value></GBQualifier>
      <GBQualifier><GBQualifier_name>translation</GBQualifier_name><GBQualifier_value>MALW</GBQualifier_value></GBQualifier>
      </GBFeature_quals></GBFeature>'''
    def xml(features):
        return f'''<GBSet><GBSeq><GBSeq_primary-accession>MK123456</GBSeq_primary-accession>
          <GBSeq_accession-version>MK123456.1</GBSeq_accession-version><GBSeq_organism>Plant species</GBSeq_organism>
          <GBSeq_feature-table>{features}</GBSeq_feature-table></GBSeq></GBSet>'''.encode()
    fasta, metadata = annotated_cds_reference('MK123456', xml(feature))
    assert fasta.startswith(b'>QTEST123.1 ')
    assert metadata['nucleotide_accession_version'] == 'MK123456.1'
    assert metadata['protein_accession'] == 'QTEST123.1'
    assert metadata['cds_location'] == '1..12'
    with pytest.raises(ValueError, match='match'):
        annotated_cds_reference('MK123456.2', xml(feature))
    with pytest.raises(ValueError, match='ambiguous'):
        annotated_cds_reference('MK123456', xml(feature + feature))
    with pytest.raises(ValueError, match='translation'):
        annotated_cds_reference('MK123456', xml(feature.replace('translation', 'unused')))


def test_published_reference_resolution_and_hit_provenance_are_consistent():
    path = Path(__file__).resolve().parents[1] / 'docs/data/phase1-targeted-protein-search.json'
    report = json.loads(path.read_text())
    for row in report['rows']:
        source_ids = set(row['source_uniprot_ids']) | set(row['source_genbank_ids'])
        resolution = row['source_reference_resolution']
        assert set(resolution) == source_ids
        assert set(row['reference_sequences_present']) <= set(resolution.values())
        assert all(hit['reference_accession'] in row['reference_sequences_present'] for hit in row['sequence_hits'])
        assert row['screened_candidate_count'] == len({hit['cannabis_accession'] for hit in row['sequence_hits'] if hit['passes_screen']})
    for record in report['genbank_retrievals']:
        if record['status'] == 'retrieved-cds':
            assert record['protein_accession'] != record['nucleotide_accession_version']
            assert len(record['source_xml_sha256']) == 64
            assert record['cds_location']
