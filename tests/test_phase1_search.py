from cannabis_carbon.phase1_search import fasta_accession


def test_archive_and_knowledgebase_accessions_keep_their_namespaces():
    assert fasta_accession('UPI003D48F166') == 'UPI003D48F166'
    assert fasta_accession('sp|Q09765|UBA3_SCHPO') == 'Q09765'
    assert fasta_accession('tr|A0A7J6DQR0|A0A7J6DQR0_CANSA') == 'A0A7J6DQR0'
