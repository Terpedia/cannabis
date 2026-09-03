import json

from cannabis_carbon.genome import build_genome_search


def test_genome_search_records_sequence_membership(tmp_path):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(">sp|P1|TEST\nMPEPTIDE\n")
    queue = tmp_path / "queue.json"
    queue.write_text(json.dumps({"items": [{"candidate_proteins": [{"proteinId": "uniprot:p1", "accession": "P1"}]}]}))
    result = build_genome_search(queue, fasta, tmp_path / "report.json")
    assert result["candidate_proteins_with_sequence"] == 1
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["candidate_proteins"][0]["sequence_search"]["length"] == 8
