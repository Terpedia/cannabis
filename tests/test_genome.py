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


def test_genome_search_records_diamond_hit_and_ec(tmp_path):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(">sp|P1|TEST\nMPEPTIDE\n")
    queue = tmp_path / "queue.json"
    queue.write_text(json.dumps({"items": [{"candidate_proteins": [{"proteinId": "uniprot:p1", "accession": "P1"}]}]}))
    hits = tmp_path / "hits.tsv"
    hits.write_text("tr|P1|TEST\tsp|R1|REF\t65.0\t8\t8\t8\t1e-20\t40.0\t100.0\t100.0\n")
    refs = tmp_path / "refs.tsv"
    refs.write_text("Entry\tEC number\tProtein names\tOrganism\nR1\t1.1.1.1\tref\tplant\n")
    result = build_genome_search(queue, fasta, tmp_path / "report.json", hits, refs)
    assert result["candidate_protein_count"] == 1
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["candidate_proteins"][0]["diamond_search"]["strong_candidate_hit"]
    assert report["candidate_proteins"][0]["diamond_search"]["hits"][0]["reference_ec_numbers"] == ["1.1.1.1"]


def test_genome_search_preserves_specialized_search_when_general_hits_are_absent(tmp_path):
    fasta = tmp_path / "proteome.fasta"
    fasta.write_text(">sp|P1|TEST\nMPEPTIDE\n")
    queue = tmp_path / "queue.json"
    queue.write_text(json.dumps({"items": [{"candidate_proteins": [{"proteinId": "uniprot:p1", "accession": "P1", "specialized_search": {"method": "DIAMOND blastp", "strong_candidate_hit": True}}]}]}))
    report_path = tmp_path / "report.json"
    build_genome_search(queue, fasta, report_path)
    report = json.loads(report_path.read_text())
    assert report["candidate_proteins"][0]["diamond_search"]["strong_candidate_hit"]
