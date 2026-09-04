import json

from cannabis_carbon.test_hypotheses import build_test_hypotheses


def test_test_hypotheses_include_assay_plan(tmp_path):
    queue = tmp_path / "queue.json"
    network = tmp_path / "network.json"
    output = tmp_path / "hypotheses.json"
    queue.write_text(json.dumps({"items": [{"id": "h1", "kind": "curated", "status": "candidate", "reaction_id": "r1", "candidate_proteins": [{"proteinId": "p1", "accession": "P1", "label": "candidate", "specialized_search": {"hits": [{"evalue": 1e-20, "bitscore": 100}]}}]}]}))
    network.write_text(json.dumps({"entities": [{"id": "r1", "type": "biochemical_reaction", "label": "A = B"}, {"id": "a", "type": "metabolite", "label": "A"}, {"id": "b", "type": "metabolite", "label": "B"}], "statements": [{"subjectId": "r1", "predicate": "has_reactant", "objectEntityId": "a"}, {"subjectId": "r1", "predicate": "has_product", "objectEntityId": "b"}]}))
    summary = build_test_hypotheses(queue, network, output)
    report = json.loads(output.read_text())
    assert summary == {"total": 1, "candidate": 1, "blocked": 0, "with_reaction": 1, "with_candidate_proteins": 1}
    assert report["hypotheses"][0]["proposed_tests"][0]["step"] == "recombinant_assay"
