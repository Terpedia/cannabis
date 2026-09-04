import json

from cannabis_carbon.test_hypotheses import build_test_hypotheses


def test_test_hypotheses_include_assay_plan(tmp_path):
    queue = tmp_path / "queue.json"
    network = tmp_path / "network.json"
    output = tmp_path / "hypotheses.json"
    lineage = tmp_path / "lineage.json"
    compounds = tmp_path / "compounds.json"
    queue.write_text(json.dumps({"items": [{"id": "h1", "kind": "curated", "status": "candidate", "reaction_id": "r1", "candidate_proteins": [{"proteinId": "p1", "accession": "P1", "label": "candidate", "specialized_search": {"hits": [{"evalue": 1e-20, "bitscore": 100}]}}]}]}))
    network.write_text(json.dumps({"entities": [{"id": "r1", "type": "biochemical_reaction", "label": "A = B"}, {"id": "a", "type": "metabolite", "label": "A"}, {"id": "b", "type": "metabolite", "label": "B"}], "statements": [{"subjectId": "r1", "predicate": "has_reactant", "objectEntityId": "a"}, {"subjectId": "r1", "predicate": "has_product", "objectEntityId": "b"}]}))
    lineage.write_text(json.dumps({"targets": [{"cannabisdb_id": "CDB1", "status": "unresolved", "reason": "no-terpedia-identity", "carbon_atom_count": 2}]}))
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "Test cannabinoid", "aliases": [], "formula": "C2H4", "smiles": "CC", "source_url": "https://cannabisdatabase.ca"}]}))
    summary = build_test_hypotheses(queue, network, output, lineage, compounds)
    report = json.loads(output.read_text())
    assert summary["total"] == 1
    assert summary["candidate"] == 1
    assert summary["with_reaction"] == 1
    assert summary["with_candidate_proteins"] == 1
    assert report["hypotheses"][0]["proposed_tests"][0]["step"] == "recombinant_assay"
    assert summary["metabolite_targets_total"] == 1
    assert report["metabolite_target_hypotheses"][0]["proposed_tests"][2]["step"] == "13CO2_lineage_validation"
    assert report["metabolite_target_hypotheses"][0]["label"] == "Test cannabinoid"
    assert report["metabolite_target_hypotheses"][0]["unresolved_carbon_atoms"] == 2
    assert report["metabolite_target_hypotheses"][0]["carbon_lineage_blocker"] == "no-terpedia-identity"


def test_test_hypotheses_include_non_enzymatic_conversion_plan(tmp_path):
    queue = tmp_path / "queue.json"
    network = tmp_path / "network.json"
    output = tmp_path / "hypotheses.json"
    queue.write_text(json.dumps({"items": []}))
    network.write_text(json.dumps({"entities": [
        {"id": "r1", "type": "biochemical_reaction", "label": "acid = neutral + CO2", "url": "https://example.test", "attributes": {"reactionClass": "non-enzymatic-decarboxylation", "conditions": "heat", "reactionSmiles": "CC(=O)O>>C=O.O=C=O"}},
        {"id": "a", "type": "metabolite", "label": "acid"}, {"id": "b", "type": "metabolite", "label": "neutral"}, {"id": "co2", "type": "metabolite", "label": "CO2"},
    ], "statements": [
        {"subjectId": "r1", "predicate": "has_reactant", "objectEntityId": "a"}, {"subjectId": "r1", "predicate": "has_product", "objectEntityId": "b"}, {"subjectId": "r1", "predicate": "has_product", "objectEntityId": "co2"},
    ]}))
    summary = build_test_hypotheses(queue, network, output)
    report = json.loads(output.read_text())
    assert summary["non_enzymatic_conversion_hypotheses"] == 1
    assert report["hypotheses"][0]["proposed_tests"][0]["step"] == "controlled_decarboxylation"
