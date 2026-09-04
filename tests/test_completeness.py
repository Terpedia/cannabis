import json
import gzip

from cannabis_carbon.completeness import compute_completeness


def test_completeness_separates_missing_metabolites_and_enzymes(tmp_path):
    network = {"entities": [
        {"id": "m:a", "type": "metabolite"}, {"id": "m:b", "type": "metabolite"},
        {"id": "m:c", "type": "metabolite"}, {"id": "r:1", "type": "biochemical_reaction"}
    ], "statements": [
        {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "m:b"}
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"carbon_atom_count": 1}]}))
    result = compute_completeness(network_path, compounds_path)
    assert result["terpedia"]["metabolites_without_reactions"] == 1
    assert result["terpedia"]["reactions_without_enzyme_association"] == 1


def test_has_catalytic_activity_counts_as_enzyme_association(tmp_path):
    network = {"entities": [{"id": "r:1", "type": "biochemical_reaction"}], "statements": [
        {"subjectId": "p:1", "predicate": "has_catalytic_activity", "objectEntityId": "r:1"}
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    result = compute_completeness(network_path, compounds_path)
    assert result["terpedia"]["reactions_with_enzyme_association"] == 1
    assert result["terpedia"]["reactions_without_enzyme_association"] == 0


def test_completeness_can_include_co2_lineage(tmp_path):
    network = {"entities": [], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    lineage_path = tmp_path / "lineage.json"
    lineage_path.write_text(json.dumps({"target_summary": {"supported": 1}, "reachable_carbon_nodes": 2, "resolved_carbon_edges": 1, "inferred_carbon_edges": 1, "candidate_carbon_edges": 0, "external_carbon_input_entity_count": 3, "carbon_source_policy": "CO2 only"}))
    result = compute_completeness(network_path, compounds_path, lineage_path=lineage_path)
    assert result["coverage"]["co2_lineage"]["target_summary"]["supported"] == 1


def test_completeness_splits_unannotated_reactions_by_candidate_proteins(tmp_path):
    network = {"entities": [
        {"id": "m:a", "type": "metabolite"}, {"id": "m:b", "type": "metabolite"},
        {"id": "r:1", "type": "biochemical_reaction"},
        {"id": "r:2", "type": "biochemical_reaction"},
    ], "statements": [
        {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "m:b"},
        {"subjectId": "r:2", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:2", "predicate": "has_product", "objectEntityId": "m:b"},
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    hypotheses_path = tmp_path / "hypotheses.json"
    hypotheses_path.write_text(json.dumps({"items": [{"reaction_id": "r:1", "candidate_proteins": [{"id": "p:1"}]}]}))
    result = compute_completeness(network_path, compounds_path, hypotheses_path=hypotheses_path)
    assert result["terpedia"]["reactions_without_enzyme_association"] == 2
    assert result["terpedia"]["reactions_without_enzyme_with_candidate_proteins"] == 1
    assert result["terpedia"]["reactions_without_enzyme_without_candidate_proteins"] == 1
