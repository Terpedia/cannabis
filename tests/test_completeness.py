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
