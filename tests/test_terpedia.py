import gzip
import json

from cannabis_carbon.terpedia import cytoscape_elements, load_network


def test_terpedia_reaction_becomes_compound_edge():
    network = {"entities": [
        {"id": "chebi:a", "type": "metabolite", "label": "A"},
        {"id": "chebi:b", "type": "metabolite", "label": "B"},
        {"id": "rhea:1", "type": "biochemical_reaction", "label": "A = B", "url": "https://rhea.example/1"},
        {"id": "protein:p", "type": "protein", "label": "enzyme P"},
    ], "statements": [
        {"subjectId": "rhea:1", "predicate": "has_reactant", "objectEntityId": "chebi:a"},
        {"subjectId": "rhea:1", "predicate": "has_product", "objectEntityId": "chebi:b"},
        {"subjectId": "protein:p", "predicate": "catalyzes", "objectEntityId": "rhea:1", "qualifiers": {"directExperimentalEvidence": False}},
    ]}
    graph = cytoscape_elements(network)
    assert graph["stats"] == {"metabolites": 2, "reaction_edges": 1, "reactions": 1, "direction_overrides": 0}
    assert graph["edges"][0]["data"]["status"] == "candidate"


def test_cytoscape_applies_curated_reverse_master_direction():
    network = {"entities": [
        {"id": "chebi:a", "type": "metabolite", "label": "A"},
        {"id": "chebi:b", "type": "metabolite", "label": "B"},
        {"id": "rhea:1", "type": "biochemical_reaction", "label": "B = A", "url": "https://rhea.example/1"},
    ], "statements": [
        {"subjectId": "rhea:1", "predicate": "has_reactant", "objectEntityId": "chebi:b"},
        {"subjectId": "rhea:1", "predicate": "has_product", "objectEntityId": "chebi:a"},
    ]}
    graph = cytoscape_elements(network, {"rhea:1": {"orientation": "reverse_master", "directional_rhea_id": "2", "source": "https://rhea.example/2"}})
    edge = graph["edges"][0]["data"]
    assert (edge["source"], edge["target"]) == ("chebi:a", "chebi:b")
    assert edge["raw_direction"] == "product-to-reactant"
    assert edge["directional_rhea_id"] == "2"


def test_load_network_merges_source_linked_reaction_additions(tmp_path):
    source = tmp_path / "network.json.gz"
    with gzip.open(source, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    (tmp_path / "reaction-additions.json").write_text(json.dumps({"entities": [{"id": "rhea:add", "type": "biochemical_reaction"}], "statements": [{"subjectId": "rhea:add", "predicate": "has_product", "objectEntityId": "m:add"}]}))
    result = load_network(source)
    assert any(entity["id"] == "rhea:add" for entity in result["entities"])
    assert result["statements"][0]["subjectId"] == "rhea:add"
